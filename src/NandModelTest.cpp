#include <windows.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "SysTypes.h"
#include "CNAND_Die.h"
#include "CVthDistribution_NormalLaplaceMixture.h"
#include "CVthDistribution_Normal.h"
#include "CConfiguration.h"
#include "Logger.h"
#include "rand.h"

static int cell_bits = 0;
static int states = 0;
static int total_dies = 0;
static int blocks = 0;
static int wordlines_per_block = 0;
static int cells_per_wordline = 0;
static int num_threads = 0;

static int cdf_mode = 0;
static int RTN_mode = 0;
static double RTN_sigma = 0.0;

static double time_param = 0.0;
static double t0 = 0.0;
static double A = 0.0;
static double C = 0.0;
static double tau = 0.0;
static double Taf = 0.0;
static double Eaf = 0.0;
static double f_param = 0.0;

static std::vector<double> standard_deviations_vector;
static std::vector<double> means_vector;
static std::vector<double> thresh_volts_vector;
static std::vector<int> state_values_vector;
static std::vector<int> vth_range_vector;
static std::vector<double> vth_step_vector;

static std::vector<double> c0_lambda;
static std::vector<double> c1_lambda;
static std::vector<double> c2_lambda;
static std::vector<double> c0_mu;
static std::vector<double> c1_mu;
static std::vector<double> c2_mu;
static std::vector<double> c0_sigma;
static std::vector<double> c0_alpha;
static std::vector<double> c1_alpha;
static std::vector<double> c0_beta;
static std::vector<double> c1_beta;

static BYTE *state_values = nullptr;
static double *means = nullptr;
static double *standard_deviations = nullptr;
static double *default_vth = nullptr;
static int *vth_range = nullptr;
static double *vth_step = nullptr;

static int total_bins = 0;
static double *bin_voltages = nullptr;
static unsigned int **bin_per_value = nullptr;
static double **bin_area_per_value = nullptr;
static unsigned int **bin_uncovered = nullptr;

static unsigned int **transition_matrix = nullptr;
static double *bit_error_rate = nullptr;
static unsigned int *bits_counter = nullptr;
static double *page_read_one_ratio = nullptr;

static double **total_one_count_ratio_per_page = nullptr;
static double **total_one_count_ratio_square_per_page = nullptr;
static int **bit_asymmetry_per_page = nullptr;
static int **data_in_bit_asymmetry = nullptr;
static int **data_in_bva_square_per_page = nullptr;

static HANDLE bin_lock = nullptr;

static void ReadWordLine(CNAND_Die* die, int block, int wordline, BYTE* buffer, BYTE* tmp_buf, int num_bits, int cells_per_wordline) {
    // read MSB first
    die->ReadPage(block, wordline, num_bits - 1, buffer);
    for (int bit_idx = num_bits - 2; bit_idx >= 0; bit_idx--) {
        die->ReadPage(block, wordline, bit_idx, tmp_buf);
        for (int i = 0; i < cells_per_wordline; ++i) {
            buffer[i] = static_cast<BYTE>((buffer[i] << 1) | tmp_buf[i]);
        }
    }
}

static void PrintCDF(unsigned int** bin_per_value, int states, int total_bins, double *bin_voltage, double norm) {
    for (int s = 0; s < states; s++) {
        int v = state_values[s];
        double cdf_sum = 0;

        // 1. print vth-range, uncovered
        for (int vth_idx = 0, bin_idx = 1; vth_idx < states - 1; ++vth_idx) {
            if (bin_uncovered[v][state_values[vth_idx]] != 0) {
                int range;
                double left_volt;
                double right_volt;
                double volt_diff;

                range = 1;
                left_volt = default_vth[vth_idx] - vth_range[vth_idx] * vth_step[vth_idx];

                if (vth_idx > 0) {
                    left_volt = default_vth[vth_idx - 1] + (vth_range[vth_idx - 1] + 1) * vth_step[vth_idx - 1];
                }

                right_volt = default_vth[vth_idx] - vth_range[vth_idx] * vth_step[vth_idx];
                volt_diff = right_volt - left_volt;
                range = static_cast<int>(lround(volt_diff / vth_step[vth_idx])) + 1;

                for (int i = 0; i < range; ++i) {
                    double bin_cdf = 1.0 / range * norm * bin_uncovered[v][state_values[vth_idx]];
                    cdf_sum += bin_cdf;
                    printf("%.4f %.4g\n", left_volt + vth_step[vth_idx] * i, bin_cdf);
                }
            }

            // 1.2 print vth-range
            for (int step_idx = -vth_range[vth_idx]; step_idx < vth_range[vth_idx]; ++step_idx, ++bin_idx) {
                if (bin_per_value[v][bin_idx] != 0) {
                    double bin_cdf = norm * bin_per_value[v][bin_idx];
                    cdf_sum += bin_cdf;
                    printf("%.4f %.4g\n", bin_voltage[bin_idx], bin_cdf);
                }
            }
        }

        // 2. print the right boundary
        if (bin_uncovered[v][state_values[states - 1]] != 0) {
            double bin_cdf = norm * bin_uncovered[v][state_values[states - 1]];
            cdf_sum += bin_cdf;
            printf("%.4f %.4g\n", default_vth[states - 2] + vth_range[states - 2] * vth_step[states - 2], bin_cdf);
        }

        printf("\n");

        if (fabs(1.0 - cdf_sum) > 1E-9) {
            LogError("CDF_SUM != 1 diff=%.9g", 1.0 - cdf_sum);
        }
    }
}

static void PrintPDF(double** bin_area_per_value, int states, int total_bins, double *bin_voltage, double norm) {
    for (int s = 0; s < states; s++) {
        int v = state_values[s];
        for (int i = 1; i < total_bins; i++) {
            if (bin_area_per_value[v][i] > 0) {
                printf("%.4f %.4g\n", (bin_voltage[i - 1] + bin_voltage[i]) * 0.5, norm * bin_area_per_value[v][i]);
            }
        }
        printf("\n");
    }
}

static void GetTransitionMatrix(CNAND_Die* die, int block, int wordline, unsigned int** transition_matrix, BYTE* data_in) {
    BYTE* buffer = new BYTE[cells_per_wordline];
    BYTE* tmp_buffer = new BYTE[cells_per_wordline];

    for (int i = 0; i < cells_per_wordline; ++i) {
        buffer[i] = 0;
        tmp_buffer[i] = 0;
    }

    ReadWordLine(die, block, wordline, buffer, tmp_buffer, cell_bits, cells_per_wordline);
    for (int i = 0; i < cells_per_wordline; ++i) {
        InterlockedIncrement(reinterpret_cast<volatile LONG*>(&transition_matrix[data_in[i]][buffer[i]]));
    }

    delete[] buffer;
    delete[] tmp_buffer;
}

static void CalculateBER(unsigned int** transition_matrix, double* bit_error_rate, double norm) {
    unsigned int* error_bits_counter = new unsigned int[cell_bits];
    for (int i = 0; i < cell_bits; ++i) {
        error_bits_counter[i] = 0;
    }

    int num_states = 1 << cell_bits;
    for (int i = 0; i < num_states; ++i) {
        for (int j = 0; j < num_states; ++j) {
            int result = i ^ j;
            unsigned int count = transition_matrix[i][j];
            for (int bit_idx = 0; bit_idx < cell_bits; ++bit_idx) {
                if (GET_VALUE_BIT(result, bit_idx)) {
                    error_bits_counter[bit_idx] += count;
                }
            }
        }
    }

    for (int i = 0; i < cell_bits; ++i) {
        bit_error_rate[i] = error_bits_counter[i] * norm;
    }

    delete[] error_bits_counter;
}

static void PrintBER(unsigned int** transition_matrix, double* bit_error_rate, double norm_factor) {
    for (int i = 0; i < states; ++i) {
        printf("\t%d", state_values[i]);
    }
    printf("\n");

    for (int i = 0; i < states; ++i) {
        printf("%d", state_values[i]);
        for (int j = 0; j < states; ++j) {
            printf("\t%.4e", transition_matrix[state_values[i]][state_values[j]] * norm_factor);
        }
        printf("\n");
    }

    for (int i = 0; i < cell_bits; ++i) {
        printf("#BER[%d]:%.4e\n", i, bit_error_rate[i]);
    }
}

static void CountBitOne(CNAND_Die* die, int block, int wordline, unsigned int* bits_counter) {
    BYTE* buf = new BYTE[cells_per_wordline];
    for (int i = 0; i < cells_per_wordline; ++i) {
        buf[i] = 0;
    }

    int id = die->Id();

    for (int bit_idx = 0; bit_idx < cell_bits; ++bit_idx) {
        die->ReadPage(block, wordline, bit_idx, buf);
        unsigned int count = 0;
        for (int i = 0; i < cells_per_wordline; ++i) {
            count += buf[i];
        }
        double ratio = static_cast<double>(count) / cells_per_wordline;
        int bit_asymmetry = 2 * static_cast<int>(count) - cells_per_wordline;

        assert(bit_idx < cell_bits);
        assert(id < num_threads);

        bit_asymmetry_per_page[bit_idx][id] += bit_asymmetry;
        total_one_count_ratio_per_page[bit_idx][id] += ratio;
        total_one_count_ratio_square_per_page[bit_idx][id] += ratio * ratio;

        (void)bits_counter;
    }

    delete[] buf;
}

static void CountDataInBVA(CNAND_Die* die, BYTE* data_in, int** data_in_bit_asymmetry, int** data_in_bva_square_per_page) {
    int id = die->Id();
    assert(id < num_threads);

    for (int bit_idx = 0; bit_idx < cell_bits; ++bit_idx) {
        int count_one = 0;
        for (int i = 0; i < cells_per_wordline; ++i) {
            count_one += GET_VALUE_BIT(data_in[i], bit_idx);
        }
        int page_bva = 2 * count_one - cells_per_wordline;
        int page_bva_square = page_bva * page_bva;

        data_in_bva_square_per_page[bit_idx][id] += page_bva_square;
        data_in_bit_asymmetry[bit_idx][id] += page_bva;
    }
}

static DWORD WINAPI ProcessDie(LPVOID arg) {
    CNAND_Die *die = (CNAND_Die*)arg;
    die->Init();

    rand_seed(die->Id());

    BYTE *data_in = new BYTE[cells_per_wordline];
    double *thresh_volts = new double[states - 1];
    BYTE *data_buf = new BYTE[2 * cells_per_wordline];

    for (int i = 0; i < 2 * cells_per_wordline; ++i) {
        data_buf[i] = static_cast<BYTE>(rand() % states);
    }

    BYTE *buffer = new BYTE[cells_per_wordline];
    BYTE *tmp_buffer = new BYTE[cells_per_wordline];
    BYTE *prev_buffer = new BYTE[cells_per_wordline];
    bool *cell_vth_found = new bool[cells_per_wordline];

    int *wordline_buf_offset = new int[blocks * wordlines_per_block];
    for (int i = 0; i < blocks; i++) {
        for (int j = 0; j < wordlines_per_block; j++) {
            wordline_buf_offset[i * wordlines_per_block + j] = static_cast<int>(rand() % cells_per_wordline);
        }
    }

    const BYTE erased_value = static_cast<BYTE>(state_values_vector[0]);
    (void)erased_value;

    for (int value = 0; value < states; ++value) {
        for (int i = 0; i < 2 * cells_per_wordline; ++i) {
            data_buf[i] = static_cast<BYTE>((data_buf[i] + 1) % states);
        }

        for (int block = 0; block < blocks; ++block) {
            LogDebug(1, "Processing value=%d block=%d", value, block);

            for (int wordline = 0; wordline < wordlines_per_block; ++wordline) {
                LogDebug(2, "Processing value=%d block=%d wordline=%d", value, block, wordline);

                int offset = wordline_buf_offset[block * wordlines_per_block + wordline];
                memcpy(data_in, data_buf + offset, cells_per_wordline);

                CountDataInBVA(die, data_in, data_in_bit_asymmetry, data_in_bva_square_per_page);

                for (int cell = 0; cell < cells_per_wordline; ++cell) {
                    cell_vth_found[cell] = false;
                }

                LogDebug(2, "Program block=%d wordline=%d ...", block, wordline);
                die->ProgramOneShot(block, wordline, data_in);

                // BER simulation (default vth)
                die->SetVth(default_vth);
                GetTransitionMatrix(die, block, wordline, transition_matrix, data_in);
                CountBitOne(die, block, wordline, bits_counter);

                // Set vth to leftmost bounds
                for (int i = 0; i < states - 1; i++) {
                    thresh_volts[i] = default_vth[i] - vth_step[i] * vth_range[i];
                }
                die->SetVth(thresh_volts);

                ReadWordLine(die, block, wordline, prev_buffer, tmp_buffer, cell_bits, cells_per_wordline);

                // Starting scan procedure
                for (int vth_idx = 0, bin_idx = 1; vth_idx < states - 1; ++vth_idx) {
                    double delta = 1.0 / vth_step[vth_idx];
                    for (int step_idx = -vth_range[vth_idx]; step_idx < vth_range[vth_idx]; ++step_idx, ++bin_idx) {
                        // update Vth
                        thresh_volts[vth_idx] += vth_step[vth_idx];
                        die->SetVth(thresh_volts);

                        ReadWordLine(die, block, wordline, buffer, tmp_buffer, cell_bits, cells_per_wordline);

                        // Check bit flip
                        for (int cell = 0; cell < cells_per_wordline; ++cell) {
                            if (cell_vth_found[cell]) {
                                continue;
                            }
                            if (prev_buffer[cell] == state_values[vth_idx + 1] && buffer[cell] == state_values[vth_idx]) {
                                cell_vth_found[cell] = true;
                                LogDebug(2, "Found cell voltage[(%d]=%6.3f", cell, thresh_volts[vth_idx]);
                                assert(bin_idx < total_bins);
                                InterlockedIncrement(reinterpret_cast<volatile LONG*>(&bin_per_value[data_in[cell]][bin_idx]));

                                // acquire lock
                                WaitForSingleObject(bin_lock, 0);
                                bin_area_per_value[data_in[cell]][bin_idx] += delta;
                                ReleaseMutex(bin_lock);
                            }
                        }

                        // Swap buffer and prev_buffer
                        BYTE *tmp = prev_buffer;
                        prev_buffer = buffer;
                        buffer = tmp;
                    }

                    // restore vth
                    thresh_volts[vth_idx] = default_vth[vth_idx];
                }

                // handle cells without flip, not covered by the scan range
                for (int cell = 0; cell < cells_per_wordline; ++cell) {
                    if (cell_vth_found[cell]) {
                        continue;
                    }
                    int v_in = data_in[cell];
                    int v_out = buffer[cell];
                    InterlockedIncrement(reinterpret_cast<volatile LONG*>(&bin_uncovered[v_in][v_out]));
                }

                LogDebug(2, "Erasing block %d ...", block);
                die->EraseBlock(block);
            }
        }
    }

    delete[] thresh_volts;
    delete[] data_in;
    delete[] buffer;
    delete[] tmp_buffer;
    delete[] prev_buffer;
    delete[] cell_vth_found;
    delete[] wordline_buf_offset;
    delete[] data_buf;

    delete die;
    return 0;
}

int main(int argc, char * argv[]) {
    UINT32 start_time = GetTickCount();

    if (argc != 2) {
        fprintf(stderr, "Usage: %s <config file>\n", argv[0]);
        exit(1);
    }

    int dbg_level = 1;
    LogInit(dbg_level, std::cerr);

    CConfiguration *model_config = new CConfiguration(argv[1]);

    model_config->GetParamValueDouble("time", time_param);
    model_config->GetParamValueDouble("t0", t0);
    model_config->GetParamValueDouble("A", A);
    model_config->GetParamValueDouble("C", C);
    model_config->GetParamValueDouble("tau", tau);
    model_config->GetParamValueDouble("Taf", Taf);
    model_config->GetParamValueDouble("Eaf", Eaf);
    model_config->GetParamValueDouble("f", f_param);

    model_config->GetParamValueInt("cdf_mode", cdf_mode);
    model_config->GetParamValueInt("RTN_mode", RTN_mode);
    model_config->GetParamValueInt("num_threads", num_threads);

    model_config->GetParamValueInt("DebugLevel", dbg_level);
    LogSetDbgLevel(dbg_level);

    model_config->GetParamValueDouble("RTN_sigma", RTN_sigma);
    model_config->GetParamValueInt("cell_type", cell_bits);
    states = 1 << cell_bits;
    CellType cell_type = static_cast<CellType>(cell_bits);

    model_config->GetParamValueDouble("PE_cycle_L", time_param);
    model_config->GetParamValueInt("total_dies", total_dies);
    model_config->GetParamValueInt("total_number_of_blocks", blocks);
    model_config->GetParamValueInt("number_of_wordlines_per_block", wordlines_per_block);
    model_config->GetParamValueInt("number_of_cells_per_wordline", cells_per_wordline);

    model_config->GetParamValueIntVec("state_values", state_values_vector, states);
    model_config->GetParamValueDoubleVec("threshold_voltages", thresh_volts_vector, states - 1);
    model_config->GetParamValueDoubleVec("means", means_vector, states);
    model_config->GetParamValueDoubleVec("standard_deviations", standard_deviations_vector, states);
    model_config->GetParamValueIntVec("vth_range", vth_range_vector, states - 1);
    model_config->GetParamValueDoubleVec("vth_step", vth_step_vector, states - 1);

    model_config->GetParamValueDoubleVec("c0_lambda", c0_lambda, states * states);
    model_config->GetParamValueDoubleVec("c1_lambda", c1_lambda, states * states);
    model_config->GetParamValueDoubleVec("c2_lambda", c2_lambda, states * states);

    model_config->GetParamValueDoubleVec("c0_mu", c0_mu, states);
    model_config->GetParamValueDoubleVec("c1_mu", c1_mu, states);
    model_config->GetParamValueDoubleVec("c2_mu", c2_mu, states);

    model_config->GetParamValueDoubleVec("c0_sigma", c0_sigma, states);
    model_config->GetParamValueDoubleVec("c0_alpha", c0_alpha, states);
    model_config->GetParamValueDoubleVec("c1_alpha", c1_alpha, states);
    model_config->GetParamValueDoubleVec("c0_beta", c0_beta, states);
    model_config->GetParamValueDoubleVec("c1_beta", c1_beta, states);

    state_values = new BYTE[states];
    for (int i = 0; i < states; ++i) {
        state_values[i] = static_cast<BYTE>(state_values_vector[i]);
    }

    vth_range = new int[states - 1];
    vth_step = new double[states - 1];
    for (int i = 0; i < states - 1; ++i) {
        vth_range[i] = vth_range_vector[i];
        vth_step[i] = vth_step_vector[i];
    }

    means = new double[states];
    for (int i = 0; i < states; ++i) {
        means[i] = means_vector[i];
    }

    standard_deviations = new double[states];
    for (int i = 0; i < states; ++i) {
        standard_deviations[i] = standard_deviations_vector[i];
    }

    transition_matrix = new unsigned int*[states];
    for (int i = 0; i < states; ++i) {
        transition_matrix[i] = new unsigned int[states];
        for (int j = 0; j < states; ++j) {
            transition_matrix[i][j] = 0;
        }
    }

    bit_error_rate = new double[cell_bits];
    for (int i = 0; i < cell_bits; ++i) {
        bit_error_rate[i] = 0;
    }

    int pages_per_block = cell_bits * wordlines_per_block;
    bits_counter = new unsigned int[pages_per_block];
    for (int i = 0; i < pages_per_block; ++i) {
        bits_counter[i] = 0;
    }

    page_read_one_ratio = new double[pages_per_block];
    for (int i = 0; i < pages_per_block; ++i) {
        page_read_one_ratio[i] = 0;
    }

    CVthDistribution_NormalLaplaceMixture *vth_dist = new CVthDistribution_NormalLaplaceMixture(cell_bits, state_values);
    vth_dist->Init(state_values, time_param, c0_sigma, c0_lambda, c1_lambda, c2_lambda,
                   c0_mu, c1_mu, c2_mu, c0_alpha, c1_alpha, c0_beta, c1_beta);

    BYTE *gray_codes = new BYTE[states];
    for (int i = 0; i < states; i++) {
        gray_codes[i] = state_values[i];
    }

    default_vth = new double[states - 1];
    for (int i = 0; i < states - 1; i++) {
        default_vth[i] = thresh_volts_vector[i];
        LogDebug(3, "vth[%d]=%6.3f", i, default_vth[i]);
    }

    total_bins = 1;
    for (int i = 0; i < states - 1; i++) {
        total_bins += 2 * vth_range[i];
    }

    bin_voltages = new double[total_bins];
    bin_voltages[0] = default_vth[0] - vth_range[0] * vth_step[0];

    for (int i = 0, bin_idx = 1; i < states - 1; i++) {
        int nbin = 2 * vth_range[i];
        double volt = default_vth[i] - (vth_range[i] - 1) * vth_step[i];
        for (int j = 0; j < nbin; j++) {
            bin_voltages[bin_idx++] = volt;
            volt += vth_step[i];
        }
    }

    for (int i = 0; i < total_bins; i++) {
        LogDebug(2, "bin_voltages[%d]=%.3f", bin_voltages[i]);
    }

    bin_per_value = new unsigned int*[states];
    bin_area_per_value = new double*[states];
    for (int value = 0; value < states; value++) {
        bin_per_value[value] = new unsigned int[total_bins];
        memset(bin_per_value[value], 0, total_bins * sizeof(unsigned int));

        bin_area_per_value[value] = new double[total_bins];
        memset(bin_area_per_value[value], 0, total_bins * sizeof(double));
    }

    bin_uncovered = new unsigned int*[states];
    for (int i = 0; i < states; ++i) {
        bin_uncovered[i] = new unsigned int[states];
        memset(bin_uncovered[i], 0, states * sizeof(unsigned int));
    }

    total_one_count_ratio_per_page = new double*[cell_bits];
    total_one_count_ratio_square_per_page = new double*[cell_bits];
    bit_asymmetry_per_page = new int*[cell_bits];
    data_in_bit_asymmetry = new int*[cell_bits];
    data_in_bva_square_per_page = new int*[cell_bits];

    for (int i = 0; i < cell_bits; i++) {
        total_one_count_ratio_per_page[i] = new double[num_threads];
        total_one_count_ratio_square_per_page[i] = new double[num_threads];
        bit_asymmetry_per_page[i] = new int[num_threads];
        data_in_bit_asymmetry[i] = new int[num_threads];
        data_in_bva_square_per_page[i] = new int[num_threads];

        for (int j = 0; j < num_threads; j++) {
            total_one_count_ratio_per_page[i][j] = 0;
            total_one_count_ratio_square_per_page[i][j] = 0;
            bit_asymmetry_per_page[i][j] = 0;
            data_in_bit_asymmetry[i][j] = 0;
            data_in_bva_square_per_page[i][j] = 0;
        }
    }

    HANDLE *thread_handles = new HANDLE[total_dies];
    bin_lock = CreateMutex(NULL, false, NULL);

    int die_iter = 0;
    int running_threads = 0;

    for (die_iter = 0; die_iter < num_threads && die_iter < total_dies; ++die_iter) {
        LogInfo("Processing die %d", die_iter);
        CNAND_Die *die = new CNAND_Die(die_iter, cell_type, state_values, blocks, wordlines_per_block,
                                      cells_per_wordline, vth_dist, default_vth, RTN_mode, RTN_sigma);
        thread_handles[die_iter] = CreateThread(NULL, 0, ProcessDie, die, 0, NULL);
        ++running_threads;
    }

    for (; die_iter < total_dies; die_iter++) {
        int idx = WaitForMultipleObjects(num_threads, thread_handles, FALSE, INFINITE) - WAIT_OBJECT_0;
        LogInfo("Processing die %d", die_iter);

        CNAND_Die *die = new CNAND_Die(idx, cell_type, state_values, blocks, wordlines_per_block,
                                      cells_per_wordline, vth_dist, default_vth, RTN_mode, RTN_sigma);
        thread_handles[idx] = CreateThread(NULL, 0, ProcessDie, die, 0, NULL);
    }

    WaitForMultipleObjects(running_threads, thread_handles, TRUE, INFINITE);

    double norm_factor = 1.0 / (total_dies * blocks * wordlines_per_block * cells_per_wordline);

    CalculateBER(transition_matrix, bit_error_rate, norm_factor / states);

    if (cdf_mode) {
        PrintCDF(bin_per_value, states, total_bins, bin_voltages, norm_factor);
    } else {
        PrintPDF(bin_area_per_value, states, total_bins, bin_voltages, norm_factor);
    }

    PrintBER(transition_matrix, bit_error_rate, norm_factor);

    // Compute Read One Ratio
    for (int i = 0; i < cell_bits; i++) {
        double total_one_count_ratio = 0;
        double total_one_count_ratio_square = 0;
        int total_bit_asymmetry_per_page = 0;
        int total_data_in_bit_asymmetry = 0;
        int total_data_in_bva_square = 0;

        for (int j = 0; j < num_threads; j++) {
            total_one_count_ratio += total_one_count_ratio_per_page[i][j];
            total_one_count_ratio_square += total_one_count_ratio_square_per_page[i][j];
            total_bit_asymmetry_per_page += bit_asymmetry_per_page[i][j];
            total_data_in_bit_asymmetry += data_in_bit_asymmetry[i][j];
            total_data_in_bva_square += data_in_bva_square_per_page[i][j];
        }

        double average_bit_asymmetry_per_page = ((double)total_bit_asymmetry_per_page) /
            (states * total_dies * blocks * wordlines_per_block);
        double average_one_count_ratio = total_one_count_ratio /
            (states * total_dies * blocks * wordlines_per_block);
        double average_one_count_ratio_square = total_one_count_ratio_square /
            (states * total_dies * blocks * wordlines_per_block);
        double read_one_ratio_sigma = sqrt(average_one_count_ratio_square -
                                           average_one_count_ratio * average_one_count_ratio);

        double average_data_in_bva = ((double)total_data_in_bit_asymmetry) /
            (states * total_dies * blocks * wordlines_per_block);
        double average_data_in_bva_square = ((double)total_data_in_bva_square) /
            (states * total_dies * blocks * wordlines_per_block);
        double data_in_bva_sigma = sqrt(average_data_in_bva_square - average_data_in_bva * average_data_in_bva);

        printf("# bit %d: pec=%.2f\t average_bit_asymmetry_per_page=%.4f\t average_one_count_ratio=%.4f\t read_one_ratio_sigma=%.6f\t average_data_in_bva=%.4f\t data_in_bva_sigma=%.4f\n",
               i, time_param, average_bit_asymmetry_per_page, average_one_count_ratio,
               read_one_ratio_sigma, average_data_in_bva, data_in_bva_sigma);
    }

    if (bin_lock) {
        CloseHandle(bin_lock);
    }

    UINT32 end_time = GetTickCount();
    LogInfo("Elapsed time: %d seconds", (end_time - start_time) / 1000);

    return 0;
}

