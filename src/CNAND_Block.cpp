#include "CNAND_Block.h"

#include <memory>
#include <random>
#include <string.h>
#include <assert.h>

#include "CVthDistribution_Normal.h"
#include "Logger.h"

#define ERASED_VOLTAGE -1.5
#define ERASED_VALUE 0

using namespace std;

CNAND_Block::CNAND_Block(CellType cell_type, BYTE* gray_codes, int wordlines, int cells_per_wordline,
                         CVthDistribution_base *vth_dist, double *thresh_volts) {
    m_cell_type = cell_type;
    m_wordlines_per_block = wordlines;
    m_cells_per_wordline = cells_per_wordline;
    m_total_cells = m_wordlines_per_block * m_cells_per_wordline;
    m_values_per_cell = (1 << cell_type);

    m_vth_dist = vth_dist;
    m_gray_codes = gray_codes;
    m_thresh_volts = thresh_volts;

    cell_values = new BYTE*[wordlines];
    cell_voltages = new double*[wordlines];

    for (int i = 0; i < wordlines; ++i) {
        cell_values[i] = new BYTE[cells_per_wordline];
        cell_voltages[i] = new double[cells_per_wordline];
    }

    cell_voltage_table = new double**[wordlines];
    for (int i = 0; i < wordlines; ++i) {
        cell_voltage_table[i] = new double*[cells_per_wordline];
        for (int j = 0; j < cells_per_wordline; ++j) {
            cell_voltage_table[i][j] = new double[m_values_per_cell];
        }
    }
}

CNAND_Block::~CNAND_Block() {
    for (int i = 0; i < m_wordlines_per_block; ++i) {
        delete[] cell_values[i];
        delete[] cell_voltages[i];
    }
    delete[] cell_values;
    delete[] cell_voltages;

    for (int i = 0; i < m_wordlines_per_block; ++i) {
        for (int j = 0; j < m_cells_per_wordline; ++j) {
            delete[] cell_voltage_table[i][j];
        }
        delete[] cell_voltage_table[i];
    }
    delete[] cell_voltage_table;
}

void CNAND_Block::Init(CVthDistribution_base *vth_dist, BYTE *gray_codes) {
    for (int i = 0; i < m_wordlines_per_block; i++) {
        BYTE *pvalue = cell_values[i];
        double *pvoltage = cell_voltages[i];

        for (int j = 0; j < m_cells_per_wordline; ++j) {
            pvalue[j] = static_cast<BYTE>(~0);
            pvoltage[j] = ERASED_VOLTAGE;

            double *pcell_voltage_table = cell_voltage_table[i][j];
            for (int k = 0; k < m_values_per_cell; ++k) {
                double v = vth_dist->GenerateVoltage(static_cast<BYTE>(k));
                pcell_voltage_table[k] = v;
            }
        }
    }
}

bool CNAND_Block::Erase() {
    for (int i = 0; i < m_wordlines_per_block; ++i) {
        double *pvoltage = cell_voltages[i];
        for (int j = 0; j < m_cells_per_wordline; ++j) {
            pvoltage[j] = ERASED_VOLTAGE;
        }
    }
    return true;
}

bool CNAND_Block::ProgramOneShot(int wordline, BYTE* data_in) {
    BYTE *pvalue = cell_values[wordline];
    memcpy(pvalue, data_in, m_cells_per_wordline);

    double *pvoltage = cell_voltages[wordline];
    for (int i = 0; i < m_cells_per_wordline; ++i) {
        pvoltage[i] = cell_voltage_table[wordline][i][pvalue[i]];
        LogDebug(2, "ProgramOneShot: wl=%d cell=%d value=%02x voltage=%6.3f", wordline, i, pvalue[i], pvoltage[i]);
    }

    return true;
}

void CNAND_Block::ReadPage(int wordline, int bit_idx, BYTE *buffer, double *thresh_volts, BYTE *gray_codes) {
    double *pwl = cell_voltages[wordline];
    for (int i = 0; i < m_cells_per_wordline; i++) {
        buffer[i] = Voltage2Bit(pwl[i], bit_idx, thresh_volts, gray_codes);
    }

    for (int i = 0; i < m_cells_per_wordline; ++i) {
        LogInfo("%1d", buffer[i]);
    }
}

void CNAND_Block::Evolve(double &cell_voltage, double v0, double t, double t0, double a,
                         double C, double tau, double Taf, double Eaf, double f) {
    cell_voltage = -a * log(C - (C - exp(-(v0 - f) / a)) * exp((-(t - t0) / tau) * Taf * Eaf)) + f;
}

void CNAND_Block::DbgDumpVoltage() {
    for (int wordline = 0; wordline < m_wordlines_per_block; wordline++) {
        for (int i = 0; i < m_cells_per_wordline; ++i) {
            double v = cell_voltages[wordline][i];
            LogInfo("wordline (%d)--cell (%d) voltage=%f, bin=%d", wordline, i, v, (int)((v + 1) / 0.025));
        }
    }
}

void CNAND_Block::DbgDumpVoltage(int wordline) {
    for (int i = 0; i < m_cells_per_wordline; ++i) {
        double v = cell_voltages[wordline][i];
        LogInfo("wordline (%d)--cell(%d) voltage=%f, bin=%d", wordline, i, v, (int)((v + 1) / 0.025));
    }
}

bool CNAND_Block::Voltage2Bit(double volt, int bit_idx, double *vth, BYTE *values) {
    // Check upper bound
    if (volt > vth[m_values_per_cell - 2]) {
        BYTE read_value = m_gray_codes[m_values_per_cell - 1];
        return GET_VALUE_BIT(read_value, bit_idx);
    }

    for (int i = 0; i < m_values_per_cell - 1; ++i) {
        if (volt < vth[i]) {
            BYTE read_value = m_gray_codes[i];
            return GET_VALUE_BIT(read_value, bit_idx);
        }
    }

    fprintf(stderr, "Voltage2Bit(): ERROR!\n");
#ifdef _WIN32
    DebugBreak();
#endif

    // Fallback for specific cell types (from original code comments)
    switch (m_cell_type) {
        case SLC_CELL:
            return volt < vth[0];
        case MLC_CELL:
            if (bit_idx == 0) {
                return (volt < vth[0] || volt > vth[2]);
            }
            if (bit_idx == 1) {
                return (volt < vth[1]);
            }
            break;
        case TLC_CELL:
        case QLC_CELL:
        default:
            break;
    }

    return false;
}

