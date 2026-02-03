#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <random>
#include <assert.h>
#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#include <windows.h>
#endif

#include "CNelderMeadAlgorithm.h"
#include "CConfiguration.h"
#include "Logger.h"
#include "SysTypes.h"

struct NMParameters {
    int fitting_mode = 0;
    int dist_function_mode = 0;
    int trace_path = 0;
    int max_iterations = 0;
    double tolerance = 1e-6;
    double reflection_coeff = 1.0;
    double expansion_coeff = 2.0;
    double contraction_coeff = 0.5;
    std::vector<double> initial_params;
};

struct NLMParams {
    int num_states = 0;

    std::vector<std::vector<double>> lambda;
    std::vector<std::vector<int>> lambda_mask;
    double lambda_step = 0.0;

    std::vector<double> mu;
    std::vector<int> mu_mask;
    double mu_step = 0.0;

    std::vector<double> sigma;
    std::vector<int> sigma_mask;
    double sigma_step = 0.0;

    std::vector<double> alpha;
    std::vector<int> alpha_mask;
    double alpha_step = 0.0;

    std::vector<double> beta;
    std::vector<int> beta_mask;
    double beta_step = 0.0;
};

static std::vector<std::vector<double>> data_x; // bin index
static std::vector<std::vector<double>> data_y; // number of cells per bin

static NMParameters method_parameters;
static NLMParams normal_laplace_mix_params;

static void LoadDataFile(const char* filename, int num_states);
static void ParamGenerateSimplex(NLMParams& normal_laplace_mix_params, std::vector<std::vector<double>>& initial_simplex, int num_vertices);
static void ParamSerialize(const NLMParams& normal_laplace_mix_params, std::vector<double>& params);
static void ParamDeserialize(NLMParams& normal_laplace_mix_params, const std::vector<double>& params);
static void PrintResult(const NLMParams& normal_laplace_mix_params);
static double TargetFunction(std::vector<double>& params);

int main(int argc, char** argv) {
    LogInit(1, std::cout);
    /*Load parameters... .*/
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <config file>\n", argv[0]);
        exit(1);
    }

    std::vector<double> optimal;

    CConfiguration* NLM_parameter_config = new CConfiguration(argv[1]);
    NLM_parameter_config->GetParamValueInt("fitting_mode", method_parameters.fitting_mode);
    NLM_parameter_config->GetParamValueInt("dist_function_mode", method_parameters.dist_function_mode);
    NLM_parameter_config->GetParamValueInt("trace_path", method_parameters.trace_path);
    NLM_parameter_config->GetParamValueInt("max_iterations", method_parameters.max_iterations);

    NLM_parameter_config->GetParamValueInt("num_states", normal_laplace_mix_params.num_states);
    int num_states = normal_laplace_mix_params.num_states;

    NLM_parameter_config->GetParamValueDoubleMatrix("lambda", normal_laplace_mix_params.lambda, num_states);
    NLM_parameter_config->GetParamMatrixMask("lambda_mask", normal_laplace_mix_params.lambda_mask, num_states);
    NLM_parameter_config->GetParamValueDouble("lambda_step", normal_laplace_mix_params.lambda_step);

    NLM_parameter_config->GetParamValueDoubleVec("mu", normal_laplace_mix_params.mu, num_states);
    NLM_parameter_config->GetParamVectorMask("mu_mask", normal_laplace_mix_params.mu_mask, num_states);
    NLM_parameter_config->GetParamValueDouble("mu_step", normal_laplace_mix_params.mu_step);

    NLM_parameter_config->GetParamValueDoubleVec("sigma", normal_laplace_mix_params.sigma, num_states);
    NLM_parameter_config->GetParamVectorMask("sigma_mask", normal_laplace_mix_params.sigma_mask, num_states);
    NLM_parameter_config->GetParamValueDouble("sigma_step", normal_laplace_mix_params.sigma_step);

    NLM_parameter_config->GetParamValueDoubleVec("alpha", normal_laplace_mix_params.alpha, num_states);
    NLM_parameter_config->GetParamVectorMask("alpha_mask", normal_laplace_mix_params.alpha_mask, num_states);
    NLM_parameter_config->GetParamValueDouble("alpha_step", normal_laplace_mix_params.alpha_step);

    NLM_parameter_config->GetParamValueDoubleVec("beta", normal_laplace_mix_params.beta, num_states);
    NLM_parameter_config->GetParamVectorMask("beta_mask", normal_laplace_mix_params.beta_mask, num_states);
    NLM_parameter_config->GetParamValueDouble("beta_step", normal_laplace_mix_params.beta_step);

    LoadDataFile(argv[2], num_states);

    ParamSerialize(normal_laplace_mix_params, method_parameters.initial_params);
    int num_params = static_cast<int>(method_parameters.initial_params.size());

    // prepare initial simplex
    std::vector<std::vector<double>> initial_simplex;
    int num_vertices = num_params + 1;
    // int num_vertices = 10 * num_params;

    // curve fitting
    CNelderMeadAlgorithm* nelder_mead_method = new CNelderMeadAlgorithm(num_params, num_vertices, method_parameters.trace_path);

    // Main Loop
    double fmin_prev = 0;
    double fmin = 0;
    int round = 1;

#if 1
    do {
        ParamGenerateSimplex(normal_laplace_mix_params, initial_simplex, num_vertices);
        nelder_mead_method->Init(initial_simplex);

        if (method_parameters.fitting_mode == 0) {
            nelder_mead_method->Curvefitting0(TargetFunction, optimal, method_parameters.max_iterations, method_parameters.tolerance, num_params,
                                              method_parameters.reflection_coeff, method_parameters.expansion_coeff, method_parameters.contraction_coeff);
        } else if (method_parameters.fitting_mode == 1) {
            nelder_mead_method->Curvefitting1(TargetFunction, optimal, method_parameters.max_iterations, method_parameters.tolerance, num_params,
                                              method_parameters.reflection_coeff, method_parameters.expansion_coeff, method_parameters.contraction_coeff);
        } else if (method_parameters.fitting_mode == 2) {
            nelder_mead_method->Curvefitting2(TargetFunction, optimal, method_parameters.max_iterations, method_parameters.tolerance, num_params,
                                              method_parameters.reflection_coeff, method_parameters.expansion_coeff, method_parameters.contraction_coeff);
        } else {
            nelder_mead_method->Curvefitting2(TargetFunction, optimal, method_parameters.max_iterations, method_parameters.tolerance, num_params,
                                              method_parameters.reflection_coeff, method_parameters.expansion_coeff, method_parameters.contraction_coeff);
        }

        fmin_prev = fmin;
        fmin = TargetFunction(optimal);
        ParamDeserialize(normal_laplace_mix_params, optimal);
        LogInfo("================ Round: %d fmin=%.9g =================", round++, fmin);
        PrintResult(normal_laplace_mix_params);
    } while (1 && fmin_prev != fmin);
#endif

#if 0 // draw B_state (mu, beta) space
    double mu = 1.90;
    while (mu < 2.30) {
        double beta_inv = 0;
        while (beta_inv < 0.1) {
            std::vector<double> params;

            normal_laplace_mix_params.mu[2] = mu;
            normal_laplace_mix_params.beta[2] = 1.0 / beta_inv;

            ParamSerialize(normal_laplace_mix_params, params);
            double F = TargetFunction(params);

            //if (isnan(F) || isinf(F))
            printf("%.9g %.9g %.9g\n", mu, beta_inv, F);

            beta_inv += 0.0002;
        }

        mu += 0.0005;
    }
#endif

#if 0 // draw C_state (mu, beta) space
    double mu = 3.20;
    while (mu < 3.40) {
        double beta_inv = 0.001;
        while (beta_inv < 0.1) {
            std::vector<double> params;

            normal_laplace_mix_params.mu[3] = mu;
            normal_laplace_mix_params.beta[3] = 1 / beta_inv;

            ParamSerialize(normal_laplace_mix_params, params);
            double F = TargetFunction(params);

            printf("%.9g %.9g %.9g\n", mu, beta_inv, F);

            beta_inv += 0.0002;
        }

        mu += 0.0004;
    }
#endif

#if 0 // draw C_state (mu, sigma) space
    double mu = 3.28;
    while (mu < 3.36) {
        double sigma = 0;
        while (sigma < 0.01) {
            std::vector<double> params;

            normal_laplace_mix_params.mu[3] = mu;
            normal_laplace_mix_params.sigma[3] = sigma;

            ParamSerialize(normal_laplace_mix_params, params);
            double F = TargetFunction(params);

            printf("%.9g %.9g %.9g\n", mu, sigma, F);

            sigma += 0.0001;
        }

        mu += 0.0001;
    }
#endif

    delete nelder_mead_method;
    delete NLM_parameter_config;

    return 0;
}

static void LoadDataFile(const char* filename, int num_states) {
    char buf[1024];

    // Open file
    FILE* fp;
    fp = fopen(filename, "r");
    if (!fp) {
        std::cout << "Error opening file";
        exit(1);
    }

    int state_idx = 0;
    data_x.resize(num_states);
    data_y.resize(num_states);
    while (fgets(buf, 1024, fp)) {
        const char* xs;
        const char* ys;

        const char* delimiters = " \t\r\n";

        xs = buf;
        while ((*xs == ' ') || (*xs == '\t')) {
            xs++;
            continue;
        }

        if (*xs == '#' || *xs == '\r' || *xs == '\n' || *xs == 0) {
            continue;
        }
        if (*xs == '&') {
            state_idx++;
            continue;
        }

        xs = strtok(buf, delimiters);
        if (!xs) {
            LogError("Operation failed!");
            exit(-1);
        }

        ys = strtok(NULL, delimiters);
        double x = strtod(xs, NULL);
        double y = strtod(ys, NULL);
        // LogInfo("x=%s y=%s", xs, ys);
        data_x[state_idx].push_back(x);
        data_y[state_idx].push_back(y);
        assert(state_idx <= num_states);
    }
}

static void ParamGenerateSimplex(NLMParams& normal_laplace_mix_params, std::vector<std::vector<double>>& initial_simplex, int num_vertices) {
#if 0
    std::vector<double> initial_params;
    std::uniform_real_distribution<double> unidist(0, 1);
    std::default_random_engine rand_engine;

    ParamSerialize(normal_laplace_mix_params, initial_params);

    initial_simplex.resize(num_vertices);
    for (int p = 0; p < num_vertices; p++) {
        initial_simplex[p] = initial_params;
        int param_idx = 0;
        int num_states = normal_laplace_mix_params.num_states;

        // lambda
        for (int i = 0; i < num_states; i++) {
            for (int j = 0; j < num_states; j++) {
                if (normal_laplace_mix_params.lambda_mask[i][j]) {
                    initial_simplex[p][param_idx] = (1.0 + 0.02 * unidist(rand_engine));
                    ++param_idx;
                }
            }
        }

        // mu
        for (int i = 0; i < num_states; ++i) {
            if (normal_laplace_mix_params.mu_mask[i]) {
                initial_simplex[p][param_idx] += normal_laplace_mix_params.mu_step * 2 * unidist(rand_engine);
                ++param_idx;
            }
        }

        // sigma
        for (int i = 0; i < num_states; ++i) {
            if (normal_laplace_mix_params.sigma_mask[i] != 0) {
                initial_simplex[p][param_idx] += normal_laplace_mix_params.sigma_step * 2 * unidist(rand_engine);
                ++param_idx;
            }
        }

        // alpha
        for (int i = 0; i < num_states; ++i) {
            if (normal_laplace_mix_params.alpha_mask[i] != 0) {
                initial_simplex[p][param_idx] *= (1.0 + 0.02 * unidist(rand_engine)); // += normal_laplace_mix_params.alpha_inv_step;
                ++param_idx;
            }
        }

        // beta_inv
        for (int i = 0; i < num_states; ++i) {
            if (normal_laplace_mix_params.beta_mask[i] != 0) {
                initial_simplex[p][param_idx] *= (1.0 + 0.02 * unidist(rand_engine)); // += normal_laplace_mix_params.beta_inv_step;
                ++param_idx;
            }
        }
    }
#else
    std::vector<double> initial_params;
    ParamSerialize(normal_laplace_mix_params, initial_params);

    initial_simplex.resize(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        initial_simplex[i] = initial_params;
    }

    int param_idx = 0;
    int num_states = normal_laplace_mix_params.num_states;

    // lambda
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            if (normal_laplace_mix_params.lambda_mask[i][j]) {
                //initial_simplex[param_idx][param_idx] *= (1.0 + normal_laplace_mix_params.lambda_step);
                initial_simplex[param_idx][param_idx] += normal_laplace_mix_params.lambda_step;
                ++param_idx;
            }
        }
    }

    // mu
    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.mu_mask[i] != 0) {
            initial_simplex[param_idx][param_idx] += normal_laplace_mix_params.mu_step;
            ++param_idx;
        }
    }

    // sigma
    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.sigma_mask[i] != 0) {
            initial_simplex[param_idx][param_idx] += normal_laplace_mix_params.sigma_step;
            ++param_idx;
        }
    }

    // alpha_inv
    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.alpha_mask[i] != 0) {
            //initial_simplex[param_idx][param_idx] *= 1.001; // += normal_laplace_mix_params.alpha_inv_step;
            initial_simplex[param_idx][param_idx] += normal_laplace_mix_params.alpha_step;
            ++param_idx;
        }
    }

    // beta_inv
    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.beta_mask[i] != 0) {
            //initial_simplex[param_idx][param_idx] *= 1.001; // += normal_laplace_mix_params.beta_inv_step;
            initial_simplex[param_idx][param_idx] += normal_laplace_mix_params.beta_step;
            ++param_idx;
        }
    }
#endif
}

static void ParamSerialize(const NLMParams& normal_laplace_mix_params, std::vector<double>& params) {
    int num_states = normal_laplace_mix_params.num_states;

    for (int i = 0; i < num_states; ++i) {
        for (int j = 0; j < num_states; ++j) {
            if (normal_laplace_mix_params.lambda_mask[i][j] != 0) {
                params.push_back(normal_laplace_mix_params.lambda[i][j]);
            }
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.mu_mask[i] != 0) {
            params.push_back(normal_laplace_mix_params.mu[i]);
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.sigma_mask[i] != 0) {
            params.push_back(normal_laplace_mix_params.sigma[i]);
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.alpha_mask[i] != 0) {
            params.push_back(normal_laplace_mix_params.alpha[i]);
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.beta_mask[i] != 0) {
            params.push_back(normal_laplace_mix_params.beta[i]);
        }
    }
}

static void ParamDeserialize(NLMParams& normal_laplace_mix_params, const std::vector<double>& params) {
    int param_idx = 0;
    int num_states = normal_laplace_mix_params.num_states;

    for (int i = 0; i < num_states; ++i) {
        for (int j = 0; j < num_states; ++j) {
            if (normal_laplace_mix_params.lambda_mask[i][j] != 0) {
                normal_laplace_mix_params.lambda[i][j] = params[param_idx++];
            }
        }
    }

    // fix diagonal elements
    for (int i = 0; i < num_states; ++i) {
        normal_laplace_mix_params.lambda[i][i] = 1;
        for (int j = 0; j < num_states; ++j) {
            if (i != j) {
                normal_laplace_mix_params.lambda[i][i] -= normal_laplace_mix_params.lambda[i][j];
            }
        }
    }

    // !!!!
    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.mu_mask[i] != 0) {
            normal_laplace_mix_params.mu[i] = params[param_idx++];
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.sigma_mask[i] != 0) {
            normal_laplace_mix_params.sigma[i] = params[param_idx++];
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.alpha_mask[i] != 0) {
            normal_laplace_mix_params.alpha[i] = params[param_idx++];
        }
    }

    for (int i = 0; i < num_states; ++i) {
        if (normal_laplace_mix_params.beta_mask[i] != 0) {
            normal_laplace_mix_params.beta[i] = params[param_idx++];
        }
    }
}

static void PrintResult(const NLMParams& normal_laplace_mix_params) {
    // print number of states.
    int states = normal_laplace_mix_params.num_states;
    printf("\nstates: %d\n", states);

    // print lambda matrix.
    for (int i = 0; i < states; ++i) {
        printf("lambda[%d]", i);
        for (int j = 0; j < states; ++j) {
            printf("%.9g ", normal_laplace_mix_params.lambda[i][j]);
        }
        printf("\n");
    }

    // print mu vector.
    printf("mu: ");
    for (int i = 0; i < states; ++i) {
        printf("%.9g ", normal_laplace_mix_params.mu[i]);
    }
    printf("\n");

    // print sigma vector.
    printf("sigma: ");
    for (int i = 0; i < states; ++i) {
        printf("%.9g ", normal_laplace_mix_params.sigma[i]);
    }
    printf("\n");

    // print alpha
    printf("alpha: ");
    for (int i = 0; i < states; ++i) {
        printf("%.9g ", normal_laplace_mix_params.alpha[i]);
    }
    printf("\n");

    // print beta
    printf("beta: ");
    for (int i = 0; i < states; ++i) {
        printf("%.9g ", normal_laplace_mix_params.beta[i]);
    }
    printf("\n");
}

const double inv_sqrt2 = 1.0 / sqrt(2);
const double inv_sqrt_2pi = 1.0 / sqrt(2 * M_PI);
const double sqrt_2pi = sqrt(2 * M_PI);

extern "C" {
double erfcx(double x); // special case for real x
}

double PDF_Normal(long double x) {
    double y = exp(-x * x * 0.5) * inv_sqrt_2pi;
    return y;
}

double CDF_Normal(long double x) {
    // double y = (1 + erf(x / sqrt(2))) / 2;
    double y = 1.0 - 0.5 * erfc(x * inv_sqrt2);
    return y;
}

double CDF_Normal(long double x, long double y) {
    double cdf = -0.5 * (erfc(y * inv_sqrt2) - erfc(x * inv_sqrt2));
    return cdf;
}

double MillsRatio(long double x) {
    double r;
    // r = (1 - CDF_Normal(x)) / PDF_Normal(x);
    r = 0.5 * erfcx(x * inv_sqrt2) * sqrt_2pi;

    if (x >= 0)
        r = 0.5 * erfcx(x * inv_sqrt2) * sqrt_2pi;
    else
        r = 0.5 * erfcx(x * inv_sqrt2) / PDF_Normal(x);
    return r;
}

double Pdf_x_MillsRatio(long double x, long double t) {
    // double r = (1 - CDF_Normal(t)) * (PDF_Normal(x) / PDF_Normal(t));
    double lam = t - x;
    double ncdf;
    if (t > 0)
        ncdf = CDF_Normal(-t);
    else
        ncdf = (1 - CDF_Normal(t));

    // double r = ncdf * exp(lam*(lam*0.5 + x));
    double r = exp(log(ncdf) + lam * (lam * 0.5 + x));
    return r;
}

double Pdf_x_MillsRatio(long double x, long double y, long double alpha_sig) {
    double tx = alpha_sig - x;
    double ty = alpha_sig - y;

    double lam = exp(-alpha_sig * (y - x));

    double ncdf = 0.5 * (lam * erfc(ty * inv_sqrt2) - erfc(tx * inv_sqrt2));

    double r;
    if (ncdf == 0)
        r = 0;
    else if (ncdf < 0)
        r = -exp(log(-ncdf) + alpha_sig * (alpha_sig * 0.5 - x));
    else
        r = exp(log(ncdf) + alpha_sig * (alpha_sig * 0.5 - x));

    return r;
}

double CDF_Normal2(double xn, double yn) {
    double pnx = PDF_Normal(xn);
    double pny = PDF_Normal(yn);
    double mrx = MillsRatio(fabs(xn));
    double mry = MillsRatio(fabs(yn));

    double cx, cy;

    double diff;

    cx = pnx * mrx;
    cy = pny * mry;

    if (xn > 0 && yn > 0)
        diff = cx - cy;
    else if (xn < 0 && yn < 0) // (1-cy) - (1-cx) ==> cx - cy
        diff = cy - cx;
    else {
        if (xn < 0) {
            cx = 1.0 - cx;
        }

        if (yn < 0) {
            cy = 1.0 - cy;
        }

        diff = cx - cy;
    }

    return diff;
}

double PDF_Exp(double x, double mu, double alpha) {
    double xn = x - mu;
    double y = alpha * exp(-xn * alpha);
    return y;
}

double CDF_Exp(double x, double mu, double alpha) {
    double y;
    double xn = x - mu;
    y = 1 - exp(-xn * alpha);
    assert(xn >= 0); // debug.
    return y;
}

double CDF_Exp(double x, double y, double mu, double alpha) {
    double cdf = exp(-alpha * (x - mu)) * (1.0 - exp(-alpha * (y - x)));
    assert(x >= mu);
    return cdf;
}

double PDF_Laplace(double x, double mu, double alpha, double beta) {
    double y;
    double xn = x - mu;
    if (xn < 0)
        y = beta * exp(beta * xn) / (1.0 + (beta / alpha));
    else
        y = alpha * exp(-alpha * xn) / (1.0 + (alpha / beta));

    return y;
}

double CDF_Laplace(double x, double mu, double alpha, double beta) {
    double y;
    double xn = x - mu;
    if (xn < 0)
        y = exp(beta * xn) / (1.0 + (beta / alpha));
    else
        y = 1 - exp(-alpha * xn) / (1.0 + (alpha / beta));

    return y;
}

double CDF_Laplace(double x, double y, double mu, double alpha, double beta) {
    double cdf;
    double xn = x - mu;
    double yn = y - mu;

    if (xn <= 0 && yn < 0) {
        cdf = exp(beta * yn) - exp(beta * xn);
        cdf *= (1.0 / (1.0 + (beta / alpha)));
    } else if (xn >= 0 && yn > 0) {
        cdf = -exp(-alpha * yn) / (1.0 + (alpha / beta));
        cdf = -(exp(-alpha * yn) - exp(-alpha * xn));
        cdf *= (1.0 / (1.0 + (alpha / beta)));
    } else {
        cdf = 1 - (exp(-alpha * yn) / (1.0 + (alpha / beta)) - exp(beta * xn) / (1.0 + (beta / alpha)));
    }

    return cdf;
}

double PDF_NormalLaplace(double x, double mu, double sigma, double alpha, double beta) {
    double xn = (x - mu) / sigma;

    double alpha_sig;
    double beta_sig;

    if (sigma == 0) {
        alpha_sig = 0;
        beta_sig = 0;
    } else {
        alpha_sig = alpha * sigma;
        beta_sig = beta * sigma;
    }

    double tmp1x = alpha_sig - xn;
    double tmp2x = beta_sig + xn;

    double pnx = PDF_Normal(xn);

    double mrlx, mr2x;

    if (tmp1x >= 0) {
        mrlx = pnx * MillsRatio(tmp1x);
    } else {
        mrlx = exp(0.5 * alpha_sig * alpha_sig - alpha * (x - mu)) - pnx * MillsRatio(-tmp1x);
    }

    if (tmp2x >= 0) {
        mr2x = pnx * MillsRatio(tmp2x);
    } else {
        mr2x = exp(0.5 * beta_sig * beta_sig + beta * (x - mu)) - pnx * MillsRatio(-tmp2x);
    }

    double pdf = (mrlx + mr2x) * (1.0 / (1.0 / alpha + 1.0 / beta));

    if (isnan(pdf) || isinf(pdf) || pdf < 0) {
        printf("ERROR!!\n");
#ifdef _WIN32
        DebugBreak();
#endif
    }

    return pdf;
}

double CDF_NormalLaplace(double x, double mu, double sigma, double alpha, double beta) {
    double xn = (x - mu) / sigma;

    double tmp1 = alpha * sigma - xn;
    double tmp2 = beta * sigma + xn;

    double cn = CDF_Normal(xn);
    double pn = PDF_Normal(xn);
    double mr1 = MillsRatio(tmp1);
    double mr2 = MillsRatio(tmp2);
    double y;

    if (isinf(mr1) || isinf(mr2))
        return cn;
    else {
        // y = cn - pn * (beta*mr1 - alpha * mr2) / (alpha + beta);
        double y1 = (1.0 / (1.0 + (alpha / beta))) * pn * mr1;
        double y2 = (1.0 / (1.0 + (beta / alpha))) * pn * mr2;
        if (isnan(y1))
            y1 = 0;
        if (isnan(y2))
            y2 = 0;
        y = cn - y1 + y2;

        if (isnan(y) || isinf(y)) {
            printf("ERROR!!\n");
#ifdef _WIN32
            DebugBreak();
#endif
        }
    }

    return y;
}

double CDF_NormalLaplace(double x, double y, double mu, double sigma, double alpha, double beta) {
    double xn = (x - mu) / sigma;
    double yn = (y - mu) / sigma;

    double alpha_sig = alpha * sigma;
    double beta_sig = beta * sigma;

    double tmp1x = alpha_sig - xn;
    double tmp2x = beta_sig + xn;
    double tmp1y = alpha_sig - yn;
    double tmp2y = beta_sig + yn;

    double delta = yn - xn;

    double diff;

    double pnx = PDF_Normal(xn);
    double pny = PDF_Normal(yn);

#if 0
    diff = CDF_Normal(xn, yn);
    if (diff == 0)
        diff = PDF_Normal(xn) * delta;
#endif

#if 1
    diff = CDF_Normal2(xn, yn);
#endif

    double mr1, mr2, mrlx, mrly, mr2x, mr2y;

    if (tmp1x >= 0) {
        mrlx = pnx * MillsRatio(tmp1x);
    } else {
        //mrlx = exp(alpha_sig*(0.5*alpha_sig - xn)) - pnx * MillsRatio(-tmp1x);
        mrlx = exp(0.5 * alpha_sig * alpha_sig - alpha * (x - mu)) - pnx * MillsRatio(-tmp1x);
    }

    if (tmp1y >= 0) {
        mrly = pny * MillsRatio(tmp1y);
    } else {
        //mrly = exp(alpha_sig*(0.5*alpha_sig - yn)) - pny * MillsRatio(-tmp1y);
        mrly = exp(0.5 * alpha_sig * alpha_sig - alpha * (y - mu)) - pny * MillsRatio(-tmp1y);
    }

    if (tmp2x >= 0) {
        mr2x = pnx * MillsRatio(tmp2x);
    } else {
        //mr2x = exp(beta_sig*(0.5*beta_sig + xn)) - pnx * MillsRatio(-tmp2x);
        mr2x = exp(0.5 * beta_sig * beta_sig + beta * (x - mu)) - pnx * MillsRatio(-tmp2x);
    }

    if (tmp2y >= 0) {
        mr2y = pny * MillsRatio(tmp2y);
    } else {
        //mr2y = exp(beta_sig*(0.5*beta_sig + yn)) - pny * MillsRatio(-tmp2y);
        mr2y = exp(0.5 * beta_sig * beta_sig + beta * (y - mu)) - pny * MillsRatio(-tmp2y);
    }

    mr1 = mrly - mrlx;
    mr2 = mr2y - mr2x;

    double d1 = (1.0 / ((alpha / beta) + 1.0)) * mr1;
    double d2 = (1.0 / ((beta / alpha) + 1.0)) * mr2;

    double cdf = diff - d1 + d2;

    if (isnan(cdf) || isinf(cdf) || cdf < 0) {
        printf("ERROR!!\n");
#ifdef _WIN32
        DebugBreak();
#endif
    }

    return cdf;
}

static double TargetFunction(std::vector<double>& params) {
    int num_states = normal_laplace_mix_params.num_states;

    ParamDeserialize(normal_laplace_mix_params, params); /// !!!!

    double sum = 0.0;
    for (int i = 0; i < num_states; ++i) {
        int num_bins = static_cast<int>(data_x[i].size());
        for (int k = 0; k < num_bins; ++k) {
            if (data_y[i][k] == 0)
                continue;

            /*TO BE FIXED*/
            double x = data_x[i][k];
            double p = data_y[i][k];

            double Qk = 0.0;
            for (int j = 0; j < num_states; ++j) {
                if (normal_laplace_mix_params.lambda[i][j] == 0)
                    continue;

                double mu = normal_laplace_mix_params.mu[j];
                double sigma = normal_laplace_mix_params.sigma[j];
                double alpha = normal_laplace_mix_params.alpha[j];
                double beta = normal_laplace_mix_params.beta[j];

                double phi;

                if (isinf(alpha) && isinf(beta)) {
                    // Normal
                    phi = PDF_Normal((x - mu) / sigma);
                } else if (sigma == 0 || sigma < 1E-3) {
                    // Exp, Laplace
                    if (isinf(alpha)) {
                        if (beta < 0)
                            phi = INFINITY;
                        else if (x > mu)
                            phi = 0;
                        else
                            phi = PDF_Exp(-x, -mu, beta);
                    } else if (isinf(beta)) {
                        if (alpha < 0)
                            phi = INFINITY;
                        else if (x < mu)
                            phi = 0;
                        else
                            phi = PDF_Exp(x, mu, alpha);
                    } else {
                        if (alpha < 0 || beta < 0)
                            phi = INFINITY;
                        else
                            phi = PDF_Laplace(x, mu, alpha, beta);
                    }
                    assert(phi >= 0);
                } else {
                    // Normal-Laplace
                    if (alpha < 0 || beta < 0 || sigma < 0)
                        phi = INFINITY;
                    else
                        phi = PDF_NormalLaplace(x, mu, sigma, alpha, beta);
                }

                if (isnan(phi))
#ifdef _WIN32
                    DebugBreak();
#else
                    ;
#endif

                if (isinf(phi) && normal_laplace_mix_params.lambda[i][j] == 0)
#ifdef _WIN32
                    DebugBreak();
#else
                    ;
#endif

                // if (isinf(phi) && normal_laplace_mix_params.lambda[i][j] != 0)
                //     DebugBreak();

                if (normal_laplace_mix_params.lambda[i][j] != 0 &&
                    (beta < 0 || alpha < 0 || sigma < 0 || normal_laplace_mix_params.lambda[i][j] < 0))
                    Qk = INFINITY;
                else if (normal_laplace_mix_params.lambda[i][j] == 0 && isinf(phi)) {
                } else
                    Qk += normal_laplace_mix_params.lambda[i][j] * phi;
            }

            if (Qk < 0) {
                printf("Qk < 0 !!!ERROR!\n");
#ifdef _WIN32
                DebugBreak();
#endif
                continue;
            }

            double diff;

            if (method_parameters.dist_function_mode == 0) {
                diff = log(p) - log(Qk);
                if (isinf(Qk))
                    diff = INFINITY;
                sum += p * diff; // Kullback-Leibler divergence (Dkl)
            } else if (method_parameters.dist_function_mode == 1) {
                diff = log(p) - log(Qk);
                sum += p * fabs(diff); // absolute of Dkl.
            } else if (method_parameters.dist_function_mode == 2) {
                diff = log(p) - log(Qk);
                sum += p * diff * diff; // plogsq
            } else if (method_parameters.dist_function_mode == 3) {
                diff = log(p) - log(Qk);
                sum += diff * diff; // logsq
            } else {
                diff = Qk - p;
                sum += diff * diff; // sq
            }

            if (isnan(sum)) {
                printf("isnan(sum) p=%.7g Qk=%.7g diff=%.7g !!!ERROR!\n", p, Qk, diff);
#ifdef _WIN32
                DebugBreak();
#endif
            }
        }
    }

    return sum;
}

// bool CheckDimension(std::vector<Point> simplex, int num_vertices, int num_params)
// {
//     double ** matrix = new double*[num_params];
//     for (int i = 0; i < num_vertices; ++i) {
//         matrix[i] = new double[num_vertices];
//     }
//
//     for (int i = 0; i < num_params; ++i) {
//         for (int j = 0; j < num_vertices; ++j) {
//             matrix[i][j] = simplex[j][i];
//         }
//     }
//
//     if (rank(matrix) == num_params) return true;
//     else return false;
// }
//
// bool VerifyMin(FUNCTION function, Point optima, double tolerance, int num_params, double min, double step_size)
// {
//     for (int i = 0; i < num_params; ++i) {
//         optima[i] += step_size;
//     }
//
//     double new_value = function(optima);
//     double diff = fabs(new_value - min);
//     if (diff < tolerance) return true;
//     else return false;
// }
