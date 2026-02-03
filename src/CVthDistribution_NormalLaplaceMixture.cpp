#include "CVthDistribution_NormalLaplaceMixture.h"

#include <cmath>
#include "Logger.h"

using namespace std;

CVthDistribution_NormalLaplaceMixture::CVthDistribution_NormalLaplaceMixture(int cell_bits, BYTE* values) {
    m_bits_per_cell = cell_bits;
    m_states = 1 << cell_bits;

    m_values = new BYTE[m_states];
    m_mu = new double[m_states];
    m_sigma = new double[m_states];
    m_alpha_inv = new double[m_states];
    m_beta_inv = new double[m_states];

    m_uniform_distribution = new std::uniform_real_distribution<double>*[m_states];
    m_normal_distribution = new std::normal_distribution<double>*[m_states];

    for (int i = 0; i < m_states; ++i) {
        int v = values[i];
        m_uniform_distribution[v] = new std::uniform_real_distribution<double>(0, 1);
        m_normal_distribution[v] = new std::normal_distribution<double>(0, 1);
    }

    m_generator = new std::default_random_engine();

    m_lambda = new double*[m_states];
    for (int i = 0; i < m_states; ++i) {
        m_lambda[i] = new double[m_states];
        for (int j = 0; j < m_states; ++j) {
            m_lambda[i][j] = 0.0;
        }
    }
}

CVthDistribution_NormalLaplaceMixture::~CVthDistribution_NormalLaplaceMixture() {
    delete[] m_values;
    delete[] m_mu;
    delete[] m_sigma;
    delete[] m_alpha_inv;
    delete[] m_beta_inv;

    delete m_generator;

    if (m_lambda) {
        for (int i = 0; i < m_states; ++i) {
            delete[] m_lambda[i];
        }
        delete[] m_lambda;
    }

    if (m_uniform_distribution) {
        for (int i = 0; i < m_states; ++i) {
            delete m_uniform_distribution[i];
        }
        delete[] m_uniform_distribution;
    }
    if (m_normal_distribution) {
        for (int i = 0; i < m_states; ++i) {
            delete m_normal_distribution[i];
        }
        delete[] m_normal_distribution;
    }
}

void CVthDistribution_NormalLaplaceMixture::Init(
    BYTE *values,
    double pe_cycle_L,
    std::vector<double>& c0_sigma,
    std::vector<double>& c0_lambda,
    std::vector<double>& c1_lambda,
    std::vector<double>& c2_lambda,
    std::vector<double>& c0_mu,
    std::vector<double>& c1_mu,
    std::vector<double>& c2_mu,
    std::vector<double>& c0_alpha_inv,
    std::vector<double>& c1_alpha_inv,
    std::vector<double>& c0_beta_inv,
    std::vector<double>& c1_beta_inv
) {
    m_values = values;

    for (int i = 0; i < m_states; ++i) {
        int vi = values[i];
        m_lambda[vi][vi] = 1.0;
        for (int j = 0; j < m_states; ++j) {
            int vj = values[j];
            m_lambda[vi][vj] = c2_lambda[i * m_states + j] *
                exp(c1_lambda[i * m_states + j] * pe_cycle_L + c0_lambda[i * m_states + j]);
            m_lambda[vi][vi] -= m_lambda[vi][vj];
        }
    }

    m_mu[values[0]] = pe_cycle_L + c0_mu[0];
    for (int i = 1; i < m_states; ++i) {
        int v = values[i];
        m_mu[v] = c1_mu[i] * pow(pe_cycle_L, c2_mu[i]) + c0_mu[i];
    }
    for (int i = 0; i < m_states; ++i) {
        LogDebug(2, "m_mu[%d] = %.3f", i, m_mu[i]);
    }

    for (int i = 0; i < m_states; ++i) {
        int v = values[i];
        m_sigma[v] = c0_sigma[i];
    }
    for (int i = 0; i < m_states; ++i) {
        LogDebug(2, "m_sigma[%d] = %.3f", i, m_sigma[i]);
    }

    for (int i = 0; i < m_states; ++i) {
        int v = values[i];
        m_alpha_inv[v] = c1_alpha_inv[i] * pe_cycle_L + c0_alpha_inv[i];
        m_beta_inv[v] = c1_beta_inv[i] * pe_cycle_L + c0_beta_inv[i];
    }
    for (int i = 0; i < m_states; ++i) {
        LogDebug(2, "m_alpha_inv[%d] = %.3f m_beta_inv[%d]=%.3f", values[i], m_alpha_inv[values[i]], m_beta_inv[values[i]]);
    }
}

double CVthDistribution_NormalLaplaceMixture::GenerateVoltageInternal(BYTE value) {
    double random_normal_var = (*m_normal_distribution[value])(*m_generator);
    double uniform_var1 = (*m_uniform_distribution[value])(*m_generator);
    double uniform_var2 = (*m_uniform_distribution[value])(*m_generator);

    while (uniform_var1 == 1.0) {
        uniform_var1 = (*m_uniform_distribution[value])(*m_generator);
    }
    while (uniform_var2 == 1.0) {
        uniform_var2 = (*m_uniform_distribution[value])(*m_generator);
    }

    double std_exp_var1 = -log(1 - uniform_var1);
    double std_exp_var2 = -log(1 - uniform_var2);

    double voltage = m_mu[value] + m_sigma[value] * random_normal_var +
        m_alpha_inv[value] * std_exp_var1 - m_beta_inv[value] * std_exp_var2;

    return voltage / 100.0;
}

double CVthDistribution_NormalLaplaceMixture::GenerateVoltage(BYTE value) {
    double u = (*m_uniform_distribution[value])(*m_generator);
    double lambda_acc = 0.0;

    for (int val = 0; val < m_states; ++val) {
        lambda_acc += m_lambda[value][val];
        if (u <= lambda_acc) {
            return GenerateVoltageInternal(static_cast<BYTE>(val));
        }
    }

    return GenerateVoltageInternal(value);
}

