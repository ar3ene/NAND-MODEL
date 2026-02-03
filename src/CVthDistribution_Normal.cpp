#include "CVthDistribution_Normal.h"

#include <cmath>
#include <cstring>

CVthDistribution_Normal::CVthDistribution_Normal(int bits) {
    m_bits_per_cell = bits;
    m_states = 1 << bits;
    m_means = new double[m_states];
    m_standard_deviations = new double[m_states];
    m_values = new BYTE[m_states];
    m_rand_generator = new std::default_random_engine();
    m_normal_distribution = nullptr;
}

CVthDistribution_Normal::~CVthDistribution_Normal() {
    if (m_normal_distribution) {
        for (int i = 0; i < m_states; ++i) {
            delete m_normal_distribution[i];
        }
        delete[] m_normal_distribution;
    }
    delete m_rand_generator;
    delete[] m_means;
    delete[] m_standard_deviations;
    delete[] m_values;
}

void CVthDistribution_Normal::Init(double *means, double *standard_deviations, BYTE *values) {
    if (m_normal_distribution) {
        for (int i = 0; i < m_states; ++i) {
            delete m_normal_distribution[i];
        }
        delete[] m_normal_distribution;
    }

    m_normal_distribution = new std::normal_distribution<double>*[m_states];

    for (int i = 0; i < m_states; ++i) {
        int v = values[i];
        double avg = means[i];
        double sigma = standard_deviations[i];
        m_values[i] = static_cast<BYTE>(v);
        m_means[i] = avg;
        m_standard_deviations[i] = sigma;
        m_normal_distribution[v] = new std::normal_distribution<double>(avg, sigma);
    }
}

double CVthDistribution_Normal::GenerateVoltage(BYTE value) {
    return (*m_normal_distribution[value])(*m_rand_generator);
}

