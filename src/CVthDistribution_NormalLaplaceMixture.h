#pragma once

#include "stdafx.h"
#include <random>
#include <vector>

#include "SysTypes.h"
#include "CVthDistribution_base.h"

class CVthDistribution_NormalLaplaceMixture : public CVthDistribution_base {
public:
    CVthDistribution_NormalLaplaceMixture(int cell_bits, BYTE* values);
    ~CVthDistribution_NormalLaplaceMixture() override;

    void Init(
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
    );

    double GenerateVoltage(BYTE value) override;

private:
    double GenerateVoltageInternal(BYTE value);

    double *m_mu = nullptr;
    double *m_sigma = nullptr;
    double *m_alpha_inv = nullptr;
    double *m_beta_inv = nullptr;
    BYTE *m_values = nullptr;
    int m_bits_per_cell = 0;
    int m_states = 0;
    double **m_lambda = nullptr;

    std::default_random_engine *m_generator = nullptr;
    std::uniform_real_distribution<double> **m_uniform_distribution = nullptr;
    std::normal_distribution<double> **m_normal_distribution = nullptr;
};

