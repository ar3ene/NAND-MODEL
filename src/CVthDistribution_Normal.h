#pragma once

#include "stdafx.h"
#include <random>

#include "SysTypes.h"
#include "CVthDistribution_base.h"

class CVthDistribution_Normal : public CVthDistribution_base {
public:
    explicit CVthDistribution_Normal(int bits);
    ~CVthDistribution_Normal() override;

    // Initialize distribution parameters.
    void Init(double *means, double *standard_deviations, BYTE *values);

    // Generate a voltage according to a given value.
    double GenerateVoltage(BYTE value) override;

private:
    int m_bits_per_cell = 0;
    int m_states = 0;
    double *m_means = nullptr;
    double *m_standard_deviations = nullptr;
    BYTE *m_values = nullptr;

    std::default_random_engine *m_rand_generator = nullptr;
    std::normal_distribution<double> **m_normal_distribution = nullptr;
};

