#pragma once
#include "stdfax.h"
#include <random>
#include "sys_types.h"
#include "CVthDistribution_Base.h"

class CVthDistribution_Normal : public CVthDistribution_base {
public:
	CVthDistribution_Normal(int bits);
	~CVthDistribution_Normal();

	void Init(double* means, double* standard_deviations, BYTE* values);

	double GenerateVoltage(BYTE value) override;

private:
	int m_bits_per_cell;
	int m_states;
	double* m_means;
	double* m_standard_deviations;
	BYTE* m_values;
	default_random_engine* m_rand_generator;
	normal_distribution<double>** m_normal_distribution;
};
