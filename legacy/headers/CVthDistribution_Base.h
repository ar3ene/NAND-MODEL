#pragma once
#include "stdafx.h"
#include "sys_types.h"
#include <random>

class CVthDistribution_base {
public:
	virtual double GenerateVoltage(BYTE value) = 0;
};