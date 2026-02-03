#pragma once

#include "stdafx.h"
#include "SysTypes.h"

class CVthDistribution_base {
public:
    virtual ~CVthDistribution_base() = default;

    // Pure virtual: generate a voltage for a given state value.
    virtual double GenerateVoltage(BYTE value) = 0;
};

