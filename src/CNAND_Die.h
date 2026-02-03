#pragma once

#include "stdafx.h"
#include "SysTypes.h"
#include "CVthDistribution_base.h"
#include "CNAND_Block.h"

class CNAND_Die {
public:
    CNAND_Die(CellType cell_type, BYTE* gray_codes, int blocks, int wordlines, int cells_per_wordline,
             CVthDistribution_base *vth_dist, double *thresh_volts);

    // Extended constructor used by NandModelTest (id + RTN parameters).
    CNAND_Die(int id, CellType cell_type, BYTE* gray_codes, int blocks, int wordlines, int cells_per_wordline,
             CVthDistribution_base *vth_dist, double *thresh_volts, int rtn_mode = 0, double rtn_sigma = 0.0);

    ~CNAND_Die();

    void Init();

    bool ProgramOneShot(int block, int wordline, BYTE* data_in);
    void ReadPage(int block, int wordline, int bit_idx, BYTE* buffer);

    bool EraseBlock(int block);
    void Evolve(int block);
    void SetVth(double *vth);

    void DbgDumpVoltage();
    void DbgDumpVoltage(int block);
    void DbgDumpVoltage(int block, int wordline);

    int Id() const { return m_id; }

private:
    int m_id = 0;
    int m_total_blocks = 0;
    CNAND_Block **m_blocks = nullptr;

    int m_bits_per_cell = 0;
    int m_num_states = 0;

    // Cell states voltage distribution.
    CVthDistribution_base *m_vth_dist = nullptr;
    double *m_thresh_volts = nullptr;

    // Input values, represented using gray coding scheme.
    BYTE *m_gray_codes = nullptr;

    // Number of program/erase cycles.
    int PEC = 0;

    int time = 0;
    int temperature = 0;
};

