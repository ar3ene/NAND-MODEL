#include "CNAND_Die.h"

#include <memory.h>
#include <random>

CNAND_Die::CNAND_Die(CellType cell_type, BYTE* gray_codes, int blocks, int wordlines, int cells_per_wordline,
                     CVthDistribution_base *vth_dist, double *thresh_volts)
    : CNAND_Die(0, cell_type, gray_codes, blocks, wordlines, cells_per_wordline, vth_dist, thresh_volts, 0, 0.0) {
}

CNAND_Die::CNAND_Die(int id, CellType cell_type, BYTE* gray_codes, int blocks, int wordlines, int cells_per_wordline,
                     CVthDistribution_base *vth_dist, double *thresh_volts, int /*rtn_mode*/, double /*rtn_sigma*/) {
    m_id = id;
    m_bits_per_cell = static_cast<int>(cell_type);
    m_num_states = (1 << m_bits_per_cell);

    m_vth_dist = vth_dist;
    m_gray_codes = gray_codes;
    m_total_blocks = blocks;

    m_blocks = new CNAND_Block*[blocks];

    int thresh_volts_size = m_num_states - 1;
    m_thresh_volts = new double[thresh_volts_size];
    memcpy(m_thresh_volts, thresh_volts, thresh_volts_size * sizeof(*m_thresh_volts));

    m_gray_codes = new BYTE[m_num_states];
    memcpy(m_gray_codes, gray_codes, m_num_states);

    for (int i = 0; i < blocks; i++) {
        m_blocks[i] = new CNAND_Block(cell_type, gray_codes, wordlines, cells_per_wordline, vth_dist, thresh_volts);
    }
}

CNAND_Die::~CNAND_Die() {
    for (int i = 0; i < m_total_blocks; i++) {
        delete m_blocks[i];
    }
    delete[] m_blocks;
    delete[] m_thresh_volts;
    delete[] m_gray_codes;
}

void CNAND_Die::Init() {
    for (int i = 0; i < m_total_blocks; i++) {
        m_blocks[i]->Init(m_vth_dist, m_gray_codes);
    }
}

bool CNAND_Die::ProgramOneShot(int block, int wordline, BYTE* data_in) {
    m_blocks[block]->ProgramOneShot(wordline, data_in);
    return true;
}

void CNAND_Die::ReadPage(int block, int wordline, int bit_idx, BYTE* buffer) {
    m_blocks[block]->ReadPage(wordline, bit_idx, buffer, m_thresh_volts, m_gray_codes);
}

bool CNAND_Die::EraseBlock(int block) {
    return m_blocks[block]->Erase();
}

void CNAND_Die::Evolve(int block) {
    // Placeholder: evolution model not fully captured in the fragments.
    (void)block;
}

void CNAND_Die::SetVth(double *vth) {
    memcpy(m_thresh_volts, vth, sizeof(*vth) * (m_num_states - 1));
}

void CNAND_Die::DbgDumpVoltage() {
    for (int i = 0; i < m_total_blocks; i++) {
        m_blocks[i]->DbgDumpVoltage();
    }
}

void CNAND_Die::DbgDumpVoltage(int block) {
    m_blocks[block]->DbgDumpVoltage();
}

void CNAND_Die::DbgDumpVoltage(int block, int wordline) {
    m_blocks[block]->DbgDumpVoltage(wordline);
}

