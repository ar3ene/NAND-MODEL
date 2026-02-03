#pragma once

#ifndef _NAND_DIE_MODEL_H_
#define _NAND_DIE_MODEL_H_

#include "stdafx.h"
#include "sys_types.h"
#include "CVthDistribution_base.h"
#include "CNAND_Block.h"

class CNAND_Die {
public:
	CNAND_Die(CellType cell_type, BYTE* gray_codes, int blocks, int wordlines, int cells_per_wordline,
		CVthDistribution_base* vth_dist, double* thresh_volts);
	~CNAND_Die();

	void Init();

	/**
	* Performing "One Shot Prgram" to the cell
	* @param block an integer argument, the specific block to be programmed.
	* @param wordline an integer argument, the specific wordline to be programmed.
	* @param data_in an argument of type BYTE, values to be programmed to the cell.
	* @see ProgramOneShot().
	* @return void.
	*/
	bool ProgramOneShot(int block, int wordline, BYTE* data_in);

	/**
	* @see void ReadPage(int block, int wordline, int bit_idx, BYTE* buffer)
	* @param block an integer argument, the specific block to be programmed.
	* @param wordline an integer argument, specifies which bit of the cell value to be read.
	* @param bit_idx an int argument, values to be programmed to the cell.
	* @param buffer an array of type BYTE, stores the bit value read from the cell.
	* @return void.
	*/
	void ReadPage(int block, int wordline, int bit_idx, BYTE* buffer);

	bool EraseBlock(int block);

	void Evolve(int block);

	void SetVth(double* vth);

	void DbgDumpVoltage();
	void DbgDumpVoltage(int block);
	void DbgDumpVoltage(int block, int wordline);

private:
	int m_total_blocks;
	CNAND_Block** m_blocks;
	int m_bits_per_cell;
	int m_num_states;

	CVthDistribution_base* m_vth_dist;
	double* m_thresh_volts;

	BYTE* m_gray_codes;

	int PEC; ///< # program-erase cycles.

	int time;

	int temperature;
};

#endif // ! _NAND_DIE_MODEL_H_