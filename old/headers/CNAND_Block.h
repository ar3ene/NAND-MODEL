#pragma once
#include "stdfax.h"
#include "sys_types.h"
#include "CVthDistribution_Base.h"

class CNAND_Block {
public:
	CNAND_Block(CellType cell_type, BYTE* gray_codes,
		int wordlines, int cells_per_wordline, 
		CVthDistribution_base* vth_dist, double* thresh_volts);
	~CNAND_Block();

	void Init(CVthDistribution_base* vth_dist, BYTE* grady_codes);

	bool Erase();

	bool ProgramOneShot(int wordline, BYTE* data_in);

	void ReadPage(int wordline, int bit_idx, BYTE* buffer,
		double* thresh_volts, BYTE* gray_codes);

	void Evolve(double& cell_voltage, double v0, double t, double t0, double A,
		double C, double tau, double Taf, double Eaf, double f);

	void DbgDumpVoltage();

	void DbgDumpVoltage(int wordline);

private:
	CellType m_cell_type;
	int m_wordlines_per_block;
	int m_cells_per_wordline;
	int m_total_cells;
	int m_values_per_cell;
	BYTE** cell_values;
	double** cell_voltages;
	double*** cell_voltage_table;

	bool Voltage2Bit(double volt, int bit_idx, double* vth, BYTE* values);
};