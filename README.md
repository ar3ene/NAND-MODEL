# NAND-MODEL
Simulation model for NAND flash controller behavior and cell threshold‑voltage (Vth) distributions.

## Overview
- Simulates program/read/erase at die/block/wordline granularity.
- Models Vth with Normal and Normal‑Laplace mixture distributions.
- Sweeps read thresholds to build CDF/PDF bins.
- Computes transition matrix, BER, read‑one ratio, bit asymmetry, and data‑in BVA.
- Includes a Nelder–Mead based parameter fitting tool.

## Repository layout
- `src/` Core NAND model + fitting utilities.
- `legacy/` Legacy files from upstream (kept for reference).
- `LICENSE` Non‑commercial license.

## Build
This code uses Windows APIs (CreateThread, Interlocked, Mutex, DebugBreak) and is easiest to build with MSVC.
Compile all sources in `src/` (including `erfcx.c`).

## Usage

### NandModelTest
```bash
NandModelTest.exe <config file>
```

### NLMPdfCurveFitting
```bash
NLMPdfCurveFitting.exe <config file> <data file>
```
Note: the usage string in code may show only one argument, but the program expects two.

## Config keys (NandModelTest)
Required keys used by `NandModelTest.cpp`:
- `time`, `t0`, `A`, `C`, `tau`, `Taf`, `Eaf`, `f`
- `cdf_mode`, `RTN_mode`, `RTN_sigma`
- `DebugLevel`, `num_threads`
- `cell_type` (bits per cell), `total_dies`, `total_number_of_blocks`,
  `number_of_wordlines_per_block`, `number_of_cells_per_wordline`
- `state_values` (size = states)
- `threshold_voltages` (size = states - 1)
- `means`, `standard_deviations` (size = states)
- `vth_range`, `vth_step` (size = states - 1)
- `c0_lambda`, `c1_lambda`, `c2_lambda` (size = states * states)
- `c0_mu`, `c1_mu`, `c2_mu` (size = states)
- `c0_sigma`, `c0_alpha`, `c1_alpha`, `c0_beta`, `c1_beta` (size = states)

## Config keys (NLMPdfCurveFitting)
Required keys used by `NLMPdfCurveFitting.cpp`:
- `fitting_mode`, `dist_function_mode`, `trace_path`, `max_iterations`
- `num_states`
- `lambda`, `lambda_mask`, `lambda_step`
- `mu`, `mu_mask`, `mu_step`
- `sigma`, `sigma_mask`, `sigma_step`
- `alpha`, `alpha_mask`, `alpha_step`
- `beta`, `beta_mask`, `beta_step`

## Data file format (NLMPdfCurveFitting)
- Each line: `x y` (whitespace separated).
- Lines starting with `#` are skipped.
- A line starting with `&` increments the state index.
- Whitespace delimiters: space / tab / CR / LF.

## Outputs
- CDF or PDF as voltage/value pairs.
- Transition matrix and per‑bit BER.
- Read‑one ratio statistics and data‑in BVA metrics.

## Notes
- Multi‑threaded per‑die processing; a mutex guards bin area accumulation.
- Some files in `src/` are placeholders (e.g., `ModelCurveFitting`, `NandDriver`, `PAConsole`) and may need completion depending on usage.
