# Classical MC Simple (Step 1)

This directory is a minimal, beginner-friendly equilibrium classical Monte Carlo implementation
without MPI.

**Quick Start**
1. Build:
```bash
cd stage1
cmake -S . -B build
cmake --build build -j
```
2. Run the provided 2D Ising sample (`L=16`):
```bash
cd samples/square_L16_Ising
../../build/MC_simple
```
3. Check results:
```bash
cat MC_simple_result.dat
```

**How To Use**
1. Default run (reads `param.def`, `lattice.def`, `interaction.def` in current directory):
```bash
./MC_simple
```
2. Specify only `param.def`:
```bash
./MC_simple my_param.def
```
3. Specify all three input files:
```bash
./MC_simple my_param.def my_lattice.def my_interaction.def
```
4. Optional reproducibility seed:
```bash
MC_SIMPLE_SEED=12345 ./MC_simple
```

**What It Computes**
1. `E_per_site`: thermal average of energy per site.
2. `C_per_site`: specific heat per site,
`C = (<E^2> - <E>^2) / (N * T^2)`.
3. `M2`: thermal average of magnetization squared per site,
`M2 = <(Mx^2 + My^2 + Mz^2) / N^2>`.
4. `acceptance`: Metropolis acceptance ratio.

Output columns in `MC_simple_result.dat`:
1. `T`
2. `E_per_site`
3. `C_per_site`
4. `M2`
5. `acceptance`

**Input Files**
1. `param.def`: MC parameters (`Burn_in`, `Total_Step`, `Sample`, `num_temp`, `Ini_T`, `Delta_T`, `spin_dim`, `lambda`, `H`).
2. `lattice.def`: lattice size (`L_x`, `L_y`, `L_z`) and orbital count (`orb_num`).
3. `interaction.def`: bond list (`i_orb j_orb dx dy dz J`).

The sample case is in `samples/square_L16_Ising/`.

**Code Structure (Detailed)**
1. `main.c`
Main control flow:
argument parsing, input reading, simulation loop, averaging, and final output.
2. `memory.c` and `memory.h`
Allocation and cleanup for runtime arrays and accumulation buffers.
3. `input_parser.c`
Parses `param.def`, `lattice.def`, and `interaction.def`.
4. `lattice.c`
Builds neighbor table (`Nr`, `J`), initializes spins, computes initial local fields and energy.
5. `mc_update.c`
One Metropolis sweep per temperature slot with incremental local-field updates.
6. `mc_def.h`
Shared structs and function prototypes used by stage1 files.
7. `dSFMT.c`, `dSFMT.h`, `dSFMT-params.h`, `dSFMT-params19937.h`
Random number generator implementation (dSFMT, MEXP=19937 only in Stage 1).
8. `CMakeLists.txt`
Build definition for the single executable `MC_simple`.

**Simulation Flow**
1. Read and validate input files.
2. Allocate arrays for spins, local fields, neighbors, and statistics buffers.
3. For each sample:
set RNG seed, initialize spins and fields, run burn-in, then run measurement steps.
4. During measurement:
accumulate `E`, `E^2`, `M2`, and acceptance.
5. Convert step averages to sample averages.
6. Average over samples and write final table.

**Deliberately Omitted In Stage 1**
1. MPI parallel execution.
2. Replica exchange MC.
3. Large output/aggregation modules from the full refactoring code.

This keeps Step 1 compact and readable for first-time users.

**License**
1. This Stage 1 code is released under the MIT License (see `LICENSE`).
2. dSFMT files (`dSFMT.c`, `dSFMT.h`, `dSFMT-params.h`, `dSFMT-params19937.h`)
are third-party code under BSD 3-Clause (see `LICENSE.txt`).
