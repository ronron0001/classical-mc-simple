# classical-mc-simple

Minimal, beginner-friendly classical spin Monte Carlo code (single process, no MPI).

The current implementation is in [`src/`](./src).

## Quick Start

1. Build:

```bash
cd src
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

## How To Use

Run `MC_simple` in a directory containing:

- `param.def`
- `lattice.def`
- `interaction.def`

CLI forms:

```bash
./MC_simple
./MC_simple my_param.def
./MC_simple my_param.def my_lattice.def my_interaction.def
```

Optional reproducibility seed:

```bash
MC_SIMPLE_SEED=12345 ./MC_simple
```

## What It Computes

Output file: `MC_simple_result.dat`

Columns:

1. `T`
2. `E_per_site`
3. `C_per_site`
4. `M2`
5. `acceptance`

Definitions:

- `E_per_site`: thermal average of energy per site
- `C_per_site`: `C = (<E^2> - <E>^2) / (N * T^2)`
- `M2`: `<(Mx^2 + My^2 + Mz^2) / N^2>`
- `acceptance`: Metropolis acceptance ratio

## Input Files

- `param.def`: simulation parameters
- `lattice.def`: lattice size and orbital count
- `interaction.def`: bond list

### `param.def` format

Use one `key = value` pair per line:

```text
Burn_in    = 2000
Total_Step = 15000
Sample     = 4
num_temp   = 30
Ini_T      = 0.75
Delta_T    = 0.01
lambda     = 1.0
H          = 0.0
spin_dim   = 2
output_spin = 0
init_state  = 0
exchange_interval = 1
```

Supported keys:

- `Burn_in` (int): burn-in MC sweeps
- `Total_Step` (int): measurement MC sweeps
- `Sample` (int): number of independent runs
- `num_temp` (int): number of temperature slots
- `Ini_T` (double): first temperature
- `Delta_T` (double): temperature increment (`T = Ini_T + i * Delta_T`)
- `lambda` (double): z-anisotropy factor in interaction term
- `H` (double): field coupled to `S_z`
- `spin_dim` (int): `1` (Ising), `2` (XY), `3` (Heisenberg)
- `exchange_interval` (int, optional): replica exchange attempt interval
  - `0`: disabled
  - `1`: every MC sweep (default)
  - `N`: every N sweeps
- `output_spin` (int, optional): `0` off, `1` on
- `init_state` (int, optional): initial spin mode
  - `0`: random
  - `1`: FM
  - `2`: AF (Neel)
  - `3`: stripe `(pi,0)`
  - `4`: stripe `(0,pi)`

Notes:

- Lines beginning with `#` and blank lines are ignored.
- Keep the `key = value` style.
- Unknown keys are ignored with a warning.

### `lattice.def` format

Use one `key = value` pair per line:

```text
L_x     = 16
L_y     = 16
L_z     = 1
orb_num = 1
```

Supported keys:

- `L_x`, `L_y`, `L_z` (int): lattice size in each direction
- `orb_num` (int): number of orbitals per unit cell

Derived quantity:

- `All_N = L_x * L_y * L_z * orb_num`

Notes:

- Lines beginning with `#` and blank lines are ignored.
- Unknown keys are ignored with a warning.

### `interaction.def` format

Each non-comment line must contain:

```text
i_orb  j_orb  dx  dy  dz  J
```

where:

- `i_orb`, `j_orb` (int): orbital indices in `[0, orb_num-1]`
- `dx`, `dy`, `dz` (int): unit-cell displacement from `i_orb` site to `j_orb`
- `J` (double): bond coupling for that displaced pair

Example (2D square nearest-neighbor bonds):

```text
# i_orb j_orb dx dy dz J
0 0 1 0 0 1.0
0 0 0 1 0 1.0
```

Notes:

- `interaction.def` should list each unique bond once; reverse bonds are generated internally.
- Sign convention in this code is `E ~ + J (S_i . S_j)`.
- Therefore:
  - `J < 0`: ferromagnetic
  - `J > 0`: antiferromagnetic (on bipartite lattices, equivalent by sublattice spin flip at zero field)
- Malformed lines or out-of-range orbital indices are skipped with warnings.
- Neighbor counts are taken from this file to allocate bond tables (`ni_max`).

## Code Structure

Main source files are in `src/`:

- `main.c`: overall simulation flow (read input, run MC, average, write output)
- `memory.c`, `memory.h`: allocation and cleanup
- `input_parser.c`: parsing of `param.def`, `lattice.def`, `interaction.def`
- `lattice.c`: neighbor table construction, spin initialization, initial fields/energy
- `mc_update.c`: one Metropolis sweep per temperature slot
- `exmc.c`: replica exchange (Exchange MC) between adjacent temperatures
- `mc_def.h`: shared structs/prototypes
- `dSFMT.c`, `dSFMT.h`, `dSFMT-params.h`, `dSFMT-params19937.h`: random number generator
- `CMakeLists.txt`: build config

## Notes

- Single-process only (no MPI)
- Replica exchange (Exchange MC / parallel tempering) is supported via `exchange_interval`
- dSFMT parameter is fixed to `MEXP=19937` in this minimal version

## License

- This code is MIT licensed: see [`src/LICENSE`](./src/LICENSE)
- Included dSFMT files are BSD 3-Clause licensed: see [`src/LICENSE.txt`](./src/LICENSE.txt)
