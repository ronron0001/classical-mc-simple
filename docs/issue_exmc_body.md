## Summary

This issue describes the addition of **Exchange Monte Carlo (EXMC, parallel tempering)** and **overlap** measurement to `classical-mc-simple`, implemented on branch `misawa_exmc`.

## Changes

### 1. EXMC (intra-process)

- **`exchange_try_pair`**: Attempts to swap replicas at adjacent temperature slots with acceptance probability  
  `P_accept = min(1, exp((1/T_i - 1/T_j) * (E_i - E_j)))`.
- **`exchange_sweep_local`**: One sweep over adjacent pairs in odd-even order (ready for future MPI boundary exchange).
- Controlled by parameter **`use_exmc`** in `param.def` (1 = on, 0 = off; default 1).

### 2. Overlap definition

Overlap is defined as the **same-temperature** overlap between the spin configuration at **MC step 0** (start of measurement) and at **step t**:

$$Q(T) = \frac{1}{N} \sum_i \mathbf{s}_i(T, \mathrm{step}=0) \cdot \mathbf{s}_i(T, \mathrm{step}=t)$$

This is computed per temperature slot, averaged over measurement steps and samples, and written to `overlap.dat`.

### 3. Comparison with / without EXMC

- **`param_no_exmc.def`**: Same as `param.def` with `use_exmc = 0` for comparison runs.
- **`scripts/run_overlap_compare.sh`**: Runs with EXMC and without, then renames `overlap.dat` to `overlap_exmc.dat` and `overlap_no_exmc.dat`.
- **`scripts/plot_overlap.py`** and **`scripts/plot_overlap.gp`**: Produce `overlap_compare.tsv` and `overlap_compare.png`.

## Result (2D Ising, L=16, 6 temperatures)

The figure below shows the overlap with step 0 at the same temperature: **with EXMC** vs **without EXMC**.

![Overlap comparison: with EXMC vs without EXMC](https://raw.githubusercontent.com/tmisawa/classical-mc-simple/misawa_exmc/docs/overlap_compare.png)

- **Without EXMC**: At low T the replica barely decorrelates, so overlap stays close to 1; at high T it drops toward 0.
- **With EXMC**: Exchange mixes configurations across temperatures, so the configuration in each slot often comes from another T; overlap with the initial (step 0) config at that slot is small or negative, as expected.

## Files touched

- **Source**: `src/mc_def.h`, `src/mc_update.c`, `src/main.c`, `src/memory.c`, `src/memory.h`, `src/input_parser.c`
- **Params**: `samples/square_L16_Ising/param.def` (added `use_exmc`), `samples/square_L16_Ising/param_no_exmc.def`
- **Scripts**: `scripts/run_overlap_compare.sh`, `scripts/plot_overlap.py`, `scripts/plot_overlap.gp`
- **Docs**: `tips/exmc_implementation_plan.md`, `tips/calculation_conditions.md`

## References

- K. Hukushima and K. Nemoto, J. Phys. Soc. Jpn. **65**, 1604 (1996).
- `tips/exchange_mc_mpi_tips.md` (implementation notes).
