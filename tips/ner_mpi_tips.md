# Tips: Implementing NER + MPI (Beginner Guide)

This note explains a practical path to add Non-Equilibrium Relaxation (NER) and MPI parallelization to `classical-mc-simple`.

## 1. What NER is

In NER, you do **not** equilibrate first.  
You prepare a controlled initial state (ordered or random), quench to a target temperature, and track time-dependent observables as a function of MC step `t`.

Typical examples:

- magnetization decay from ordered initial state
- autocorrelation or overlap with the initial state
- fluctuation growth from random initial state

Near criticality, many observables show power-law-like time dependence before finite-size saturation.

## 2. Current status in `mc_simple`

Useful NER-related pieces already exist in `src`:

- `init_state` is parsed in `param.def`
- `initial_ner(...)` exists in `lattice.c`
- struct fields for initial-spin overlap and NER diagnostics exist in `mc_def.h`

But the current executable path is still equilibrium-only:

- `main.c` calls `initial(...)` (not `initial_ner(...)`)
- no dedicated NER time-series output loop
- minimal allocator (`memory.c`) does not allocate NER-specific arrays such as `Ini_sx/Ini_sy/Ini_sz`

So NER support is currently partial scaffolding, not a complete runtime mode.

## 3. Minimal NER implementation plan

1. Add a runtime switch in `param.def` (for example `run_mode = ner` or `enable_ner = 1`).
2. In `main.c`, branch initialization:
   - equilibrium mode: keep `initial(...)`
   - NER mode: call `initial_ner(...)`
3. In `memory.c`, allocate/free NER arrays needed by `initial_ner(...)`:
   - `Ini_sx`, `Ini_sy`, `Ini_sz` (size `All_N`)
4. For NER mode, use a single target temperature first (`num_temp = 1`) to simplify validation.
5. Replace equilibrium averaging loop with a time-series loop:
   - at each step, run `MC(...)`
   - measure and write `m(t)`, `m2(t)`, energy, and optional overlap `q(t)`
6. Write one line per MC step to `NER_result.dat`, e.g.:
   - `t  T  m  m2  e_per_site  q`

This is the smallest path that produces a usable NER dataset.

## 4. Recommended NER observables

For a first implementation, keep observables simple and robust:

- `m(t)`: order parameter magnitude (or component for Ising)
- `m2(t)`: second moment of order parameter
- `e(t)`: energy per site
- `q(t) = (1/N) sum_i S_i(t) · S_i(0)` (overlap with initial configuration)

The overlap is especially useful because `Ini_s*` is already part of the data model.

## 5. MPI strategy for beginners

For NER, the easiest MPI parallelization is **independent-run parallelism** (not domain decomposition):

1. Every rank simulates the full lattice at the same temperature.
2. Each rank uses a different seed.
3. All ranks run the same number of MC steps.
4. At each output step (or at the end), combine observables with `MPI_Reduce`/`MPI_Allreduce`.

Why this is best first:

- very small code changes
- no boundary communication of spins
- nearly ideal parallel efficiency
- straightforward statistical error estimation across ranks

## 6. MPI seeding and reproducibility

Do not reuse the same seed on all ranks.

A safe pattern is:

- `seed = base_seed + stride * rank + sample_offset`

Choose `stride` large enough (for example `1000003`) and print `(rank, seed)` once at startup.

## 7. Error bars for NER curves

A practical first method:

1. Each rank outputs local time series.
2. Compute rank-wise mean at each `t`.
3. Estimate error bar as standard error over ranks:
   - `stderr(t) = stddev_rank(x_t) / sqrt(n_rank)`

After this works, move to binning/jackknife if needed.

## 8. Common pitfalls

- mixing equilibrium burn-in logic into NER mode
- not allocating `Ini_s*` before `initial_ner(...)`
- combining MPI data with inconsistent step counters
- comparing different ranks with different temperatures by mistake
- interpreting finite-size saturation region as critical power law

## 9. Representative NER references

1. N. Ito and Y. Ozeki, "Nonequilibrium Relaxation Method And Its Applications," *Int. J. Mod. Phys. C* **10**, 1495-1502 (1999). DOI: https://doi.org/10.1142/S0129183199001273
2. N. Ito, K. Hukushima, K. Ogawa, and Y. Ozeki, "Nonequilibrium Relaxation of Fluctuations of Physical Quantities," *J. Phys. Soc. Jpn.* **69**, 1931-1934 (2000). DOI: https://doi.org/10.1143/jpsj.69.1931
3. Y. Ozeki and N. Ito, "Nonequilibrium relaxation study of Ising spin glass models," *Phys. Rev. B* **64**, 024416 (2001). DOI: https://doi.org/10.1103/PhysRevB.64.024416
4. Y. Ozeki, K. Ogawa, and N. Ito, "Nonequilibrium relaxation analysis of Kosterlitz-Thouless phase transition," *Phys. Rev. E* **67**, 026702 (2003). DOI: https://doi.org/10.1103/PhysRevE.67.026702

These references are a good starting set for both method basics and concrete analysis patterns.
