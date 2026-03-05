# Tips: Implementing Exchange MC + MPI (Beginner Guide)

This note gives practical tips for adding Exchange Monte Carlo (EXMC, parallel tempering) and MPI to `classical-mc-simple`.

The current `mc_simple` stage is:

- single-process execution
- Metropolis updates only
- no EXMC
- no MPI exchange

So the goal here is to extend that baseline step by step with a clean and testable implementation path.

## 1. EXMC in one paragraph

Run many replicas at different temperatures, update each replica normally, then attempt swaps between neighboring temperatures.

For replicas `i` and `j`, the standard acceptance is:

`P_accept = min(1, exp((1/T_i - 1/T_j) * (E_i - E_j)))`

This keeps detailed balance and improves sampling at low temperatures.

## 2. Keep the current `mc_simple` data flow

`mc_simple` already has per-temperature arrays (`sx/sy/sz`, `env_*`, `Energy`) and loops over `int_T`.
Use that structure directly when adding EXMC:

1. run one Metropolis sweep for each `int_T`
2. attempt exchanges between adjacent temperatures
3. measure observables

## 3. Temperature indexing policy for MPI

In single process, `mc_simple` uses:

`T(int_T) = Ini_T + Delta_T * int_T`

For MPI extension, use the same idea with a global index:

`gT = rank * num_temp + int_T`  
`T(rank, int_T) = Ini_T + Delta_T * gT`

This keeps indexing simple and naturally extends the current one-process setup.

## 4. Minimal staged implementation plan

1. Add intra-process EXMC first (one rank only).
2. Confirm energy/observable continuity versus temperature.
3. Add MPI process topology (`rank`, `total_proc`) and global temperature mapping.
4. Add boundary exchange between neighboring ranks (odd-even pairing to avoid deadlock).
5. Re-test acceptance and thermodynamic curves.

Implementing intra-process exchange first reduces debugging difficulty.

## 5. Common pitfalls

- temperature index offset mistakes (`rank`, `int_T`, global index)
- deadlock from symmetric send/recv ordering
- stale local fields after accepted state exchange
- using identical RNG seeds on all ranks
- temperature spacing too wide, giving near-zero exchange acceptance

## 6. References (EXMC / Parallel Tempering)

1. K. Hukushima and K. Nemoto, "Exchange Monte Carlo Method and Application to Spin Glass Simulations," *J. Phys. Soc. Jpn.* **65**, 1604-1608 (1996). DOI: https://doi.org/10.1143/jpsj.65.1604
2. D. J. Earl and M. W. Deem, "Parallel tempering: Theory, applications, and new perspectives," *Phys. Chem. Chem. Phys.* **7**, 3910-3916 (2005). DOI: https://doi.org/10.1039/B509983H
3. R. H. Swendsen and J.-S. Wang, "Replica Monte Carlo simulation of spin-glasses," *Phys. Rev. Lett.* **57**, 2607-2609 (1986). DOI: https://doi.org/10.1103/PhysRevLett.57.2607
