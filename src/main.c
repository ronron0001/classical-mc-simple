#include "mc_def.h"
#include "memory.h"

/*
 * Overlap at temperature slot k: (1/N) sum_i s_i(step=0) · s_i(step=t).
 * Same temperature, initial config vs current config.
 */
static double overlap_with_initial(struct BindStruct *X, struct SimpleWorkArrays *W, int k) {
    int all_i, All_N;
    double sum;

    All_N = X->Def.All_N;
    sum = 0.0;
    for (all_i = 0; all_i < All_N; all_i++) {
        sum += W->sx0[k][all_i] * X->Def.sx[k][all_i] +
               W->sy0[k][all_i] * X->Def.sy[k][all_i] +
               W->sz0[k][all_i] * X->Def.sz[k][all_i];
    }
    return sum / (double)All_N;
}

/*
 * Return instantaneous M^2 per site:
 *   M^2 = (Mx^2 + My^2 + Mz^2) / N^2
 * This is computed at each MC step, then averaged in time.
 */
static double magnetization_sq(struct BindStruct *X, int int_T) {
    int all_i, All_N;
    double mx, my, mz;

    mx = 0.0;
    my = 0.0;
    mz = 0.0;
    All_N = X->Def.All_N;

    for (all_i = 0; all_i < All_N; all_i++) {
        mx += X->Def.sx[int_T][all_i];
        my += X->Def.sy[int_T][all_i];
        mz += X->Def.sz[int_T][all_i];
    }

    return (mx * mx + my * my + mz * mz) / ((double)All_N * (double)All_N);
}

/*
 * Accepted CLI forms:
 *   1) ./MC_simple
 *   2) ./MC_simple param.def
 *   3) ./MC_simple param.def lattice.def interaction.def
 */
static int parse_input_files(int argc, char **argv, const char **param_file,
                             const char **lattice_file,
                             const char **interaction_file) {
    *param_file = "param.def";
    *lattice_file = "lattice.def";
    *interaction_file = "interaction.def";

    if (argc == 2) {
        *param_file = argv[1];
        return 0;
    }

    if (argc == 4) {
        *param_file = argv[1];
        *lattice_file = argv[2];
        *interaction_file = argv[3];
        return 0;
    }

    if (argc == 1)
        return 0;

    return -1;
}

int main(int argc, char **argv) {
    struct MCMainCalStruct X;
    struct SimpleWorkArrays W;
    dsfmt_t dsfmt;
    const char *param_file;
    const char *lattice_file;
    const char *interaction_file;
    const char *seed_env;
    char *endptr;
    FILE *fp;

    int num_temp, All_N, ni_max;
    int Burn_in, Total_Step, Sample;
    int int_samp, int_T, step, k, all_i;
    int n_ex_accept, total_ex_accept, total_ex_attempt;
    unsigned long base_seed;
    unsigned long seed;
    double Ini_T, Delta_T;
    int rc;

    memset(&X, 0, sizeof(X));
    memset(&W, 0, sizeof(W));
    rc = 1;

    if (parse_input_files(argc, argv, &param_file, &lattice_file,
                          &interaction_file) != 0) {
        fprintf(stderr,
                "Usage: %s [param.def] or %s [param.def lattice.def "
                "interaction.def]\n",
                argv[0], argv[0]);
        return 1;
    }

    X.Bind.Def.myrank = MASTER;
    X.Bind.Def.total_proc = 1;

    /* Read model/simulation settings from input files. */
    if (read_param(param_file, &X.Bind.Def) != 0)
        return 1;
    if (read_lattice(lattice_file, &X.Bind.Def) != 0)
        return 1;
    if (read_interaction(interaction_file, &X.Bind.Def) < 0)
        return 1;

    num_temp = X.Bind.Def.num_temp;
    All_N = X.Bind.Def.All_N;
    ni_max = X.Bind.Def.ni_max;
    Burn_in = X.Bind.Def.Burn_in;
    Total_Step = X.Bind.Def.Total_Step;
    Sample = X.Bind.Def.Sample;
    Ini_T = X.Bind.Def.Ini_T;
    Delta_T = X.Bind.Def.Delta_T;
    int use_exmc = X.Bind.Def.use_exmc;

    if (num_temp <= 0 || All_N <= 0 || ni_max <= 0) {
        fprintf(stderr, "Error: invalid lattice/temperature setting\n");
        goto cleanup;
    }
    if (Burn_in < 0 || Total_Step <= 0 || Sample <= 0) {
        fprintf(
            stderr,
            "Error: Burn_in >= 0, Total_Step > 0, Sample > 0 are required\n");
        goto cleanup;
    }

    if (allocate_minimal_mc_arrays(&X) != 0) {
        fprintf(stderr, "Error: memory allocation failed\n");
        goto cleanup;
    }

    if (allocate_work_arrays(num_temp, All_N, &W) != 0) {
        fprintf(stderr, "Error: memory allocation failed\n");
        goto cleanup;
    }

    base_seed = 11456UL;
    seed_env = getenv("MC_SIMPLE_SEED");
    if (seed_env != NULL && seed_env[0] != '\0') {
        base_seed = strtoul(seed_env, &endptr, 10);
        if (endptr == seed_env || *endptr != '\0') {
            fprintf(stderr,
                    "Warning: invalid MC_SIMPLE_SEED='%s'; using default %lu\n",
                    seed_env, base_seed);
            base_seed = 11456UL;
        }
    }

    printf("MC_simple %s (intra-process exchange, no MPI)\n",
           use_exmc ? "with EXMC" : "without EXMC");
    printf("L=(%d,%d,%d) orb=%d spin_dim=%d All_N=%d ni_max=%d\n",
           X.Bind.Def.L_x, X.Bind.Def.L_y, X.Bind.Def.L_z, X.Bind.Def.orb_num,
           X.Bind.Def.spin_dim, All_N, ni_max);
    printf("Burn_in=%d Total_Step=%d Sample=%d num_temp=%d\n", Burn_in,
           Total_Step, Sample, num_temp);

    total_ex_accept = 0;
    total_ex_attempt = 0;

    /*
     * Outer loop over independent samples (different RNG seeds).
     * We average observables over these samples at the end.
     */
    for (int_samp = 0; int_samp < Sample; int_samp++) {
        seed = base_seed + 8945UL * (unsigned long)int_samp;
        dsfmt_init_gen_rand(&dsfmt, seed);

        /* Initialize spins and neighbor-based effective fields. */
        initial(&dsfmt, &(X.Bind));

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));

        /* Burn-in: evolve the system without collecting measurements. */
        for (step = 0; step < Burn_in; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            if (use_exmc)
                exchange_sweep_local(&dsfmt, &(X.Bind));
        }

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));
        memset(W.sample_E, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_E2, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_M2, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_overlap, 0, (size_t)num_temp * sizeof(double));

        /* Save spin config at step 0 (start of measurement) for overlap. */
        for (int_T = 0; int_T < num_temp; int_T++) {
            for (all_i = 0; all_i < All_N; all_i++) {
                W.sx0[int_T][all_i] = X.Bind.Def.sx[int_T][all_i];
                W.sy0[int_T][all_i] = X.Bind.Def.sy[int_T][all_i];
                W.sz0[int_T][all_i] = X.Bind.Def.sz[int_T][all_i];
            }
        }

        /* Measurement phase: MC sweep -> (EXMC exchange if use_exmc) -> record. */
        for (step = 0; step < Total_Step; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            if (use_exmc) {
                n_ex_accept = exchange_sweep_local(&dsfmt, &(X.Bind));
                total_ex_accept += n_ex_accept;
                total_ex_attempt += (num_temp > 1) ? (num_temp - 1) : 0;
            }

            for (int_T = 0; int_T < num_temp; int_T++) {
                double e_total = X.Bind.Phys.Energy[int_T];
                W.sample_E[int_T] += e_total / (double)All_N;
                W.sample_E2[int_T] += e_total * e_total;
                W.sample_M2[int_T] += magnetization_sq(&(X.Bind), int_T);
            }
            for (k = 0; k < num_temp; k++)
                W.sample_overlap[k] += overlap_with_initial(&(X.Bind), &W, k);
        }

        /*
         * Convert step sums to per-sample averages.
         * Specific heat per site:
         *   C = ( <E^2> - <E>^2 ) / (N * T^2)
         */
        for (int_T = 0; int_T < num_temp; int_T++) {
            double T_cur = Ini_T + Delta_T * (double)int_T;
            double e_avg = W.sample_E[int_T] / (double)Total_Step;
            double e2_avg = W.sample_E2[int_T] / (double)Total_Step;
            double m2_avg = W.sample_M2[int_T] / (double)Total_Step;
            double a_avg = (double)X.Bind.Phys.ratio_1[int_T] /
                           ((double)All_N * (double)Total_Step);
            double e_total_avg = e_avg * (double)All_N;
            double c_avg = 0.0;
            if (T_cur > 0.0) {
                c_avg = (e2_avg - e_total_avg * e_total_avg) /
                        ((double)All_N * T_cur * T_cur);
            }
            W.accum_E[int_T] += e_avg;
            W.accum_C[int_T] += c_avg;
            W.accum_M2[int_T] += m2_avg;
            W.accum_A[int_T] += a_avg;
        }
        for (k = 0; k < num_temp; k++)
            W.accum_overlap[k] += W.sample_overlap[k] / (double)Total_Step;
        printf("sample %d/%d done (seed=%lu)\n", int_samp + 1, Sample, seed);
    }

    fp = fopen("MC_simple_result.dat", "w");
    if (fp != NULL) {
        fprintf(fp, "# T  E_per_site  C_per_site  M2  acceptance\n");
    } else {
        fprintf(stderr,
                "Warning: cannot open MC_simple_result.dat for write\n");
    }

    if (total_ex_attempt > 0) {
        printf("exchange_acceptance = %.6f (%d / %d)\n",
               (double)total_ex_accept / (double)total_ex_attempt,
               total_ex_accept, total_ex_attempt);
        if (fp != NULL)
            fprintf(fp, "# exchange_acceptance = %.6f (%d / %d)\n",
                    (double)total_ex_accept / (double)total_ex_attempt,
                    total_ex_accept, total_ex_attempt);
    }

    printf("\n# T  E_per_site  C_per_site  M2  acceptance\n");
    /* Final average over independent samples. */
    for (int_T = 0; int_T < num_temp; int_T++) {
        double T = Ini_T + Delta_T * (double)int_T;
        double E = W.accum_E[int_T] / (double)Sample;
        double C = W.accum_C[int_T] / (double)Sample;
        double M2 = W.accum_M2[int_T] / (double)Sample;
        double A = W.accum_A[int_T] / (double)Sample;
        printf("%.8f  %.10f  %.10f  %.10f  %.10f\n", T, E, C, M2, A);
        if (fp != NULL) {
            fprintf(fp, "%.8f  %.10f  %.10f  %.10f  %.10f\n", T, E, C, M2, A);
        }
    }

    if (fp != NULL)
        fclose(fp);

    /* Write overlap with initial (step 0) at same temperature: T, overlap. */
    {
        FILE *fp_ov = fopen("overlap.dat", "w");
        if (fp_ov != NULL) {
            fprintf(fp_ov, "# use_exmc = %d\n", use_exmc);
            fprintf(fp_ov, "# Overlap = (1/N) sum_i s_i(step=0) . s_i(step=t), same T\n");
            fprintf(fp_ov, "# int_T   T   overlap\n");
            for (k = 0; k < num_temp; k++) {
                double T = Ini_T + Delta_T * (double)k;
                double ov = W.accum_overlap[k] / (double)Sample;
                fprintf(fp_ov, "%d  %.6f  %.10f\n", k, T, ov);
            }
            fclose(fp_ov);
        }
    }

    rc = 0;

cleanup:
    free_work_arrays(num_temp, All_N, &W);
    free_minimal_mc_arrays(&X);
    return rc;
}
