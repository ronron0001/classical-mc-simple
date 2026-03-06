#include "mc_def.h"
#include "memory.h"

/* Allocate a[n1][n2], zero-initialized. */
static double **alloc_double_2d_local(int n1, int n2) {
    int i;
    double **a = (double **)malloc((size_t)n1 * sizeof(double *));
    if (a == NULL)
        return NULL;
    for (i = 0; i < n1; i++) {
        a[i] = (double *)calloc((size_t)n2, sizeof(double));
        if (a[i] == NULL) {
            int j;
            for (j = 0; j < i; j++)
                free(a[j]);
            free(a);
            return NULL;
        }
    }
    return a;
}

static void free_double_2d_local(double **a, int n1) {
    int i;
    if (a == NULL)
        return;
    for (i = 0; i < n1; i++)
        free(a[i]);
    free(a);
}

/*
 * Overlap with reference configuration at one temperature slot:
 *   q(t) = (1/N) * sum_i S_i(0) · S_i(t)
 */
static double spin_overlap(struct BindStruct *X, double **ref_sx, double **ref_sy,
                           double **ref_sz, int int_T) {
    int all_i;
    int All_N = X->Def.All_N;
    int spin_dim = X->Def.spin_dim;
    double q = 0.0;

    for (all_i = 0; all_i < All_N; all_i++) {
        if (spin_dim == 1) {
            q += ref_sz[int_T][all_i] * X->Def.sz[int_T][all_i];
        } else if (spin_dim == 2) {
            q += ref_sx[int_T][all_i] * X->Def.sx[int_T][all_i] +
                 ref_sy[int_T][all_i] * X->Def.sy[int_T][all_i];
        } else {
            q += ref_sx[int_T][all_i] * X->Def.sx[int_T][all_i] +
                 ref_sy[int_T][all_i] * X->Def.sy[int_T][all_i] +
                 ref_sz[int_T][all_i] * X->Def.sz[int_T][all_i];
        }
    }

    return q / (double)All_N;
}

/*
 * Return instantaneous M^2 per site:
 *   M^2 = (Mx^2 + My^2 + Mz^2) / N^2
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
    FILE *fp_overlap;

    int num_temp, All_N, ni_max, num_pairs;
    int Burn_in, Total_Step, Sample, exchange_interval;
    int int_samp, int_T, step;
    unsigned long base_seed;
    unsigned long seed;
    double Ini_T, Delta_T;
    double **ref_sx, **ref_sy, **ref_sz;
    double *accum_overlap;
    double *sample_overlap;
    double *sample_bin_E, *sample_bin_C, *sample_bin_M2, *sample_bin_A;
    int rc;

    memset(&X, 0, sizeof(X));
    memset(&W, 0, sizeof(W));
    ref_sx = NULL;
    ref_sy = NULL;
    ref_sz = NULL;
    accum_overlap = NULL;
    sample_overlap = NULL;
    sample_bin_E = NULL;
    sample_bin_C = NULL;
    sample_bin_M2 = NULL;
    sample_bin_A = NULL;
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
    exchange_interval = X.Bind.Def.exchange_interval;
    num_pairs = (num_temp > 1) ? (num_temp - 1) : 0;

    if (exchange_interval < 0) {
        fprintf(stderr, "Error: exchange_interval must be >= 0 (0 = disabled)\n");
        goto cleanup;
    }
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
    if (allocate_work_arrays(num_temp, &W) != 0) {
        fprintf(stderr, "Error: memory allocation failed\n");
        goto cleanup;
    }

    ref_sx = alloc_double_2d_local(num_temp, All_N);
    ref_sy = alloc_double_2d_local(num_temp, All_N);
    ref_sz = alloc_double_2d_local(num_temp, All_N);
    accum_overlap =
        (double *)calloc((size_t)num_temp * (size_t)(Total_Step + 1), sizeof(double));
    sample_overlap = (double *)calloc((size_t)Sample * (size_t)num_temp *
                                          (size_t)(Total_Step + 1),
                                      sizeof(double));
    sample_bin_E = (double *)calloc((size_t)Sample * (size_t)num_temp, sizeof(double));
    sample_bin_C = (double *)calloc((size_t)Sample * (size_t)num_temp, sizeof(double));
    sample_bin_M2 = (double *)calloc((size_t)Sample * (size_t)num_temp, sizeof(double));
    sample_bin_A = (double *)calloc((size_t)Sample * (size_t)num_temp, sizeof(double));
    if (ref_sx == NULL || ref_sy == NULL || ref_sz == NULL ||
        accum_overlap == NULL || sample_overlap == NULL ||
        sample_bin_E == NULL || sample_bin_C == NULL || sample_bin_M2 == NULL ||
        sample_bin_A == NULL) {
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

    printf("MC_simple with Exchange Monte Carlo (no MPI)\n");
    printf("L=(%d,%d,%d) orb=%d spin_dim=%d All_N=%d ni_max=%d\n",
           X.Bind.Def.L_x, X.Bind.Def.L_y, X.Bind.Def.L_z, X.Bind.Def.orb_num,
           X.Bind.Def.spin_dim, All_N, ni_max);
    printf("Burn_in=%d Total_Step=%d Sample=%d num_temp=%d exchange_interval=%d\n",
           Burn_in, Total_Step, Sample, num_temp, exchange_interval);

    for (int_samp = 0; int_samp < Sample; int_samp++) {
        seed = base_seed + 8945UL * (unsigned long)int_samp;
        dsfmt_init_gen_rand(&dsfmt, seed);

        if (initial(&dsfmt, &(X.Bind), interaction_file) != 0) {
            fprintf(stderr, "Error: failed to initialize lattice from '%s'\n",
                    interaction_file);
            goto cleanup;
        }

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));

        for (step = 0; step < Burn_in; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            if (exchange_interval > 0 && num_pairs > 0 &&
                step % exchange_interval == 0) {
                exchange_step(&dsfmt, &(X.Bind), NULL, NULL);
            }
        }

        memset(X.Bind.Phys.ratio_1, 0, (size_t)num_temp * sizeof(int));
        memset(W.sample_E, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_E2, 0, (size_t)num_temp * sizeof(double));
        memset(W.sample_M2, 0, (size_t)num_temp * sizeof(double));
        if (num_pairs > 0) {
            memset(W.exchange_accept, 0, (size_t)num_pairs * sizeof(int));
            memset(W.exchange_attempt, 0, (size_t)num_pairs * sizeof(int));
        }

        for (int_T = 0; int_T < num_temp; int_T++) {
            memcpy(ref_sx[int_T], X.Bind.Def.sx[int_T], (size_t)All_N * sizeof(double));
            memcpy(ref_sy[int_T], X.Bind.Def.sy[int_T], (size_t)All_N * sizeof(double));
            memcpy(ref_sz[int_T], X.Bind.Def.sz[int_T], (size_t)All_N * sizeof(double));
            accum_overlap[(size_t)int_T * (size_t)(Total_Step + 1)] += 1.0;
            sample_overlap[(size_t)int_samp * (size_t)num_temp * (size_t)(Total_Step + 1) +
                           (size_t)int_T * (size_t)(Total_Step + 1)] = 1.0;
        }

        for (step = 0; step < Total_Step; step++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                X.Bind.Def.int_T = int_T;
                X.Bind.Phys.T = Ini_T + Delta_T * (double)int_T;
                MC(&dsfmt, &(X.Bind));
            }
            if (exchange_interval > 0 && num_pairs > 0 &&
                step % exchange_interval == 0) {
                exchange_step(&dsfmt, &(X.Bind), W.exchange_accept,
                              W.exchange_attempt);
            }
            for (int_T = 0; int_T < num_temp; int_T++) {
                double e_total = X.Bind.Phys.Energy[int_T];
                double q_val =
                    spin_overlap(&(X.Bind), ref_sx, ref_sy, ref_sz, int_T);
                W.sample_E[int_T] += e_total / (double)All_N;
                W.sample_E2[int_T] += e_total * e_total;
                W.sample_M2[int_T] += magnetization_sq(&(X.Bind), int_T);
                accum_overlap[(size_t)int_T * (size_t)(Total_Step + 1) +
                              (size_t)(step + 1)] += q_val;
                sample_overlap[(size_t)int_samp * (size_t)num_temp *
                                   (size_t)(Total_Step + 1) +
                               (size_t)int_T * (size_t)(Total_Step + 1) +
                               (size_t)(step + 1)] = q_val;
            }
        }

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

            sample_bin_E[(size_t)int_samp * (size_t)num_temp + (size_t)int_T] = e_avg;
            sample_bin_C[(size_t)int_samp * (size_t)num_temp + (size_t)int_T] = c_avg;
            sample_bin_M2[(size_t)int_samp * (size_t)num_temp + (size_t)int_T] = m2_avg;
            sample_bin_A[(size_t)int_samp * (size_t)num_temp + (size_t)int_T] = a_avg;
        }

        for (int_T = 0; int_T < num_pairs; int_T++) {
            W.accum_exchange_accept[int_T] += W.exchange_accept[int_T];
            W.accum_exchange_attempt[int_T] += W.exchange_attempt[int_T];
        }

        printf("sample %d/%d done (seed=%lu)\n", int_samp + 1, Sample, seed);
    }

    fp = fopen("MC_simple_result.dat", "w");
    if (fp != NULL) {
        fprintf(fp, "# T  E_per_site  C_per_site  M2  acceptance  overlap\n");
    } else {
        fprintf(stderr,
                "Warning: cannot open MC_simple_result.dat for write\n");
    }

    printf("\n# T  E_per_site  C_per_site  M2  acceptance  overlap\n");
    for (int_T = 0; int_T < num_temp; int_T++) {
        double T = Ini_T + Delta_T * (double)int_T;
        double E = W.accum_E[int_T] / (double)Sample;
        double C = W.accum_C[int_T] / (double)Sample;
        double M2 = W.accum_M2[int_T] / (double)Sample;
        double A = W.accum_A[int_T] / (double)Sample;
        /* overlap = <q(t)>_t = time-averaged (1/N) S(t)·S(0) over measurement phase */
        double overlap_sum = 0.0;
        for (step = 1; step <= Total_Step; step++) {
            overlap_sum += accum_overlap[(size_t)int_T * (size_t)(Total_Step + 1) +
                                         (size_t)step];
        }
        double Ov = overlap_sum / ((double)Sample * (double)Total_Step);
        printf("%.8f  %.10f  %.10f  %.10f  %.10f  %.10f\n", T, E, C, M2, A, Ov);
        if (fp != NULL) {
            fprintf(fp, "%.8f  %.10f  %.10f  %.10f  %.10f  %.10f\n",
                    T, E, C, M2, A, Ov);
        }
    }

    if (num_pairs > 0) {
        printf("\n# Exchange acceptance per pair (T_i, T_{i+1}):\n");
        if (fp != NULL) {
            fprintf(fp, "\n# Exchange acceptance per pair (T_i, T_{i+1}):\n");
        }
        for (int_T = 0; int_T < num_pairs; int_T++) {
            double T_lo = Ini_T + Delta_T * (double)int_T;
            double T_hi = Ini_T + Delta_T * (double)(int_T + 1);
            double acc = (W.accum_exchange_attempt[int_T] > 0)
                             ? (double)W.accum_exchange_accept[int_T] /
                                   (double)W.accum_exchange_attempt[int_T]
                             : 0.0;
            printf("# pair (%d,%d): T=%.4f-%.4f  accept=%d/%d  rate=%.4f\n",
                   int_T, int_T + 1, T_lo, T_hi,
                   W.accum_exchange_accept[int_T],
                   W.accum_exchange_attempt[int_T], acc);
            if (fp != NULL) {
                fprintf(fp,
                        "# pair (%d,%d): T=%.4f-%.4f  accept=%d/%d  rate=%.4f\n",
                        int_T, int_T + 1, T_lo, T_hi,
                        W.accum_exchange_accept[int_T],
                        W.accum_exchange_attempt[int_T], acc);
            }
        }
    }

    if (fp != NULL)
        fclose(fp);

    fp_overlap = fopen("MC_simple_overlap.dat", "w");
    if (fp_overlap != NULL) {
        fprintf(fp_overlap, "# step");
        for (int_T = 0; int_T < num_temp; int_T++) {
            double T = Ini_T + Delta_T * (double)int_T;
            fprintf(fp_overlap, "  q_T%.6f", T);
        }
        fprintf(fp_overlap, "\n");
        for (step = 0; step <= Total_Step; step++) {
            fprintf(fp_overlap, "%d", step);
            for (int_T = 0; int_T < num_temp; int_T++) {
                double q = accum_overlap[(size_t)int_T * (size_t)(Total_Step + 1) +
                                         (size_t)step] /
                           (double)Sample;
                fprintf(fp_overlap, "  %.10f", q);
            }
            fprintf(fp_overlap, "\n");
        }
        fclose(fp_overlap);
    } else {
        fprintf(stderr, "Warning: cannot open MC_simple_overlap.dat for write\n");
    }

    fp_overlap = fopen("MC_simple_overlap_samples.dat", "w");
    if (fp_overlap != NULL) {
        fprintf(fp_overlap, "# sample  step");
        for (int_T = 0; int_T < num_temp; int_T++) {
            double T = Ini_T + Delta_T * (double)int_T;
            fprintf(fp_overlap, "  q_T%.6f", T);
        }
        fprintf(fp_overlap, "\n");
        for (int_samp = 0; int_samp < Sample; int_samp++) {
            for (step = 0; step <= Total_Step; step++) {
                fprintf(fp_overlap, "%d  %d", int_samp, step);
                for (int_T = 0; int_T < num_temp; int_T++) {
                    size_t idx = (size_t)int_samp * (size_t)num_temp *
                                     (size_t)(Total_Step + 1) +
                                 (size_t)int_T * (size_t)(Total_Step + 1) +
                                 (size_t)step;
                    fprintf(fp_overlap, "  %.10f", sample_overlap[idx]);
                }
                fprintf(fp_overlap, "\n");
            }
        }
        fclose(fp_overlap);
    } else {
        fprintf(stderr,
                "Warning: cannot open MC_simple_overlap_samples.dat for write\n");
    }

    fp_overlap = fopen("MC_simple_sample_bins.dat", "w");
    if (fp_overlap != NULL) {
        fprintf(fp_overlap, "# sample  T  E_per_site  C_per_site  M2  acceptance\n");
        for (int_samp = 0; int_samp < Sample; int_samp++) {
            for (int_T = 0; int_T < num_temp; int_T++) {
                double T = Ini_T + Delta_T * (double)int_T;
                size_t idx = (size_t)int_samp * (size_t)num_temp + (size_t)int_T;
                fprintf(fp_overlap, "%d  %.8f  %.10f  %.10f  %.10f  %.10f\n",
                        int_samp, T, sample_bin_E[idx], sample_bin_C[idx],
                        sample_bin_M2[idx], sample_bin_A[idx]);
            }
        }
        fclose(fp_overlap);
    } else {
        fprintf(stderr, "Warning: cannot open MC_simple_sample_bins.dat for write\n");
    }

    rc = 0;

cleanup:
    free(sample_bin_E);
    free(sample_bin_C);
    free(sample_bin_M2);
    free(sample_bin_A);
    free(sample_overlap);
    free(accum_overlap);
    free_double_2d_local(ref_sx, num_temp);
    free_double_2d_local(ref_sy, num_temp);
    free_double_2d_local(ref_sz, num_temp);
    free_work_arrays(&W);
    free_minimal_mc_arrays(&X);
    return rc;
}
