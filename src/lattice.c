#include "mc_def.h"

/*
 * lattice.c (stage1 minimal runtime)
 *
 * Responsibilities:
 *   1) Build neighbor table (Nr) and bond strengths (J) from interaction.def.
 *   2) Initialize spin configurations.
 *   3) Compute initial local fields (env_*) and total energy.
 *
 * Important convention:
 *   - all_i is a flattened site index including orbital:
 *       all_i = orb + site_index * orb_num
 *   - site_index itself is flattened from (ix, iy, iz).
 */
int build_lattice(const char *filename, struct BindStruct *X) {
    FILE *fp;
    char line[1024];
    int i_orb, j_orb, dx, dy, dz;
    double J_val;
    int ix, iy, iz, all_i, all_j, jx, jy, jz, site_i, site_j;
    int L_x, L_y, L_z, orb_num, All_N;
    /*
     * count[all_i] tracks how many neighbor slots are already filled
     * for each site. This allows us to append bonds directly into Nr/J.
     */
    int *count;

    L_x = X->Def.L_x;
    L_y = X->Def.L_y;
    L_z = X->Def.L_z;
    orb_num = X->Def.orb_num;
    All_N = X->Def.All_N;

    count = (int *)calloc(All_N, sizeof(int));

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open %s\n", filename);
        free(count);
        return -1;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        /* skip comments and blank lines */
        char *p = line;
        while (*p == ' ' || *p == '\t')
            p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0')
            continue;

        if (sscanf(line, "%d %d %d %d %d %lf", &i_orb, &j_orb, &dx, &dy, &dz,
                   &J_val) != 6)
            continue;

        /*
         * One interaction.def row defines a bond pattern:
         *   (i_orb in cell R) --J--> (j_orb in cell R + (dx,dy,dz)).
         * We replicate this pattern over all unit cells.
         *
         * We insert BOTH directions (i->j and j->i), so each site has
         * a complete neighbor list for fast local-field updates.
         */
        for (iz = 0; iz < L_z; iz++) {
            for (iy = 0; iy < L_y; iy++) {
                for (ix = 0; ix < L_x; ix++) {
                    /* Source site in flattened indexing. */
                    site_i = ix + iy * L_x + iz * L_x * L_y;
                    all_i = i_orb + site_i * orb_num;

                    /*
                     * Target cell with periodic boundary conditions:
                     * wrap each coordinate into [0, L_dim-1].
                     */
                    jx = (ix + dx % L_x + L_x) % L_x;
                    jy = (iy + dy % L_y + L_y) % L_y;
                    jz = (iz + dz % L_z + L_z) % L_z;
                    site_j = jx + jy * L_x + jz * L_x * L_y;
                    all_j = j_orb + site_j * orb_num;

                    /* Forward: i -> j */
                    X->Def.Nr[count[all_i]][all_i] = all_j;
                    X->Def.J[count[all_i]][all_i] = J_val;
                    count[all_i]++;

                    /* Reverse: j -> i */
                    X->Def.Nr[count[all_j]][all_j] = all_i;
                    X->Def.J[count[all_j]][all_j] = J_val;
                    count[all_j]++;
                }
            }
        }
    }
    fclose(fp);
    free(count);

    /*
     * Precompute sign factors used by order parameters.
     * These are inexpensive to store and avoid recomputing parity every step.
     */
    for (all_i = 0; all_i < All_N; all_i++) {
        int site_i_l = all_i / orb_num;
        int iz_l = site_i_l / (L_x * L_y);
        int iy_l = (site_i_l % (L_x * L_y)) / L_x;
        int ix_l = site_i_l % L_x;
        X->Def.stag_sign[all_i] = ((ix_l + iy_l + iz_l) % 2 == 0) ? 1 : -1;
        X->Def.stripe_sign_1[all_i] = (ix_l % 2 == 0) ? 1 : -1;
        X->Def.stripe_sign_2[all_i] = (iy_l % 2 == 0) ? 1 : -1;
    }

    return 0;
}

void initial_ner(dsfmt_t *dsfmt, struct BindStruct *X) {
    int all_i, all_j, n_i;
    int All_N, ni_max, spin_dim, init_state;
    double tmp_E, theta, x, y, z, tmp, tmp_2;
    double lambda;

    All_N = X->Def.All_N;
    ni_max = X->Def.ni_max;
    lambda = X->Def.lambda;
    spin_dim = X->Def.spin_dim;
    init_state = X->Def.init_state;

    /*
     * NER uses exactly one temperature slot (int_T=0).
     * We initialize from a chosen ordered/random state, then monitor
     * relaxation as MC steps proceed.
     */
    for (all_i = 0; all_i < All_N; all_i++) {
        if (init_state == 1) {
            /* FM: all spins aligned */
            if (spin_dim == 2) {
                x = 1.0;
                y = 0.0;
                z = 0.0; /* XY: along +x */
            } else {
                x = 0.0;
                y = 0.0;
                z = 1.0; /* Ising/Heisenberg: along +z */
            }
        } else if (init_state == 2) {
            /* AF (Neel): staggered sign (-1)^(ix+iy+iz) */
            if (spin_dim == 2) {
                x = (double)X->Def.stag_sign[all_i];
                y = 0.0;
                z = 0.0;
            } else {
                x = 0.0;
                y = 0.0;
                z = (double)X->Def.stag_sign[all_i];
            }
        } else if (init_state == 3) {
            /* Stripe1 (pi,0): sign (-1)^ix */
            if (spin_dim == 2) {
                x = (double)X->Def.stripe_sign_1[all_i];
                y = 0.0;
                z = 0.0;
            } else {
                x = 0.0;
                y = 0.0;
                z = (double)X->Def.stripe_sign_1[all_i];
            }
        } else if (init_state == 4) {
            /* Stripe2 (0,pi): sign (-1)^iy */
            if (spin_dim == 2) {
                x = (double)X->Def.stripe_sign_2[all_i];
                y = 0.0;
                z = 0.0;
            } else {
                x = 0.0;
                y = 0.0;
                z = (double)X->Def.stripe_sign_2[all_i];
            }
        } else {
            /* RANDOM */
            if (spin_dim == 3) {
                /*
                 * Heisenberg random unit vector on S^2 via Marsaglia method:
                 * sample (x,y) uniformly in unit disk, then map to sphere.
                 */
                x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                tmp = x * x + y * y;
                while (tmp > 1.0) {
                    x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                    y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                    tmp = x * x + y * y;
                }
                tmp_2 = sqrt(1 - tmp);
                x = 2 * x * tmp_2;
                y = 2 * y * tmp_2;
                z = 1 - 2 * tmp;
            } else if (spin_dim == 2) {
                /* XY random unit vector on circle S^1. */
                theta = 2.0 * PI * dsfmt_genrand_close_open(dsfmt);
                x = cos(theta);
                y = sin(theta);
                z = 0.0;
            } else {
                /* Ising random z = +-1. */
                x = 0.0;
                y = 0.0;
                z = (dsfmt_genrand_close_open(dsfmt) < 0.5) ? 1.0 : -1.0;
            }
        }
        X->Def.sx[0][all_i] = x;
        X->Def.sy[0][all_i] = y;
        X->Def.sz[0][all_i] = z;
    }

    /* Save initial configuration */
    for (all_i = 0; all_i < All_N; all_i++) {
        X->Def.Ini_sx[all_i] = X->Def.sx[0][all_i];
        X->Def.Ini_sy[all_i] = X->Def.sy[0][all_i];
        X->Def.Ini_sz[all_i] = X->Def.sz[0][all_i];
    }

    /*
     * Build local effective field env_* for each site:
     *   env = sum_j J_ij * S_j  (+ anisotropy on z, + external field H).
     * MC update will reuse env to evaluate local energy change quickly.
     */
    for (all_i = 0; all_i < All_N; all_i++) {
        X->Def.env_sx[0][all_i] = 0.0;
        X->Def.env_sy[0][all_i] = 0.0;
        X->Def.env_sz[0][all_i] = 0.0;
        for (n_i = 0; n_i < ni_max; n_i++) {
            all_j = X->Def.Nr[n_i][all_i];
            X->Def.env_sx[0][all_i] +=
                X->Def.J[n_i][all_i] * X->Def.sx[0][all_j];
            X->Def.env_sy[0][all_i] +=
                X->Def.J[n_i][all_i] * X->Def.sy[0][all_j];
            X->Def.env_sz[0][all_i] +=
                lambda * X->Def.J[n_i][all_i] * X->Def.sz[0][all_j];
        }
        X->Def.env_sz[0][all_i] += X->Def.H;
    }

    /*
     * Compute total energy from pair interactions and field term.
     * The neighbor table contains both i->j and j->i, so divide by 2.
     */
    X->Phys.Energy[0] = 0.0;
    for (all_i = 0; all_i < All_N; all_i++) {
        tmp_E = 0.0;
        for (n_i = 0; n_i < ni_max; n_i++) {
            all_j = X->Def.Nr[n_i][all_i];
            tmp_E += X->Def.J[n_i][all_i] *
                     (X->Def.sx[0][all_i] * X->Def.sx[0][all_j] +
                      X->Def.sy[0][all_i] * X->Def.sy[0][all_j] +
                      lambda * X->Def.sz[0][all_i] * X->Def.sz[0][all_j]);
        }
        tmp_E += 2 * X->Def.H * X->Def.sz[0][all_i];
        X->Phys.Energy[0] += tmp_E;
    }
    X->Phys.Energy[0] = X->Phys.Energy[0] * 0.5;
}

int initial(dsfmt_t *dsfmt, struct BindStruct *X, const char *interaction_file) {
    int all_i, all_j, n_i;
    int int_T;
    int All_N, ni_max, spin_dim;
    double tmp, tmp_2, theta, x, y, z;
    double tmp_E;
    double lambda;

    int num_temp = X->Def.num_temp;

    All_N = X->Def.All_N;
    ni_max = X->Def.ni_max;
    lambda = X->Def.lambda;
    spin_dim = X->Def.spin_dim;

    /* Build neighbor list from the interaction file selected by CLI/default. */
    if (interaction_file == NULL) {
        fprintf(stderr, "Error: interaction file path is NULL\n");
        return -1;
    }
    if (build_lattice(interaction_file, X) != 0) {
        return -1;
    }

    /*
     * Initialize each temperature replica independently with random spins.
     * (Replica exchange is removed in stage1, but the data layout keeps
     * the original multi-temperature structure.)
     */
    for (int_T = 0; int_T < num_temp; int_T++) {
        for (all_i = 0; all_i < All_N; all_i++) {
            if (spin_dim == 3) {
                /* Heisenberg: Marsaglia method for uniform S^2 sampling */
                x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                tmp = x * x + y * y;
                while (tmp > 1.0) {
                    x = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                    y = 2 * (dsfmt_genrand_close_open(dsfmt) - 0.5);
                    tmp = x * x + y * y;
                }
                tmp_2 = sqrt(1 - tmp);
                x = 2 * x * tmp_2;
                y = 2 * y * tmp_2;
                z = 1 - 2 * tmp;
            } else if (spin_dim == 2) {
                /* XY: uniform angle on S^1 */
                theta = 2.0 * PI * dsfmt_genrand_close_open(dsfmt);
                x = cos(theta);
                y = sin(theta);
                z = 0.0;
            } else {
                /* Ising: random +1 or -1 */
                x = 0.0;
                y = 0.0;
                z = (dsfmt_genrand_close_open(dsfmt) < 0.5) ? 1.0 : -1.0;
            }

            X->Def.sx[int_T][all_i] = x;
            X->Def.sy[int_T][all_i] = y;
            X->Def.sz[int_T][all_i] = z;
        }

        /*
         * Compute local effective field for this temperature slot.
         * This allows O(1) proposal evaluation per spin in MC().
         */
        for (all_i = 0; all_i < All_N; all_i++) {
            X->Def.env_sx[int_T][all_i] = 0.0;
            X->Def.env_sy[int_T][all_i] = 0.0;
            X->Def.env_sz[int_T][all_i] = 0.0;
            for (n_i = 0; n_i < ni_max; n_i++) {
                all_j = X->Def.Nr[n_i][all_i];
                X->Def.env_sx[int_T][all_i] +=
                    X->Def.J[n_i][all_i] * X->Def.sx[int_T][all_j];
                X->Def.env_sy[int_T][all_i] +=
                    X->Def.J[n_i][all_i] * X->Def.sy[int_T][all_j];
                X->Def.env_sz[int_T][all_i] +=
                    lambda * X->Def.J[n_i][all_i] * X->Def.sz[int_T][all_j];
            }
            X->Def.env_sz[int_T][all_i] += X->Def.H;
        }

        /*
         * Compute initial total energy for this temperature.
         * Factor 1/2 removes double counting from bidirectional bonds.
         */
        X->Phys.Energy[int_T] = 0.0;
        for (all_i = 0; all_i < All_N; all_i++) {
            tmp_E = 0.0;
            for (n_i = 0; n_i < ni_max; n_i++) {
                all_j = X->Def.Nr[n_i][all_i];
                tmp_E += X->Def.J[n_i][all_i] *
                         (X->Def.sx[int_T][all_i] * X->Def.sx[int_T][all_j] +
                          X->Def.sy[int_T][all_i] * X->Def.sy[int_T][all_j] +
                          lambda * X->Def.sz[int_T][all_i] *
                              X->Def.sz[int_T][all_j]);
            }
            tmp_E += 2 * X->Def.H * X->Def.sz[int_T][all_i];
            X->Phys.Energy[int_T] += tmp_E;
        }
        X->Phys.Energy[int_T] = (X->Phys.Energy[int_T]) * 0.5;
    }

    return 0;
}
