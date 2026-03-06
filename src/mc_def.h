#ifndef MC_DEF_H
#define MC_DEF_H

#include "dSFMT.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#ifndef NO_MPI
#include "mpi.h"
#endif

/*=== Constants ===*/
#define MASTER 0
#define PI 3.141592654

/*=== Struct definitions ===*/
struct DefineList {
    /* Lattice geometry */
    int L_x, L_y, L_z;
    int orb_num;
    int spin_dim; /* 1=Ising, 2=XY, 3=Heisenberg */
    int All_N, ni_max;

    /* Simulation parameters */
    int Burn_in, Total_Step, Sample, num_temp;
    int use_exmc; /* 1=enable EXMC (default), 0=disable */
    double Ini_T, Delta_T;

    /* Hamiltonian parameters */
    double lambda;
    double H;

    /* Output control */
    int output_spin; /* 0: skip spin output (default), 1: enable */

    /* NER (Non-Equilibrium Relaxation) */
    int init_state; /* 0=RANDOM(default), 1=FM, 2=AF(Neel), 3=Stripe1(pi,0),
                       4=Stripe2(0,pi) */
    double *Ini_sx, *Ini_sy, *Ini_sz; /* initial spin config [All_N] */

    /* Internal state */
    int int_T, myrank, total_proc;
    int *Mloop, *trace;

    double *send_buffer, *recv_buffer;
    double **sx, **sy, **sz;
    double *All_sx, *All_sy, *All_sz;
    double *AG_sx, *AG_sy, *AG_sz;
    double **env_sx, **env_sy, **env_sz;
    double **dist_loop, *All_dist_loop, *AG_dist_loop;
    double **Ave_dist_loop, **Err_dist_loop;
    double *All_Ave_dist_loop, *AG_Ave_dist_loop;
    double *All_Err_dist_loop, *AG_Err_dist_loop;

    int *stag_sign;     /* (-1)^(ix+iy+iz) for staggered magnetization */
    int *stripe_sign_1; /* (-1)^ix for stripe order (pi,0) */
    int *stripe_sign_2; /* (-1)^iy for stripe order (0,pi) */

    int **Nr;
    double **J;
};

struct PhysList {
    int step;
    int *ratio_1;

    double tmp_E, tmp_Mx, tmp_My, tmp_Mz;
    double T;
    double *Energy;

    double *Sum_Ene, *Sum_Ene_2, *Sum_Ene_4;
    double *Sum_Spc;

    double *Ave_Ene, *Err_Ene;
    double *All_Ene, *All_Err_Ene;
    double *Ave_Ene_2, *Err_Ene_2;
    double *Ave_EB, *Err_EB;
    double *Ave_Spc, *Err_Spc;
    double *All_Spc, *All_Err_Spc;
    double *All_Ene_4, *All_Ene_2, *All_M_4;
    double *All_EB, *All_MB;

    double **M, **Quad, **Octa;
    double **Sum_M, **Sum_Quad, **Sum_Octa;
    double **Sum_M_2, **Sum_Quad_2, **Sum_Octa_2;
    double **Sum_M_4;

    double **Ave_M, **Ave_Quad, **Ave_Octa;
    double **Err_M, **Err_Quad, **Err_Octa;

    double **Ave_M_2, **Ave_Quad_2, **Ave_Octa_2;
    double **Err_M_2, **Err_Quad_2, **Err_Octa_2;
    double **Ave_MB, **Err_MB;

    double *All_M, *All_Quad, *All_Octa;
    double *All_Err_M, *All_Err_Quad, *All_Err_Octa;
    double *All_Err_EB, *All_Err_MB;

    double *All_M_2, *All_Quad_2, *All_Octa_2;
    double *All_Err_M_2, *All_Err_Quad_2, *All_Err_Octa_2;

    /* Staggered magnetization */
    double **Ms;
    double **Sum_Ms, **Sum_Ms_2, **Sum_Ms_4;
    double **Ave_Ms, **Ave_Ms_2, **Ave_MsB;
    double **Err_Ms, **Err_Ms_2, **Err_MsB;
    double *All_Ms, *All_Ms_2, *All_MsB;
    double *All_Err_Ms, *All_Err_Ms_2, *All_Err_MsB;

    /* Stripe magnetization: Mst1 = (pi,0), Mst2 = (0,pi) */
    double **Mst1, **Mst2;
    double **Sum_Mst1, **Sum_Mst1_2, **Sum_Mst1_4;
    double **Sum_Mst2, **Sum_Mst2_2, **Sum_Mst2_4;
    double **Ave_Mst1, **Ave_Mst1_2, **Ave_Mst1B;
    double **Ave_Mst2, **Ave_Mst2_2, **Ave_Mst2B;
    double **Err_Mst1, **Err_Mst1_2, **Err_Mst1B;
    double **Err_Mst2, **Err_Mst2_2, **Err_Mst2B;
    double *All_Mst1, *All_Mst1_2, *All_Mst1B;
    double *All_Mst2, *All_Mst2_2, *All_Mst2B;
    double *All_Err_Mst1, *All_Err_Mst1_2, *All_Err_Mst1B;
    double *All_Err_Mst2, *All_Err_Mst2_2, *All_Err_Mst2B;

    double ***Sum_Sqx, ***Sum_Sqy, ***Sum_Sqz;
    double *All_Sqx, *All_Sqy, *All_Sqz;
    double *AG_Sqx, *AG_Sqy, *AG_Sqz;

    /* Z2 chirality (bond nematic order parameter) — scalar */
    double *KappaZ2; /* [num_temp] instantaneous raw sum K */
    double *Sum_KappaZ2;
    double *Sum_KappaZ2_2, *Sum_KappaZ2_4; /* [num_temp] sum(K^2), sum(K^4) */

    double *Ave_KappaZ2, *Err_KappaZ2;     /* <kappa> */
    double *Ave_KappaZ2_2, *Err_KappaZ2_2; /* <kappa^2> */
    double *Ave_KappaZ2B, *Err_KappaZ2B;   /* Binder */

    double *All_KappaZ2, *All_Err_KappaZ2;
    double *All_KappaZ2_2, *All_Err_KappaZ2_2;
    double *All_KappaZ2B, *All_Err_KappaZ2B;

    /* NER kappa_Z2 */
    double NER_kappa_Z2; /* normalized kappa = K/N_p */

    /* NER overlap with initial configuration */
    double NER_overlap[3];    /* S(t)·S(0) per component (x,y,z) */
    double NER_overlap_total; /* sum of all components */
};

struct TimeList {
    time_t start, mid1, mid2, end;
};

struct BindStruct {
    struct DefineList Def;
    struct PhysList Phys;
    struct TimeList Time;
};

struct MCMainCalStruct {
    struct BindStruct Bind;
};

/*=== Memory macros (from mfmemory08.c) ===*/
#define d_malloc1(X, N1) X = (double *)malloc((N1) * sizeof(double));

#define d_malloc2(X, N1, N2)                                                   \
    X = (double **)malloc((N1) * sizeof(double *));                            \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        X[mfint[0]] = (double *)malloc((N2) * sizeof(double));                 \
    }

#define d_malloc3(X, N1, N2, N3)                                               \
    X = (double ***)malloc((N1) * sizeof(double **));                          \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        X[mfint[0]] = (double **)malloc((N2) * sizeof(double *));              \
        for (mfint[1] = 0; mfint[1] < (N2); mfint[1]++) {                      \
            X[mfint[0]][mfint[1]] = (double *)malloc((N3) * sizeof(double));   \
        }                                                                      \
    }

#define d_free1(X, N1) free(X);

#define d_free2(X, N1, N2)                                                     \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        free(X[mfint[0]]);                                                     \
    }                                                                          \
    free(X);

#define d_free3(X, N1, N2, N3)                                                 \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        for (mfint[1] = 0; mfint[1] < (N2); mfint[1]++) {                      \
            free(X[mfint[0]][mfint[1]]);                                       \
        }                                                                      \
        free(X[mfint[0]]);                                                     \
    }                                                                          \
    free(X);

#define i_malloc1(X, N1) X = (int *)malloc((N1) * sizeof(int));

#define i_malloc2(X, N1, N2)                                                   \
    X = (int **)malloc((N1) * sizeof(int *));                                  \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        X[mfint[0]] = (int *)malloc((N2) * sizeof(int));                       \
    }

#define i_free1(X, N1) free(X);

#define i_free2(X, N1, N2)                                                     \
    for (mfint[0] = 0; mfint[0] < (N1); mfint[0]++) {                          \
        free(X[mfint[0]]);                                                     \
    }                                                                          \
    free(X);

/*=== Function prototypes (stage1 minimal set) ===*/
/* input_parser.c */
int read_param(const char *filename, struct DefineList *Def);
int read_lattice(const char *filename, struct DefineList *Def);
int read_interaction(const char *filename, struct DefineList *Def);

/* lattice.c */
int build_lattice(const char *filename, struct BindStruct *X);
void initial(dsfmt_t *dsfmt, struct BindStruct *X);

/* mc_update.c */
void MC(dsfmt_t *dsfmt, struct BindStruct *X);
/* EXMC: try exchange between temperature slots i and j; return 1 if accepted */
int exchange_try_pair(dsfmt_t *dsfmt, struct BindStruct *X, int i, int j);
/* EXMC: one sweep of adjacent-pair exchanges (odd-even); return number accepted */
int exchange_sweep_local(dsfmt_t *dsfmt, struct BindStruct *X);

#endif /* MC_DEF_H */
