#include "memory.h"

/* Allocate 2D double array: a[n1][n2], zero-initialized. */
static double **alloc_double_2d(int n1, int n2) {
    int i;
    double **a;

    a = (double **)malloc((size_t)n1 * sizeof(double *));
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

/* Allocate 2D int array: a[n1][n2], initialized with init_val. */
static int **alloc_int_2d(int n1, int n2, int init_val) {
    int i, j;
    int **a;

    a = (int **)malloc((size_t)n1 * sizeof(int *));
    if (a == NULL)
        return NULL;

    for (i = 0; i < n1; i++) {
        a[i] = (int *)malloc((size_t)n2 * sizeof(int));
        if (a[i] == NULL) {
            int k;
            for (k = 0; k < i; k++)
                free(a[k]);
            free(a);
            return NULL;
        }
        for (j = 0; j < n2; j++)
            a[i][j] = init_val;
    }
    return a;
}

/* Free helper for alloc_double_2d(). */
static void free_double_2d(double **a, int n1) {
    int i;
    if (a == NULL)
        return;
    for (i = 0; i < n1; i++)
        free(a[i]);
    free(a);
}

/* Free helper for alloc_int_2d(). */
static void free_int_2d(int **a, int n1) {
    int i;
    if (a == NULL)
        return;
    for (i = 0; i < n1; i++)
        free(a[i]);
    free(a);
}

/*
 * Allocate only the arrays required by stage1 runtime:
 * spins, effective fields, neighbor table, energy, acceptance counters.
 */
int allocate_minimal_mc_arrays(struct MCMainCalStruct *X) {
    int num_temp, All_N, ni_max;

    num_temp = X->Bind.Def.num_temp;
    All_N = X->Bind.Def.All_N;
    ni_max = X->Bind.Def.ni_max;

    if (num_temp <= 0 || All_N <= 0 || ni_max <= 0)
        return -1;

    X->Bind.Def.sx = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.sy = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.sz = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.env_sx = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.env_sy = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.env_sz = alloc_double_2d(num_temp, All_N);
    X->Bind.Def.Nr = alloc_int_2d(ni_max, All_N, -1);
    X->Bind.Def.J = alloc_double_2d(ni_max, All_N);
    X->Bind.Def.stag_sign = (int *)calloc((size_t)All_N, sizeof(int));
    X->Bind.Def.stripe_sign_1 = (int *)calloc((size_t)All_N, sizeof(int));
    X->Bind.Def.stripe_sign_2 = (int *)calloc((size_t)All_N, sizeof(int));
    X->Bind.Phys.Energy = (double *)calloc((size_t)num_temp, sizeof(double));
    X->Bind.Phys.ratio_1 = (int *)calloc((size_t)num_temp, sizeof(int));

    if (X->Bind.Def.sx == NULL || X->Bind.Def.sy == NULL ||
        X->Bind.Def.sz == NULL || X->Bind.Def.env_sx == NULL ||
        X->Bind.Def.env_sy == NULL || X->Bind.Def.env_sz == NULL ||
        X->Bind.Def.Nr == NULL || X->Bind.Def.J == NULL ||
        X->Bind.Def.stag_sign == NULL || X->Bind.Def.stripe_sign_1 == NULL ||
        X->Bind.Def.stripe_sign_2 == NULL || X->Bind.Phys.Energy == NULL ||
        X->Bind.Phys.ratio_1 == NULL) {
        free_minimal_mc_arrays(X);
        return -1;
    }

    return 0;
}

/* Free all arrays allocated by allocate_minimal_mc_arrays(). */
void free_minimal_mc_arrays(struct MCMainCalStruct *X) {
    int num_temp, ni_max;

    num_temp = X->Bind.Def.num_temp;
    ni_max = X->Bind.Def.ni_max;

    free_double_2d(X->Bind.Def.sx, num_temp);
    free_double_2d(X->Bind.Def.sy, num_temp);
    free_double_2d(X->Bind.Def.sz, num_temp);
    free_double_2d(X->Bind.Def.env_sx, num_temp);
    free_double_2d(X->Bind.Def.env_sy, num_temp);
    free_double_2d(X->Bind.Def.env_sz, num_temp);
    free_int_2d(X->Bind.Def.Nr, ni_max);
    free_double_2d(X->Bind.Def.J, ni_max);

    free(X->Bind.Def.stag_sign);
    free(X->Bind.Def.stripe_sign_1);
    free(X->Bind.Def.stripe_sign_2);

    free(X->Bind.Phys.Energy);
    free(X->Bind.Phys.ratio_1);
}

/*
 * Allocate temporary/accumulation arrays for outputs:
 * E, E^2, specific heat, M^2, and acceptance during step/sample averaging.
 */
int allocate_work_arrays(int num_temp, struct SimpleWorkArrays *W) {
    int num_pairs;

    if (num_temp <= 0 || W == NULL)
        return -1;

    num_pairs = (num_temp > 1) ? (num_temp - 1) : 0;

    W->accum_E = (double *)calloc((size_t)num_temp, sizeof(double));
    W->accum_C = (double *)calloc((size_t)num_temp, sizeof(double));
    W->accum_M2 = (double *)calloc((size_t)num_temp, sizeof(double));
    W->accum_A = (double *)calloc((size_t)num_temp, sizeof(double));
    W->sample_E = (double *)calloc((size_t)num_temp, sizeof(double));
    W->sample_E2 = (double *)calloc((size_t)num_temp, sizeof(double));
    W->sample_M2 = (double *)calloc((size_t)num_temp, sizeof(double));
    W->exchange_accept = (num_pairs > 0)
                             ? (int *)calloc((size_t)num_pairs, sizeof(int))
                             : NULL;
    W->exchange_attempt = (num_pairs > 0)
                              ? (int *)calloc((size_t)num_pairs, sizeof(int))
                              : NULL;
    W->accum_exchange_accept = (num_pairs > 0)
                                   ? (int *)calloc((size_t)num_pairs, sizeof(int))
                                   : NULL;
    W->accum_exchange_attempt = (num_pairs > 0)
                                    ? (int *)calloc((size_t)num_pairs, sizeof(int))
                                    : NULL;

    if (W->accum_E == NULL || W->accum_C == NULL || W->accum_M2 == NULL ||
        W->accum_A == NULL || W->sample_E == NULL || W->sample_E2 == NULL ||
        W->sample_M2 == NULL) {
        free_work_arrays(W);
        return -1;
    }
    if (num_pairs > 0 &&
        (W->exchange_accept == NULL || W->exchange_attempt == NULL ||
         W->accum_exchange_accept == NULL || W->accum_exchange_attempt == NULL)) {
        free_work_arrays(W);
        return -1;
    }

    return 0;
}

/* Free arrays allocated by allocate_work_arrays(). */
void free_work_arrays(struct SimpleWorkArrays *W) {
    if (W == NULL)
        return;
    free(W->accum_E);
    free(W->accum_C);
    free(W->accum_M2);
    free(W->accum_A);
    free(W->sample_E);
    free(W->sample_E2);
    free(W->sample_M2);
    free(W->exchange_accept);
    free(W->exchange_attempt);
    free(W->accum_exchange_accept);
    free(W->accum_exchange_attempt);
}
