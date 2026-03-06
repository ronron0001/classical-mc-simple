#ifndef MEMORY_H
#define MEMORY_H

#include "mc_def.h"

struct SimpleWorkArrays {
    double *accum_E;   /* Sum over samples of per-sample <E/N> */
    double *accum_C;   /* Sum over samples of per-sample specific heat */
    double *accum_M2;  /* Sum over samples of per-sample <M^2> */
    double *accum_A;   /* Sum over samples of per-sample acceptance */
    double *sample_E;  /* Sum over MC steps within one sample, using E/N */
    double *sample_E2; /* Sum over MC steps within one sample, using E^2 */
    double *sample_M2; /* Sum over MC steps within one sample */
    /* Exchange MC: per-pair (i,i+1) counts, length num_temp-1 */
    int *exchange_accept; /* current sample */
    int *exchange_attempt;
    int *accum_exchange_accept; /* sum over samples */
    int *accum_exchange_attempt;
};

int allocate_minimal_mc_arrays(struct MCMainCalStruct *X);
void free_minimal_mc_arrays(struct MCMainCalStruct *X);

int allocate_work_arrays(int num_temp, struct SimpleWorkArrays *W);
void free_work_arrays(struct SimpleWorkArrays *W);

#endif
