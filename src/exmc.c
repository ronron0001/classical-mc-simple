#include "mc_def.h"

static void swap_state(struct BindStruct *X, int i, int j) {
    double *tmp_d;
    double tmp;

    tmp_d = X->Def.sx[i];
    X->Def.sx[i] = X->Def.sx[j];
    X->Def.sx[j] = tmp_d;

    tmp_d = X->Def.sy[i];
    X->Def.sy[i] = X->Def.sy[j];
    X->Def.sy[j] = tmp_d;

    tmp_d = X->Def.sz[i];
    X->Def.sz[i] = X->Def.sz[j];
    X->Def.sz[j] = tmp_d;

    tmp_d = X->Def.env_sx[i];
    X->Def.env_sx[i] = X->Def.env_sx[j];
    X->Def.env_sx[j] = tmp_d;

    tmp_d = X->Def.env_sy[i];
    X->Def.env_sy[i] = X->Def.env_sy[j];
    X->Def.env_sy[j] = tmp_d;

    tmp_d = X->Def.env_sz[i];
    X->Def.env_sz[i] = X->Def.env_sz[j];
    X->Def.env_sz[j] = tmp_d;

    tmp = X->Phys.Energy[i];
    X->Phys.Energy[i] = X->Phys.Energy[j];
    X->Phys.Energy[j] = tmp;
}

void exchange_step(dsfmt_t *dsfmt, struct BindStruct *X, int *exchange_accept,
                   int *exchange_attempt) {
    int num_temp = X->Def.num_temp;
    double Ini_T = X->Def.Ini_T;
    double Delta_T = X->Def.Delta_T;
    int parity, i;

    if (num_temp < 2) {
        return;
    }

    parity = (int)(dsfmt_genrand_close_open(dsfmt) * 2.0);
    for (i = parity; i < num_temp - 1; i += 2) {
        double Ti = Ini_T + Delta_T * (double)i;
        double Tj = Ini_T + Delta_T * (double)(i + 1);
        double beta_i = 1.0 / Ti;
        double beta_j = 1.0 / Tj;
        double Ei = X->Phys.Energy[i];
        double Ej = X->Phys.Energy[i + 1];
        double logw = (beta_i - beta_j) * (Ei - Ej);
        double r = dsfmt_genrand_close_open(dsfmt);
        int accept = (logw >= 0.0) || (r < exp(logw));

        if (exchange_attempt != NULL) {
            exchange_attempt[i] += 1;
        }
        if (accept) {
            swap_state(X, i, i + 1);
            if (exchange_accept != NULL) {
                exchange_accept[i] += 1;
            }
        }
    }
}
