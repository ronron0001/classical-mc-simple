/*-------------------------------------------------------------
 * input_parser.c
 *
 * Functions for reading simulation input files:
 *   read_param       -- reads param.def       (simulation parameters)
 *   read_lattice     -- reads lattice.def      (lattice geometry)
 *   read_interaction -- reads interaction.def  (bond counting only)
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Takahiro Misawa
 *-------------------------------------------------------------*/

#include "mc_def.h"

#define LINE_BUF 1024
#define MAX_ORB 100

/*-------------------------------------------------------------
 * is_comment_or_blank
 *   Returns 1 if the line is a comment (starts with '#')
 *   or is blank (only whitespace / newline). Returns 0 otherwise.
 *-------------------------------------------------------------*/
static int is_comment_or_blank(const char *line) {
    int i = 0;
    /* skip leading whitespace */
    while (line[i] == ' ' || line[i] == '\t') {
        i++;
    }
    if (line[i] == '#' || line[i] == '\n' || line[i] == '\r' ||
        line[i] == '\0') {
        return 1;
    }
    return 0;
}

/*-------------------------------------------------------------
 * read_param
 *
 * Reads a parameter file with "key = value" format.
 * Sets: Def->Burn_in, Def->Total_Step, Def->Sample,
 *       Def->num_temp, Def->Ini_T, Def->Delta_T,
 *       Def->lambda, Def->H, Def->exchange_interval
 *
 * Returns 0 on success, -1 on error.
 *-------------------------------------------------------------*/
int read_param(const char *filename, struct DefineList *Def) {
    FILE *fp;
    char line[LINE_BUF];
    char key[LINE_BUF];
    char val[LINE_BUF];

    /* Default values */
    Def->spin_dim = 3;    /* Heisenberg */
    Def->output_spin = 0; /* skip spin output by default */
    Def->init_state = 0;  /* RANDOM */
    Def->exchange_interval = 1; /* enable exchange every step by default */

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file '%s'\n", filename);
        return -1;
    }

    while (fgets(line, LINE_BUF, fp) != NULL) {
        if (is_comment_or_blank(line)) {
            continue;
        }
        if (sscanf(line, "%s = %s", key, val) != 2) {
            continue;
        }

        if (strcmp(key, "Burn_in") == 0) {
            Def->Burn_in = atoi(val);
        } else if (strcmp(key, "Total_Step") == 0) {
            Def->Total_Step = atoi(val);
        } else if (strcmp(key, "Sample") == 0) {
            Def->Sample = atoi(val);
        } else if (strcmp(key, "num_temp") == 0) {
            Def->num_temp = atoi(val);
        } else if (strcmp(key, "Ini_T") == 0) {
            Def->Ini_T = atof(val);
        } else if (strcmp(key, "Delta_T") == 0) {
            Def->Delta_T = atof(val);
        } else if (strcmp(key, "lambda") == 0) {
            Def->lambda = atof(val);
        } else if (strcmp(key, "H") == 0) {
            Def->H = atof(val);
        } else if (strcmp(key, "exchange_interval") == 0) {
            Def->exchange_interval = atoi(val);
        } else if (strcmp(key, "spin_dim") == 0) {
            Def->spin_dim = atoi(val);
        } else if (strcmp(key, "output_spin") == 0) {
            Def->output_spin = atoi(val);
        } else if (strcmp(key, "init_state") == 0) {
            Def->init_state = atoi(val);
        } else {
            fprintf(stderr, "Warning: unknown key '%s' in '%s'\n", key,
                    filename);
        }
    }

    fclose(fp);

    if (Def->myrank == MASTER) {
        printf("read_param: Burn_in=%d Total_Step=%d Sample=%d num_temp=%d\n",
               Def->Burn_in, Def->Total_Step, Def->Sample, Def->num_temp);
        printf("read_param: exchange_interval=%d\n", Def->exchange_interval);
        printf(
            "read_param: Ini_T=%lf Delta_T=%lf lambda=%lf H=%lf spin_dim=%d\n",
            Def->Ini_T, Def->Delta_T, Def->lambda, Def->H, Def->spin_dim);
    }

    return 0;
}

/*-------------------------------------------------------------
 * read_lattice
 *
 * Reads a lattice definition file with "key = value" format.
 * Sets: Def->L_x, Def->L_y, Def->L_z, Def->orb_num
 * Computes: Def->All_N = L_x * L_y * L_z * orb_num
 *
 * Returns 0 on success, -1 on error.
 *-------------------------------------------------------------*/
int read_lattice(const char *filename, struct DefineList *Def) {
    FILE *fp;
    char line[LINE_BUF];
    char key[LINE_BUF];
    char val[LINE_BUF];

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file '%s'\n", filename);
        return -1;
    }

    /* Initialize to safe defaults */
    Def->L_x = 0;
    Def->L_y = 0;
    Def->L_z = 1;
    Def->orb_num = 1;

    while (fgets(line, LINE_BUF, fp) != NULL) {
        if (is_comment_or_blank(line)) {
            continue;
        }
        if (sscanf(line, "%s = %s", key, val) != 2) {
            continue;
        }

        if (strcmp(key, "L_x") == 0) {
            Def->L_x = atoi(val);
        } else if (strcmp(key, "L_y") == 0) {
            Def->L_y = atoi(val);
        } else if (strcmp(key, "L_z") == 0) {
            Def->L_z = atoi(val);
        } else if (strcmp(key, "orb_num") == 0) {
            Def->orb_num = atoi(val);
        } else {
            fprintf(stderr, "Warning: unknown key '%s' in '%s'\n", key,
                    filename);
        }
    }

    fclose(fp);

    Def->All_N = Def->L_x * Def->L_y * Def->L_z * Def->orb_num;

    if (Def->myrank == MASTER) {
        printf("read_lattice: L_x=%d L_y=%d L_z=%d orb_num=%d All_N=%d\n",
               Def->L_x, Def->L_y, Def->L_z, Def->orb_num, Def->All_N);
    }

    return 0;
}

/*-------------------------------------------------------------
 * read_interaction
 *
 * Reads an interaction definition file (bond list).
 * Each data line has: i_orb  j_orb  dx  dy  dz  J
 *
 * This function performs a COUNTING PASS ONLY:
 *   - Counts the number of neighbors for each orbital
 *     (each bond contributes +1 to i_orb and +1 to j_orb)
 *   - Determines ni_max = max neighbor count across all orbitals
 *   - Sets Def->ni_max
 *
 * The actual Nr/J arrays are built later in lattice.c's initial()
 * after xsetmem has allocated them using ni_max.
 *
 * Returns the number of unique bonds on success, -1 on error.
 *-------------------------------------------------------------*/
int read_interaction(const char *filename, struct DefineList *Def) {
    FILE *fp;
    char line[LINE_BUF];
    int i_orb, j_orb, dx, dy, dz;
    double J_val;
    int count[MAX_ORB];
    int orb_num;
    int num_bonds;
    int ni_max;
    int i;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open file '%s'\n", filename);
        return -1;
    }

    orb_num = Def->orb_num;
    if (orb_num <= 0 || orb_num > MAX_ORB) {
        fprintf(stderr, "Error: orb_num=%d out of range (1..%d)\n", orb_num,
                MAX_ORB);
        fclose(fp);
        return -1;
    }

    /* Initialize neighbor count per orbital */
    for (i = 0; i < orb_num; i++) {
        count[i] = 0;
    }

    num_bonds = 0;

    while (fgets(line, LINE_BUF, fp) != NULL) {
        if (is_comment_or_blank(line)) {
            continue;
        }
        if (sscanf(line, "%d %d %d %d %d %lf", &i_orb, &j_orb, &dx, &dy, &dz,
                   &J_val) != 6) {
            fprintf(stderr, "Warning: skipping malformed line in '%s': %s",
                    filename, line);
            continue;
        }

        if (i_orb < 0 || i_orb >= orb_num || j_orb < 0 || j_orb >= orb_num) {
            fprintf(stderr, "Warning: orbital index out of range in '%s': %s",
                    filename, line);
            continue;
        }

        /* Forward direction: i_orb gets one more neighbor */
        count[i_orb]++;
        /* Reverse direction: j_orb gets one more neighbor */
        count[j_orb]++;

        num_bonds++;
    }

    fclose(fp);

    /* Determine ni_max = maximum neighbor count across all orbitals */
    ni_max = 0;
    for (i = 0; i < orb_num; i++) {
        if (count[i] > ni_max) {
            ni_max = count[i];
        }
    }

    /* Verify all orbitals have the same neighbor count; warn if not */
    for (i = 0; i < orb_num; i++) {
        if (count[i] != ni_max) {
            if (Def->myrank == MASTER) {
                fprintf(stderr,
                        "Warning: orb %d has %d neighbors, but ni_max=%d\n", i,
                        count[i], ni_max);
            }
        }
    }

    Def->ni_max = ni_max;

    if (Def->myrank == MASTER) {
        printf("read_interaction: Read %d unique bonds, ni_max = %d\n",
               num_bonds, ni_max);
    }

    return num_bonds;
}
