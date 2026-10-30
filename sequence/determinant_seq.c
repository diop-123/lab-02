/*
 * Variant 9 — Calcul du déterminant d'une matrice par définition
 * (formule de Leibniz)
 *
 *   det(A) = Σ_{σ ∈ S_n} sign(σ) · Π_{i=0}^{n-1} a[i][σ(i)]
 *
 * Version SÉQUENTIELLE uniquement.
 *
 * Compilation :
 *   gcc -O2 -Wall -Wextra -o determinant_seq determinant_seq.c -lm
 *
 * Utilisation :
 *   ./determinant_seq --n 8
 *   ./determinant_seq --input matrix.txt
 *   ./determinant_seq --n 7 --verbose
 *   ./determinant_seq --n 9 --repeat 5
 *
 * Format du fichier --input :
 *   ligne 1 : n
 *   puis n*n valeurs (ligne par ligne)
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 13   /* au-delà de 12-13, n! devient trop grand */

/* ------------------------------------------------------------------ */
/* Table des factorielles (exacte jusqu'à 20!)                        */
/* ------------------------------------------------------------------ */
static unsigned long long FACT[21];

static void init_factorials(void) {
    FACT[0] = 1ULL;
    for (int i = 1; i <= 20; i++)
        FACT[i] = FACT[i - 1] * (unsigned long long)i;
}

/* ------------------------------------------------------------------ */
/* Dé-codage du rang -> permutation (système factoradique)            */
/*                                                                    */
/* Le rang r ∈ [0, n!) est écrit en base factorielle :               */
/*   r = c_0·(n-1)! + c_1·(n-2)! + ... + c_{n-1}·0!                  */
/* On pioche ensuite dans la liste des éléments disponibles.          */
/* ------------------------------------------------------------------ */
static void unrank_perm(unsigned long long rank, int n, int *perm) {
    int avail[MAX_N];
    for (int i = 0; i < n; i++) avail[i] = i;

    for (int i = 0; i < n; i++) {
        int rem = n - i;                 /* éléments encore disponibles */
        unsigned long long f = FACT[rem - 1];
        int idx = (int)(rank / f);
        rank %= f;

        perm[i] = avail[idx];

        /* retire avail[idx] en décalant */
        for (int j = idx; j < rem - 1; j++)
            avail[j] = avail[j + 1];
    }
}

/* ------------------------------------------------------------------ */
/* Signature d'une permutation = (-1)^(nombre d'inversions)          */
/* ------------------------------------------------------------------ */
static int perm_sign(const int *perm, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (perm[i] > perm[j]) inv++;
    return (inv & 1) ? -1 : 1;
}

/* ------------------------------------------------------------------ */
/* Produit  a[0][perm[0]] · a[1][perm[1]] · ...                      */
/* La matrice est stockée en 1D (row-major) : A[i*n + j]             */
/* ------------------------------------------------------------------ */
static long double perm_product(const long double *A, int n, const int *perm) {
    long double p = 1.0L;
    for (int i = 0; i < n; i++)
        p *= A[i * n + perm[i]];
    return p;
}

/* ------------------------------------------------------------------ */
/* Déterminant séquentiel (formule de Leibniz)                        */
/* ------------------------------------------------------------------ */
static long double det_seq(const long double *A, int n) {
    unsigned long long total = FACT[n];
    int perm[MAX_N];
    long double sum = 0.0L;

    for (unsigned long long r = 0; r < total; r++) {
        unrank_perm(r, n, perm);
        int s = perm_sign(perm, n);
        long double p = perm_product(A, n, perm);
        sum += (long double)s * p;
    }
    return sum;
}

/* ------------------------------------------------------------------ */
/* Chronomètre monotone en millisecondes                              */
/* ------------------------------------------------------------------ */
static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv) {
    init_factorials();

    int n = 0;
    int verbose = 0;
    int repeat = 1;
    const char *input = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--n")      && i + 1 < argc) n       = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--input")  && i + 1 < argc) input   = argv[++i];
        else if (!strcmp(argv[i], "--repeat") && i + 1 < argc) repeat  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--verbose"))                verbose = 1;
        else {
            fprintf(stderr, "Argument inconnu : %s\n", argv[i]);
            return 1;
        }
    }

    if (!input && n <= 0) {
        fprintf(stderr,
            "Usage : %s --n N [--repeat R] [--verbose]\n"
            "     ou %s --input fichier [--repeat R] [--verbose]\n",
            argv[0], argv[0]);
        return 1;
    }

    long double *A = NULL;

    /* ---- Lecture depuis un fichier ------------------------------- */
    if (input) {
        FILE *f = fopen(input, "r");
        if (!f) { perror("fopen"); return 1; }

        if (fscanf(f, "%d", &n) != 1 || n < 1 || n > MAX_N) {
            fprintf(stderr, "Entête de matrice invalide (n attendu 1..%d)\n", MAX_N);
            fclose(f);
            return 1;
        }

        A = malloc(sizeof(long double) * (size_t)n * (size_t)n);
        if (!A) { perror("malloc"); fclose(f); return 1; }

        for (int i = 0; i < n * n; i++) {
            double v;
            if (fscanf(f, "%lf", &v) != 1) {
                fprintf(stderr, "Valeur manquante à l'indice %d\n", i);
                free(A); fclose(f);
                return 1;
            }
            A[i] = (long double)v;
        }
        fclose(f);
    }
    /* ---- Ou génération pseudo-aléatoire ------------------------- */
    else {
        if (n < 1 || n > MAX_N) {
            fprintf(stderr, "n hors bornes (1..%d)\n", MAX_N);
            return 1;
        }
        A = malloc(sizeof(long double) * (size_t)n * (size_t)n);
        if (!A) { perror("malloc"); return 1; }

        srand(42);   /* graine fixe : mêmes entrées à chaque exécution */
        for (int i = 0; i < n * n; i++)
            A[i] = (long double)(rand() % 19 - 9);   /* valeurs entières dans [-9,9] */
    }

    /* ---- Affichage optionnel de la matrice ----------------------- */
    if (verbose) {
        printf("Matrice %dx%d :\n", n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                printf("%8.2Lf", A[i * n + j]);
            printf("\n");
        }
    }

    if (repeat < 1) repeat = 1;

    /* ---- Mesure ------------------------------------------------ */
    long double det = 0.0L;
    double best_ms = 1e30;

    for (int k = 0; k < repeat; k++) {
        double t0 = now_ms();
        det = det_seq(A, n);
        double t1 = now_ms();
        double dt = t1 - t0;
        if (dt < best_ms) best_ms = dt;
        if (verbose && repeat > 1)
            printf("  run %d : %.3f ms\n", k + 1, dt);
    }

    /* ---- Résultat ---------------------------------------------- */
    printf("mode=seq n=%d repeat=%d time_ms=%.3f det=%.6Lf\n",
           n, repeat, best_ms, det);

    if (n >= 10) {
        fprintf(stderr,
            "Rappel : n=%d => %llu permutations (%.3e), temps en O(n·n!).\n",
            n, FACT[n], (double)FACT[n]);
    }

    free(A);
    return 0;
}