/*
 * Variant 9 — Déterminant d'une matrice par définition (Leibniz)
 *
 *   det(A) = Σ_{σ ∈ S_n} sign(σ) · Π_{i=0}^{n-1} a[i][σ(i)]
 *
 * Version PARALLÈLE avec POSIX Threads (pthreads), Linux / Ubuntu.
 *
 * Contraintes respectées :
 *   - uniquement pthreads (pas d'OpenMP, pas de std::thread, pas de C11 threads)
 *   - nombre maximal de threads passé en argument (--threads)
 *   - mesure du temps, répétable (--repeat), vérification (--verify)
 *   - option --hold pour inspecter les threads avec ps -L / top -H
 *
 * Répartition : découpage STATIQUE en tranches égales de l'intervalle
 *   [0, n!)  (chaque permutation a le même coût, donc le découpage
 *    statique est parfaitement équilibré : pas besoin de mutex pour un
 *    compteur dynamique, ce qui évite tout surcoût de synchronisation).
 *
 * Compilation :
 *   gcc -O2 -Wall -Wextra -pthread -o determinant_par determinant_par.c -lm
 *
 * Utilisation :
 *   ./determinant_par --n 9 --threads 8
 *   ./determinant_par --n 10 --threads 4 --repeat 5
 *   ./determinant_par --input matrix.txt --threads 8 --verify
 *   ./determinant_par --n 11 --threads 8 --hold 20     # inspection OS
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#define MAX_N 13   /* n! devient trop grand au-delà */

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
/* Dé-codage rang -> permutation (système factoradique)               */
/* ------------------------------------------------------------------ */
static void unrank_perm(unsigned long long rank, int n, int *perm) {
    int avail[MAX_N];
    for (int i = 0; i < n; i++) avail[i] = i;

    for (int i = 0; i < n; i++) {
        int rem = n - i;
        unsigned long long f = FACT[rem - 1];
        int idx = (int)(rank / f);
        rank %= f;

        perm[i] = avail[idx];
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
/* Matrice 1D row-major : A[i*n + j]                                  */
/* ------------------------------------------------------------------ */
static long double perm_product(const long double *A, int n, const int *perm) {
    long double p = 1.0L;
    for (int i = 0; i < n; i++)
        p *= A[i * n + perm[i]];
    return p;
}

/* ------------------------------------------------------------------ */
/* Référence séquentielle (uniquement pour --verify)                  */
/* ------------------------------------------------------------------ */
static long double det_seq(const long double *A, int n) {
    unsigned long long total = FACT[n];
    int perm[MAX_N];
    long double sum = 0.0L;
    for (unsigned long long r = 0; r < total; r++) {
        unrank_perm(r, n, perm);
        sum += (long double)perm_sign(perm, n) * perm_product(A, n, perm);
    }
    return sum;
}

/* ------------------------------------------------------------------ */
/* Structure passée à chaque thread                                   */
/* ------------------------------------------------------------------ */
typedef struct {
    int tid;
    int n;
    unsigned long long start;   /* inclus */
    unsigned long long end;     /* exclu  */
    const long double *A;
    long double partial;
} thread_arg_t;

/* ------------------------------------------------------------------ */
/* Fonction exécutée par chaque thread                                */
/* ------------------------------------------------------------------ */
static void *worker(void *arg) {
    thread_arg_t *ta = (thread_arg_t *)arg;
    int n = ta->n;
    int perm[MAX_N];
    long double sum = 0.0L;

    for (unsigned long long r = ta->start; r < ta->end; r++) {
        unrank_perm(r, n, perm);
        sum += (long double)perm_sign(perm, n) * perm_product(ta->A, n, perm);
    }
    ta->partial = sum;
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Déterminant parallèle                                              */
/* ------------------------------------------------------------------ */
static long double det_par(const long double *A, int n, int nthreads, int verbose) {
    unsigned long long total = FACT[n];

    /* on ne crée jamais plus de threads que de permutations à traiter */
    if ((unsigned long long)nthreads > total)
        nthreads = (int)total;
    if (nthreads < 1) nthreads = 1;

    thread_arg_t *args = calloc((size_t)nthreads, sizeof(*args));
    pthread_t    *tids = calloc((size_t)nthreads, sizeof(*tids));
    if (!args || !tids) { perror("calloc"); exit(EXIT_FAILURE); }

    /* Découpage statique : chunks de taille total/nthreads, +1 pour les
       (total % nthreads) premiers threads.  Chaque permutation coûtant
       exactement le même temps, ce découpage est parfaitement équilibré. */
    unsigned long long chunk = total / (unsigned long long)nthreads;
    unsigned long long extra = total % (unsigned long long)nthreads;
    unsigned long long pos   = 0;

    for (int t = 0; t < nthreads; t++) {
        unsigned long long sz = chunk + ((unsigned long long)t < extra ? 1ULL : 0ULL);

        args[t].tid     = t;
        args[t].n       = n;
        args[t].start   = pos;
        args[t].end     = pos + sz;
        args[t].A       = A;
        args[t].partial = 0.0L;
        pos += sz;

        int rc = pthread_create(&tids[t], NULL, worker, &args[t]);
        if (rc != 0) {
            fprintf(stderr, "pthread_create(t=%d) : %s\n", t, strerror(rc));
            exit(EXIT_FAILURE);
        }
    }

    long double sum = 0.0L;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(tids[t], NULL);
        sum += args[t].partial;
        if (verbose)
            fprintf(stderr,
                "  thread %2d : perms [%llu .. %llu)  partial=%.6Lf\n",
                t, args[t].start, args[t].end, args[t].partial);
    }

    free(args);
    free(tids);
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

    int n        = 0;
    int nthreads = 1;
    int verbose  = 0;
    int repeat   = 1;
    int hold_s   = 0;
    int verify   = 0;
    const char *input = NULL;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--n")       && i + 1 < argc) n        = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--input")   && i + 1 < argc) input    = argv[++i];
        else if (!strcmp(argv[i], "--repeat")  && i + 1 < argc) repeat   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--hold")    && i + 1 < argc) hold_s   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--verify"))                  verify   = 1;
        else if (!strcmp(argv[i], "--verbose"))                 verbose  = 1;
        else {
            fprintf(stderr, "Argument inconnu : %s\n", argv[i]);
            return 1;
        }
    }

    if (!input && n <= 0) {
        fprintf(stderr,
            "Usage : %s --n N --threads T [--repeat R] [--verify] [--verbose] [--hold S]\n"
            "     ou %s --input fichier --threads T [--repeat R] [--verify] [--verbose] [--hold S]\n",
            argv[0], argv[0]);
        return 1;
    }

    /* --------- Chargement de la matrice ------------------------- */
    long double *A = NULL;

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
    } else {
        if (n < 1 || n > MAX_N) { fprintf(stderr, "n hors bornes (1..%d)\n", MAX_N); return 1; }
        A = malloc(sizeof(long double) * (size_t)n * (size_t)n);
        if (!A) { perror("malloc"); return 1; }

        srand(42);  /* graine fixe : matrice identique à chaque exécution */
        for (int i = 0; i < n * n; i++)
            A[i] = (long double)(rand() % 19 - 9);
    }

    if (verbose) {
        printf("Matrice %dx%d :\n", n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                printf("%8.2Lf", A[i * n + j]);
            printf("\n");
        }
    }

    /* --------- Attente optionnelle pour inspecter les threads --- */
    if (hold_s > 0) {
        printf("PID=%d — pause %d s. Inspection depuis un autre terminal :\n"
               "  ps -L -p %d | tail -n +2 | wc -l     # nb threads\n"
               "  top -H -p %d                        # vue temps réel\n",
               getpid(), hold_s, getpid(), getpid());
        fflush(stdout);
        sleep(hold_s);
    }

    if (repeat < 1) repeat = 1;

    /* --------- Mesures ----------------------------------------- */
    long double det = 0.0L;
    double best_ms = 1e30;

    for (int k = 0; k < repeat; k++) {
        double t0 = now_ms();
        det = det_par(A, n, nthreads, verbose && (k == 0));
        double t1 = now_ms();
        double dt = t1 - t0;
        if (dt < best_ms) best_ms = dt;
        if (verbose && repeat > 1)
            fprintf(stderr, "  run %d : %.3f ms\n", k + 1, dt);
    }

    printf("mode=par n=%d threads=%d repeat=%d time_ms=%.3f det=%.6Lf\n",
           n, nthreads, repeat, best_ms, det);

    /* --------- Vérification contre la version séquentielle ----- */
    if (verify) {
        double t0 = now_ms();
        long double ref = det_seq(A, n);
        double t1 = now_ms();
        printf("verify_seq n=%d time_ms=%.3f det=%.6Lf match=%s\n",
               n, t1 - t0, ref,
               (fabsl(ref - det) < 1e-6L) ? "YES" : "NO");
    }

    if (n >= 10) {
        fprintf(stderr,
            "Info : n=%d => %llu permutations (%.3e). Temps O(n·n!).\n",
            n, FACT[n], (double)FACT[n]);
    }

    free(A);
    return 0;
}