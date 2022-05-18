#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Pablo Chen 1530349
// icc -Ofast -fno-alias -unroll4 -fp-model fast -march=ivybridge -mtune=ivybridge CH3.c -o CH3.out
// time 2.575969092 s

// Define here as constant for easy change
typedef float real_t;

void printCheck(real_t V[], int N, int T) {
    int x;

    real_t S = 0;
    for (x = 0; x <= N + 1; x++)
        S = S + V[x];

    printf("\nCheckSum = %1.10e\nSome values: ", S);

    for (x = 0; x < 10; x++)
        printf("(%d) = %1.10f, ", x * N / 10, V[x * N / 10]);

    printf("(%d) = %1.10f\n", x * N / 10, V[x * N / 10]);
}

/*
 Simutaltion of the heat equation using finite differences.
 The naive algorithm had a bad memory acces pattern, it looped over two large arrays (140 MB) 
 causing a lot of cache misses.
 In order to improve cache hit, one array was eliminated and the algorithm was changed to acces elements
 that we just calculated and were cached.
*/

void Simulation(real_t *__restrict U, real_t L, const unsigned int N, unsigned int T) {

    const real_t L2 = 1 - 2 * L;

    int rem = N % T;
    if (rem != 0) {
        rem = T;
        // Find the largest T' such that T' <= T  && N % T' == 0
        while (rem-- && N % rem)
            ;
        int end = (rem == N || rem == 1 ? T : T - rem);

        for (int t = 0; t < end; t++) {
            real_t prevj = U[0];
            for (int i = 1; i <= N; i++) {
                real_t nextT = L2 * U[i] + L * (U[i + 1] + prevj);
                prevj = U[i];
                U[i] = nextT;
            }
        }

        if (rem == N || rem == 1)
            return;
        T = rem;
    }

    real_t *tprev = (real_t *)malloc(sizeof(real_t) * T);

    const real_t t0 = U[0];
    for (int i = 0; i < T; i++) {
        tprev[i] = t0;
    }

    // start
    for (int t = 0; t < T; t++) {
        int start = 1, end = T - t;
        real_t prevj = tprev[t];
        tprev[t] = U[end]; // save last result for the next loop
        for (; start <= end; start++) {
            real_t thisT = U[start];
            real_t nextj = U[start + 1];
            real_t nextT = L2 * thisT + L * (nextj + prevj);
            prevj = thisT;
            U[start] = nextT;
        }
    }

    // mid, main loop
    real_t *aux2 = malloc(sizeof(real_t) * (T + 1));

    for (int z = T + 1; z < N; z += T) {
        for (int t = 0; t < T; t++) {
            int start = z - t, end = start + T - 1;
            // copy part of the vector so it can be vectorized
            memcpy(aux2 + 1, U + start, sizeof(real_t) * (T - 1) /*(end - start)*/);

            aux2[0] = tprev[t];
            tprev[t] = U[end]; // save last result for the next loop

            for (int i = 0; i <= end - start; i++) {
                real_t thisT = U[start + i];
                real_t nextj = U[start + 1 + i];
                real_t prevj = aux2[i];
                real_t nextT = L2 * thisT + L * (nextj + prevj);
                U[start + i] = nextT;
            }
        }
    }

    free(aux2);

    // end
    for (int t = 0; t < T - 1; t++) {
        int start = N - t, end = N;
        real_t prevj = tprev[t + 1];
        // no need to save last result
        for (; start <= end; start++) {
            real_t thisT = U[start];
            real_t nextj = U[start + 1];
            real_t nextT = L2 * thisT + L * (nextj + prevj);
            prevj = thisT;
            U[start] = nextT;
        }
    }

    free(tprev);
}

int main(int argc, char **argv) {
    int N = 10000000, T = 1000;
    real_t L = 0.123456;
    real_t *U1;

    if (argc > 1) {
        T = atoi(argv[1]);
    }
    if (argc > 2) {
        N = atoi(argv[2]);
    }
    if (argc > 3) {
        L = atof(argv[3]);
    }

    if (N < 1 || T < 1 || L >= 0.5) {
        printf("arguments: T N L (T: steps, N: vector size, L < 0.5)\n");
        return 1;
    }

    U1 = (real_t *)malloc(sizeof(real_t) * (N + 2));
    if (!U1) {
        printf("Cannot allocate vectors\n");
        exit(1);
    }

    // initialize temperatures at time t=0
    for (int x = 0; x <= N + 1; x++) {
        // U1[x] = x * x;
        U1[x] = x * 3.1416;
    }

    // initialize fixed boundary conditions on U
    U1[0] = 1.2345678e+12;
    U1[N + 1] = -1.2345678e+16;

    Simulation(U1, L, N, T);

    printf("Challenge #3: Simulate %d steps on 1-D vector of %d elements with L=%1.10e\n", T, N, L);
    printCheck(U1, N, T);
    free(U1);
}

