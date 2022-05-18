#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Pablo Chen 1530349
// The main idea is to increase cache hits by reusing the calculations.
// 2.657412179 seconds time elapsed with float
// 3.533459400 seconds time elapsed with double
// module load gcc/9.2.0
// gcc -Ofast -fstrict-aliasing -fopenmp -funroll-loops -ffast-math -march=native -lm CH4.c -o CH4

#define CHUNK 1990           // This number must satisfy (CHUNK/2) % 10 == 5
#define TIMESTEP (CHUNK / 2) // Number of time steps that will be calculated each time simulation is executed

// We need to use double type because there are values that are really near to 0.
// If we use float we can get more than 0.01 % error, e.g T = 9927 N = 1234412
// float Checksum = 1.5926073438e+05
// double CheckSum = 1.5921428362e+05
// error 0.0292%
typedef double real_t;

unsigned int seed;

inline int SIZE(int N) { return (N + 1) / (CHUNK + 1); }

inline real_t OP(real_t *in, real_t f, real_t l, real_t l2, int n) {
    return f * (l2 * in[n] + l * (in[n + 1] + in[n - 1]));
}

inline void SWAP(real_t **a, real_t **b) {
    real_t *tmp = *a;
    *a = *b;
    *b = tmp;
}

inline int myRandom() {
    seed = (214013 * seed + 2531011);
    return (seed >> 3);
}

void printCheck(real_t V[], int N) {
    int x;

    real_t S = 0;
#pragma omp parallel
#pragma omp for reduction(+ \
                          : S)
    for (x = 0; x <= N + 1; x++)
        S = S + V[x];

    printf("\nCheckSum = %1.10e\n", S);

    for (x = 0; x <= 8; x++)
        printf("%1.10e  ", V[x * N / 8]);
    printf("\n");
}

real_t RandomFactor(int K) {
    int max, min;
    max = myRandom();
    min = max;
    for (int x = 1; x < K; x++) {
        int c = myRandom();
        if (c > max)
            max = c;

        if (c < min)
            min = c;
    }

    return (real_t)(max - min) / (real_t)(max + min);
}

#define EQ(a, b) (fabs((a) - (b)) < __DBL_EPSILON__)

void triangleUp(real_t *out, real_t *in, real_t *left, real_t *right, real_t *factor, real_t l, int nStart) {
    real_t l2 = 1 - 2.0f * l;

    for (int t = 0; t < TIMESTEP; t++) {
        int start = nStart + t, end = nStart + CHUNK - t - 1;
        left[t] = in[start];
        right[t] = in[end];
        real_t f = factor[t];
        for (int n = start; n <= end; n++) {
            out[n] = OP(in, f, l, l2, n);
        }

        out[start - 1] = in[start - 1];
        out[end + 1] = in[end + 1];
        SWAP(&in, &out);
    }
}

void triangleDown(real_t *out, real_t *in, real_t *left, real_t *right, real_t *factor, real_t l, int nStart) {
    real_t l2 = 1 - 2.0f * l;
    for (int t = 0; t < TIMESTEP; t++) {
        int start = nStart - t, end = nStart + t;
        real_t f = factor[t];
        in[start - 1] = left[t];
        in[end + 1] = right[t];
        for (int n = start; n <= end; n++) {
            out[n] = OP(in, f, l, l2, n);
        }

        SWAP(&in, &out);
    }
}

void Prologue(real_t *out, real_t *in, real_t *left, real_t *factor, real_t l, int nStart) {
    SWAP(&in, &out);
    real_t l2 = 1 - 2.0f * l;
    for (int t = 0; t < TIMESTEP - 1; t++) {
        real_t f = factor[t + 1];
        int end = nStart + t;
        in[end + 1] = left[t + 1];
        for (int n = nStart; n <= end; n++) {
            out[n] = OP(in, f, l, l2, n);
        }
        SWAP(&in, &out);
    }
}

void Epilogue(real_t *out, real_t *in, real_t *right, real_t *factor, real_t l, int nEnd) {
    real_t l2 = 1 - 2.0f * l;

    // Check if last point of the last triangle overlaps with the last position
    int c = (SIZE(nEnd) - 1) * (CHUNK + 1) + 1 + CHUNK;
    int b = c > nEnd; // overlaps, we need to calculate "TIMESTEP - 1" steps

    real_t tmp = in[nEnd - (TIMESTEP - b - 1) - 1];
    int start;
    for (int t = 0; t < TIMESTEP - b; t++) {
        real_t f = factor[t];
        start = nEnd - t;
        in[start - 1] = right[t + b];
        for (int n = start; n <= nEnd; n++) {
            out[n] = OP(in, f, l, l2, n);
        }
        SWAP(&in, &out);
    }

    if (b) {
        memcpy(out + start, in + start, (nEnd - start + 1) * sizeof(real_t));
        out[start - 1] = tmp;
    } else {
        SWAP(&in, &out);
    }
}

void Epilogue2(real_t *out, real_t *in, real_t *right, real_t *factor, real_t l, int nStart, const int nEnd) {
    real_t l2 = 1 - 2.0f * l;
    for (int t = 0; t < TIMESTEP; t++) {
        real_t f = factor[t];
        int start = nStart - t;
        in[start - 1] = right[t];
        for (int n = start; n <= nEnd; n++) {
            out[n] = OP(in, f, l, l2, n);
        }
        SWAP(&in, &out);
    }
}

void Simulation(real_t *out, real_t *in, real_t *factor, real_t l, int N, real_t *left, real_t *right) {

#pragma omp parallel num_threads(4)
#pragma omp for schedule(static, 8)
    for (int i = 0; i < SIZE(N); i++) {
        int nStart = i * (CHUNK + 1) + 1;
        triangleUp(out, in, &left[i * TIMESTEP], &right[i * TIMESTEP], factor, l, nStart);
    }

#pragma omp parallel num_threads(4)
#pragma omp for schedule(static, 8)
    for (int i = 1; i < SIZE(N); i++) {
        int nStart = i * (CHUNK + 1);
        triangleDown(out, in, &right[(i - 1) * TIMESTEP], &left[i * TIMESTEP], factor, l, nStart);
    }

    // coordinate of the last point of the last triangle
    int a = (SIZE(N) - 1) * (CHUNK + 1) + 1 + CHUNK;

    if (a >= N) { // check if the last position (N) was calculated
        Epilogue(out, in, right + (SIZE(N) - 1) * TIMESTEP, factor, l, N);
    } else {
        Epilogue2(out, in, right + (SIZE(N) - 1) * TIMESTEP, factor, l, a, N);
    }

    Prologue(out, in, left, factor, l, 1);
}

int main(int argc, char **argv) {
    int N = 10000000, T = 1000;

    real_t *out, *in, l = 0.123456789;

    seed = 12345;

    if (argc > 1) T = atoi(argv[1]);
    if (argc > 2) N = atoi(argv[2]);
    if (argc > 3) l = atof(argv[3]);
    if (argc > 4) seed = atof(argv[4]);

    printf("Challenge #4: Simulate %d steps on 1-D vector of %d elements with L=%1.10e\n", T, N, l);

    if (N < 1 || l >= 0.5) {
        printf("arguments: T N L (T: steps, N: vector size, L < 0.5)\n");
        return 1;
    }

    out = malloc((N + 2) * sizeof(real_t) * 2);
    in = out + (N + 2);

    if (!out) {
        printf("Cannot allocate vectors\n");
        exit(1);
    }

    {
        real_t t1, t2, t3;
        t1 = 0.1234;
        t2 = -0.9456;
        t3 = 0.9789;
        for (int x = 0; x <= N + 1; x++) {
            in[x] = t1 * t2 * t3;
            out[x] = t1 * t2 * t3;
            t1 = t1 > 1.3 ? 0.1234 : t1 + 0.00567;
            t2 = t2 > 1.8 ? -0.9456 : t2 + 0.00987;
            t3 = t3 < 0.007 ? 0.9789 : t3 - 0.00321;
        }

        in[0] = 1.2345678901;
        in[N + 1] = 1.2345678901;

        out[0] = 1.2345678901;
        out[N + 1] = 1.2345678901;
    }

    int G = N <= CHUNK ? T : T % TIMESTEP;

    if (G) {
        real_t l2 = 1 - 2.0f * l;
        for (int t = 0; t < G; t++) {
            real_t f = RandomFactor(N / 10);
            for (int n = 1; n <= N; n++) {
                real_t nextT = l2 * in[n] + l * (in[n + 1] + in[n - 1]);
                out[n] = nextT * f;
            }
            memcpy(in, out, (N + 2) * sizeof(real_t));
        }

        T = T - G;
    }

    if (T > 0) {
        real_t *factor = malloc((T + TIMESTEP * SIZE(N) * 2) * sizeof(real_t)); // factor + (right + left)
        real_t *left = factor + T;
        real_t *right = left + (TIMESTEP * SIZE(N));

        if (!factor) {
            printf("Cannot allocate vectors\n");
            exit(1);
        }

        for (int i = 0; i < T; i++)
            factor[i] = RandomFactor(N / 10);

        for (int i = 0; i < T / TIMESTEP; i++) {
            Simulation(out, in, factor + i * TIMESTEP, l, N, left, right);
            memcpy(in, out, (N + 2) * sizeof(real_t)); // in and out must have the same data
        }

        free(factor);
    }

    printCheck(out, N);

    free(out);
}
