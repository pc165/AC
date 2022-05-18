#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Pablo Chen 1530349
// gcc -Ofast -funroll-loops CH2.c
// 0,17608 +- 0,00440 seconds time elapsed  ( +-  2,50% )


#define LEN 16                          // Max number of digits in a ull_t
#define MAX 10000000000000000           // 10E16
#define MAXINT 100000000                // 10E8
#define CEIL(A) (((A) + LEN - 1) / LEN) // Minimum array size required to hold N digits

// Variable used to generate pseudo-random numbers
unsigned int seed;
typedef unsigned long long ull_t;
const ull_t POWER[LEN + 1] =
    {
        1,
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        100000000,
        1000000000,
        10000000000,
        100000000000,
        1000000000000,
        10000000000000,
        100000000000000,
        1000000000000000,
        10000000000000000,
};

// Function to generate pseudo-random numbers
int myRandom() {
    seed = (214013 * seed + 2531011);
    return (seed >> 13);
}

void LongNumPrint(ull_t *a, unsigned N, char b[]) {
    printf("%s", b);
    for (int i = CEIL(N) - 1; i >= 0; i--) {
        printf("%016llu", a[i]); // pad numbers with zeros to print them with a fixed width
    }
    printf("\n");
}

void LongNumInit(ull_t *vec, const unsigned N) {
    for (int i = 0; i < N / LEN; i++) {
        ull_t c = 0;
        for (int z = 0; z < LEN; z++) {
            ull_t e = POWER[z];
            int s = (myRandom() % 10);
            c += s * e;
        }
        vec[i] = c;
    }

    if (N % LEN) {
        ull_t c = 0;
        for (int z = 0; z < N % LEN; z++) {
            ull_t e = POWER[z];
            int s = (myRandom() % 10);
            c += s * e;
        }
        vec[CEIL(N) - 1] = c;
    }

    return;
}

void LongNumTruncate(ull_t *a, const unsigned N) {
    if (N % LEN) {
        a[CEIL(N) - 1] = a[CEIL(N) - 1] % POWER[N % LEN];
    }
}

void LongNumAddDigit(ull_t *a, const char digit, const unsigned N) {
    ull_t R = a[0] + digit;
    if (R < MAX) {
        a[0] = R;
        LongNumTruncate(a, N);
        return;
    }

    a[0] = R - MAX;
    char carry = 1;
    int i = 1;

    while (carry && i < CEIL(N)) {
        a[i] += carry;
        carry = a[i] / MAX;
        a[i] %= MAX;
    }

    LongNumTruncate(a, N);

    return;
}

void LongNumAddition(ull_t *a, ull_t *b, ull_t *c, const unsigned N) {
    char carry = 0;
    for (int z = 0; z < CEIL(N); z++) {
        c[z] = a[z] + b[z] + carry;
        carry = c[z] / MAX;
        c[z] %= MAX;
    }

    LongNumTruncate(c, N);

    return;
}

void LongNumHorizAdd(ull_t *a, ull_t *out, const unsigned N) {
    memset(out, 0, sizeof(ull_t) * CEIL(N));
    for (int i = 0; i < CEIL(N); i++) {
        ull_t t = a[i];
        while (t) {
            out[0] += t % 10;
            t /= 10;
        }
    }

    LongNumTruncate(out, N);
}

void LongNumMultiply64Bits(ull_t v1, ull_t v2, ull_t carry, ull_t *hi, ull_t *lo) {
    /*  
        product = (v1 * v2) mod 2^64 =  high * 2^64 + low = high << 64 | low
        
        product / MAX = (high * 2^64 + low) / MAX
        product / MAX = high * (ULLONG_MAX / MAX) + low / MAX

        high16 = high * (ULLONG_MAX / MAX) + low / MAX
        low16 = product - high16 * MAX 
    */
    __int128_t product = (__int128_t)v1 * v2;
    ull_t high = (ull_t)(product >> 64);
    ull_t low = (ull_t)product;

    // get high and low in base MAX
    ull_t high16 = (ull_t)(high * ((double)ULLONG_MAX / MAX)) + low / MAX;
    ull_t low16 = (ull_t)((product - (__int128_t)high16 * MAX));

    // Fix overflow
    high16 -= low16 >> 63;
    low16 = (ull_t)((product - (__int128_t)high16 * MAX));

    *hi = high16;
    *lo += low16 + carry;

    *hi += *lo / MAX;
    *lo %= MAX;

    // check
    // if (((__int128_t)high16 * MAX + low16) != product) {
    //     __int128_t w = (__int128_t)high16 * MAX;
    //     printf("%llu%llu\n", high16, low16);
    // }
}

void LongNumMultiply(ull_t *Vin1, ull_t *Vin2, ull_t *VoutH, ull_t *VoutL, ull_t *Res, const unsigned N) {
    size_t s = sizeof(ull_t) * CEIL(N);
    memset(Res, 0, 2 * s); // clear Res
    for (int i = 0; i < CEIL(N); i++) {
        ull_t carry = 0;
        for (int z = 0; z < CEIL(N); z++) {
            LongNumMultiply64Bits(Vin1[i], Vin2[z], carry, &carry, &Res[i + z]);
        }
        Res[i + CEIL(N)] += carry;
    }

    // Split the array in two vectors
    if (N < LEN) {
        VoutL[0] = Res[0] % POWER[N];
        VoutH[0] = Res[0] / POWER[N] + Res[1] * POWER[LEN - N];
    } else if (CEIL(N) == N / LEN) { // for all N which are powers of 2
        memcpy(VoutL, Res, s);
        memcpy(VoutH, Res + CEIL(N), s);
    } else { // other cases
        memset(VoutL, 0, s);
        memset(VoutH, 0, s);
        memcpy(VoutL, Res, sizeof(ull_t) * (N / LEN));
        memcpy(VoutH, Res + CEIL(N), sizeof(ull_t) * (N / LEN));
        int len = N - LEN * (N / LEN);
        ull_t midL = Res[N / LEN] % POWER[len];
        ull_t midH = Res[N / LEN] / POWER[len];
        VoutL[N / LEN] = midL;
        while (len < LEN) { // shift left vector LEN - len positions
            char carry = 0;
            for (int z = 0; z < CEIL(N); z++) {
                ull_t r = VoutH[z] * 10 + carry;
                carry = r / MAX;
                r = r % MAX;
                VoutH[z] = r;
            }
            len++;
        }
        VoutH[0] += midH;
    }
    return;
}

int main(char argc, char **argv) {
    int i;
    int N = 10000, Rep = 50;

    if (argc > 1) {
        N = atoi(argv[1]);
    }

    if (argc > 2) {
        Rep = atoi(argv[2]);
    }

    printf("Challenge #2: Vector size is %d. Repeat %d times\n", N, Rep);

    seed = 12345;

    size_t s = sizeof(ull_t) * CEIL(N);
    ull_t *V1 = malloc(s * 8);
    memset(V1, 0, s * 8);
    ull_t *V2 = V1 + CEIL(N);
    ull_t *V3 = V2 + CEIL(N);
    ull_t *V4 = V3 + CEIL(N);
    ull_t *Res = V4 + 2 * CEIL(N);
    LongNumInit(V1, N);
    LongNumInit(V2, N);
    LongNumInit(V3, N);

    for (int i = 0; i < Rep; i++) {
        LongNumAddition(V1, V2, V4, N);
        LongNumMultiply(V3, V4, V2, V1, Res, N);
        LongNumHorizAdd(V1, V2, N);
        LongNumAddDigit(V3, V2[0] % 10, N);
    }

    LongNumPrint(V1, 32, "V1:");
    LongNumPrint(V2, 32, "V2:");
    LongNumPrint(V3, 32, "V3:");
    LongNumPrint(V4, 32, "V4:");

    free(V1);
    return 0;
}