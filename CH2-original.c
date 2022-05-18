#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <expat_config

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Function to generate pseudo-random numbers
int myRandom() {
    seed = (214013 * seed + 2531011);
    return (seed >> 13);
}

void LongNumInit(char *L, unsigned N) {
    for (int i = 0; i < N; i++) {
        L[i] = myRandom() % 10; // digito decimal
    }
}

void LongNumPrint(char *L, char *__restrict L2, unsigned N, char *__restrict Name) {
    printf("%s:", Name);
    for (int i = N; i > 0; i--) {
        printf("%d", L[i - 1]);
        if (L[i - 1] != L2[N - i]) {
            printf("Error");
            exit(EXIT_FAILURE);
        }
    }

    printf("\nOK\n");
}

char LongNumAddition(char *Vin1, char *__restrict Vin2, char *__restrict Vout, unsigned N) {
    char CARRY = 0;
    for (int i = 0; i < N; i++) {
        char R = Vin1[i] + Vin2[i] + CARRY;
        if (R <= 9) {
            Vout[i] = R;
            CARRY = 0;
        } else {
            Vout[i] = R - 10;
            CARRY = 1;
        }
    }
    return CARRY;
}

char LongNumAddDigit(char *V, char digit, unsigned N) {
    int i = 0;
    char R = V[0] + digit;
    if (R < 10) {
        V[0] = R;
        return 0; // No carry
    }

    V[0] = R - 10;
    // add carry, maybe iteratively for all digits
    char CARRY = 1;
    i = 1;
    while (CARRY && i < N) {
        if (V[i] < 9) {
            V[i] = V[i] + 1;
            CARRY = 0;
        } else {
            V[i] = 0;
            i++; // CARRY remains set to 1
        }
    }
    return CARRY;
}

char LongNumHorizAdd(char *Vin, char *__restrict Vout, unsigned N) {
    char CARRY = 0;
    memset(Vout, 0, N);
    for (int i = 0; i < N; i++) {
        LongNumAddDigit(Vout, Vin[i], N);
    }
    return 0; // CARRY can never be set
}

char LongNumConstMult(char *V, unsigned N, char digit) {
    char CARRY = 0;
    for (int i = 0; i < N; i++) {
        char R = V[i] * digit + CARRY;
        CARRY = R / 10;
        R = R - CARRY * 10;
        V[i] = R;
    }
    return CARRY; // may be from 0 to 9
}

void LongNumMultiply(char *Vin1, char *__restrict Vin2, char *__restrict VoutH, char *__restrict VoutL, unsigned N) {

    // Create Temporal Long Integer with double size
    unsigned char *TEMP = (unsigned char *)malloc(2 * N * sizeof(unsigned char));
    unsigned char *RES = (unsigned char *)malloc(2 * N * sizeof(unsigned char));

    memset(RES, 0, 2 * N); // Set RES to 0

    for (int i = 0; i < N; i++) {
        memset(TEMP, 0, 2 * N);                 // Set TEMP to 0
        memcpy(TEMP + i, Vin1, N);              // Copy Vin1 -> TEMP, with offset i
        LongNumConstMult(TEMP, 2 * N, Vin2[i]); // TEMP * Vin2[i] -> TEMP
        LongNumAddition(TEMP, RES, RES, 2 * N); // TEMP + RES -> RES
    }

    // Result goes to VoutH-VoutL
    memcpy(VoutL, RES, N);     // Copy RES   -> VoutL
    memcpy(VoutH, RES + N, N); // Copy RES+N -> VoutH
}

int main(int argc, char **argv) {
    int i, sum1, sum2, sum3, N = 10000, Rep = 50;

    seed = 12345;

    // obtain parameters at run time
    if (argc > 1) {
        N = atoi(argv[1]);
    }
    if (argc > 2) {
        Rep = atoi(argv[2]);
    }

    printf("Challenge #2: Vector size is %d. Repeat %d times\n", N, Rep);

    // Create Long Nums
    unsigned char *__restrict V1 = (unsigned char *)malloc(N * sizeof(unsigned char));
    unsigned char *__restrict V2 = (unsigned char *)malloc(N * sizeof(unsigned char));
    unsigned char *__restrict V3 = (unsigned char *)malloc(N * sizeof(unsigned char));
    unsigned char *__restrict V4 = (unsigned char *)malloc(N * sizeof(unsigned char));

    LongNumInit(V1, N);
    LongNumInit(V2, N);
    LongNumInit(V3, N);

    // Repeat
    for (i = 0; i < Rep; i++) {
        LongNumAddition(V1, V2, V4, N);
        LongNumMultiply(V3, V4, V2, V1, N);
        LongNumHorizAdd(V1, V2, N);
        LongNumAddDigit(V3, V2[0], N);
    }

    // Print last 32 digits of Long Numbers
    char V1T[] = {3, 2, 1, 0, 3, 9, 9, 7, 7, 3, 1, 0, 2, 9, 8, 9, 5, 2, 5, 1, 7, 3, 9, 7, 6, 6, 1, 0, 8, 9, 8, 0},
         V2T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0, 8, 7},
         V3T[] = {6, 3, 3, 9, 2, 7, 7, 1, 3, 2, 4, 9, 5, 3, 8, 1, 8, 0, 8, 9, 0, 3, 8, 2, 8, 0, 6, 5, 7, 0, 9, 7},
         V4T[] = {2, 1, 9, 2, 9, 0, 1, 0, 5, 3, 8, 2, 2, 2, 3, 8, 3, 5, 0, 1, 9, 8, 2, 1, 0, 9, 9, 4, 1, 7, 2, 2};

    LongNumPrint(V1, V1T, 32, "V1");
    LongNumPrint(V2, V2T, 32, "V2");
    LongNumPrint(V3, V3T, 32, "V3");
    LongNumPrint(V4, V4T, 32, "V4");

    free(V1);
    free(V2);
    free(V3);
    free(V4);
    return 0;
}
