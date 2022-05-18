#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
// gcc -Ofast -funroll-loops -lm CH1.c -o CH1.out

typedef uint_fast8_t uint8;
typedef int_fast8_t int8;
typedef uint32_t uint32;

struct pixel {
    uint8 *__restrict x;
    uint8 *__restrict y;
    uint8 *__restrict z;
};

//struct pixel32 {
//    uint32 *__restrict x;
//    uint32 *__restrict y;
//    uint32 *__restrict z;
//};

typedef struct {
    uint32 *__restrict x;
    uint32 xlen;
    uint32 *__restrict y;
    uint32 ylen;
    uint32 *__restrict z;
    uint32 zlen;
} vec3;

typedef struct {
    uint32 *__restrict x;
    uint32 xlen;
    uint32 *__restrict xy;
    uint32 xylen;
    uint32 *__restrict xz;
    uint32 xzlen;
    uint32 *__restrict xyz;
    uint32 xyzlen;
    uint32 *__restrict y;
    uint32 ylen;
    uint32 *__restrict yz;
    uint32 yzlen;
    uint32 *__restrict z;
    uint32 zlen;
} vecKernel;

uint32 seed;

inline uint32 myRandom() {
    seed = (214013 * seed + 2531011);
    return (seed >> 13);
}

void allocatePixel(struct pixel *__restrict image, uint32 N) {
    image->x = (uint8 *)malloc(N * sizeof(uint8));
    image->y = (uint8 *)malloc(N * sizeof(uint8));
    image->z = (uint8 *)malloc(N * sizeof(uint8));
}

//void allocatePixel32(struct pixel32 *__restrict image, uint32 N) {
//    image->x = (uint32 *)malloc(N * sizeof(uint32));
//    image->y = (uint32 *)malloc(N * sizeof(uint32));
//    image->z = (uint32 *)malloc(N * sizeof(uint32));
//}

void allocateVec(vec3 *i, uint32 N) {
    i->x = (uint32 *)malloc(N * sizeof(uint32));
    i->y = (uint32 *)malloc(N * sizeof(uint32));
    i->z = (uint32 *)malloc(N * sizeof(uint32));
}

void allocateKernel(vecKernel *i, uint32 N) {
    i->x = (uint32 *)malloc(N * sizeof(uint32));
    i->y = (uint32 *)malloc(N * sizeof(uint32));
    i->z = (uint32 *)malloc(N * sizeof(uint32));
    i->xy = (uint32 *)malloc(N * sizeof(uint32));
    i->xz = (uint32 *)malloc(N * sizeof(uint32));
    i->yz = (uint32 *)malloc(N * sizeof(uint32));
    i->xyz = (uint32 *)malloc(N * sizeof(uint32));
}

void InitImage(struct pixel I, uint32 N) {
    for(uint32 j = 0; j < N; j++) {
        I.x[j] = myRandom() % 31; // 5 bits
        I.y[j] = myRandom() % 31;
        I.z[j] = myRandom() % 31;
    }
}

vecKernel *InitKernel(uint32 Klen, uint32 N) {
    vecKernel *K = (vecKernel *)malloc(sizeof(vecKernel));
    allocateKernel(K, Klen / 7 + 7); // magic
    K->xlen = 0, K->ylen = 0, K->zlen = 0;
    K->xyzlen = 0, K->xylen = 0, K->xzlen = 0;
    K->yzlen = 0;
    for(uint32 i = 0; i < Klen; i++) {
        uint8 rand = myRandom() & 7; // 8 bits
        if(i >= N) break;
        if(!rand) continue;
        /*
        x 100 4; 4,5,6,7
        y 010 2; 2,3,6,7
        z 001 1; 1,3,5,7
        */
        switch(rand) {
            case 1: K->z[K->zlen++] = i; break;
            case 2: K->y[K->ylen++] = i; break;
            case 3: K->yz[K->yzlen++] = i; break;
            case 4: K->x[K->xlen++] = i; break;
            case 5: K->xz[K->xzlen++] = i; break;
            case 6: K->xy[K->xylen++] = i; break;
            case 7: K->xyz[K->xyzlen++] = i; break;
        }
    }
    return K;
}

inline void forloop1(uint32 *__restrict v,uint8 *__restrict T,uint32 *K, uint32 len, uint32 N) {
    for(uint32 k = 0; k < len; k++) {
        for(uint32 i = 0; i + K[k] < N; i++) {
            v[i] += T[i + K[k]];
        }
    }
}

inline void forloop2(uint32 *__restrict v, uint32 *__restrict v2,
                     uint8 *__restrict T, uint8 *__restrict T2,
                     uint32 *K, uint32 len, uint32 N) {
    for(uint32 k = 0; k < len; k++) {
        for(uint32 i = 0; i + K[k] < N; i++) {
            v[i] += T[i + K[k]];
            v2[i] += T2[i + K[k]];
        }
    }
}

inline void forloop3(uint32 *__restrict v, uint32 *__restrict v2, uint32 *__restrict v3,
                     uint8 *__restrict T, uint8 *__restrict T2, uint8 *__restrict T3,
                     uint32 *K, uint32 len, uint32 N) {
    for(uint32 k = 0; k < len; k++) {
        for(uint32 i = 0; i + K[k] < N; i++) {
            v[i] += T[i + K[k]];
            v2[i] += T2[i + K[k]];
            v3[i] += T3[i + K[k]];
        }
    }
}

void TransfImage(struct pixel I, uint32 N, vecKernel *K) {
    vec3 v;
    {
        allocateVec(&v, N);
        struct pixel T;
        allocatePixel(&T, N);
        for(uint32 i = 0; i < N; i++) {
            T.x[i] = I.x[i];
            T.y[i] = I.y[i];
            T.z[i] = I.z[i];
            v.x[i] = 0;
            v.y[i] = 0;
            v.z[i] = 0;
        }
        forloop1(v.x, T.x, K->x, K->xlen, N);
        forloop1(v.y, T.y, K->y, K->ylen, N);
        forloop1(v.z, T.z, K->z, K->zlen, N);
        forloop2(v.x, v.y, T.x, T.y, K->xy, K->xylen, N);
        forloop2(v.x, v.z, T.x, T.z, K->xz, K->xzlen, N);
        forloop2(v.y, v.z, T.y, T.z, K->yz, K->yzlen, N);
        forloop3(v.x, v.y, v.z, T.x, T.y, T.z, K->xyz, K->xyzlen, N);
        free(T.x);
        free(T.y);
        free(T.z);
    }

    for(uint32 i = 0; i < N; i++) {
        const uint8 x = v.x[i], y = v.y[i], z = v.z[i];
        const float c = 31.0 / (x + y + z + 1);
        I.x[i] = (uint8)(x * c);
        I.y[i] = (uint8)(y * c);
        I.z[i] = (uint8)(z * c);
    }
    free(v.x);
    free(v.y);
    free(v.z);
}

int main(int argc, char **argv) {
    uint32 i, sumX = 0, sumY = 0, sumZ = 0, N = 10000, Klen = 2000, Rep = 1000;

    seed = 12345;

    if(argc > 1) { N = atoi(argv[1]); }
    if(argc > 2) { Klen = atoi(argv[2]); }
    if(argc > 3) { Rep = atoi(argv[3]); }

    printf("Build time %s\n", __TIMESTAMP__);
    printf("Challenge #1: Vector size is %d. Kernel size is %d. Repeat %d times\n", N, Klen, Rep);

    struct pixel image;
    allocatePixel(&image, N);
    InitImage(image, N);
    vecKernel *K = InitKernel(Klen, N);

    for(i = 0; i < Rep; i++) {
        TransfImage(image, N, K);
        int ii;
        ii = myRandom() % N; image.x[ii] = myRandom() % 31;
        ii = myRandom() % N; image.y[ii] = myRandom() % 31;
        ii = myRandom() % N; image.z[ii] = myRandom() % 31;
    }

    for(i = 0, sumX = 0; i < N; i++) {
        sumX += image.x[i];
        sumY += image.y[i];
        sumZ += image.z[i];
    }

    printf("Result: sumX= %d, sumY= %d, sumZ= %d\n", sumX, sumY, sumZ);
    printf("%s\n", (sumX + sumY + sumZ) == 294167 ? "OK" : "ERROR");
    // Free K ...
    return 0;
}
