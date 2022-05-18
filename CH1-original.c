#include <stdio.h>
#include <stdlib.h>

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Function to generate pseudo-random numbers
inline int myRandom() {
  seed = (214013*seed+2531011);
  return (seed>>13);
}

void InitKernel( int *k, unsigned Klen )
{
  for ( int i=0; i< Klen; i++ ) {
    k[i] = myRandom() & 7;
  }
}

struct pixel { unsigned char x; unsigned char y; unsigned char z; };

void InitImage( struct pixel *I, unsigned N )
{
  for ( int j=0; j< N; j++ ) 
  {
    I[j].x = myRandom() % 31;
    I[j].y = myRandom() % 31;
    I[j].z = myRandom() % 31;
  }
}

void TransfImage( struct pixel * I, unsigned N, const int *K, unsigned Klen )
{
  struct pixel *T= (struct pixel*) malloc( N*sizeof(struct pixel) );

  // copy I to T in order to prevent data dependencies
  for ( int i=0; i< N; i++ )  T[i].x = I[i].x;
  for ( int i=0; i< N; i++ )  T[i].y = I[i].y;
  for ( int i=0; i< N; i++ )  T[i].z = I[i].z;

  for ( int i=0; i< N; i++ ) 
  {
    unsigned vx=0, vy=0, vz=0, sum;
    for ( int k=0; k < Klen; k++ )
    {
      int image_i = i + k;
      if (image_i >=0 && image_i < N)
      {
        if (K[k] >= 4)
          vx = vx + T[ image_i ].x;
        if ((K[k] & 3) == 2 || (K[k] & 3) == 3) 
          vy = vy + T[ image_i ].y;
        if (K[k] == 1 || K[k] == 3 || K[k] == 5 || K[k] == 7) 
          vz = vz + T[ image_i ].z;
      }
    }

    I[i].x = vx; 
    I[i].y = vy; 
    I[i].z = vz; 
    sum   = I[i].x+I[i].y+I[i].z+1;
    I[i].x = I[i].x*31 / sum;
    I[i].y = I[i].y*31 / sum;
    I[i].z = I[i].z*31 / sum;
  }
  free(T);
}

int main (int argc, char **argv)
{
  int i, sumX, sumY, sumZ, N=10000, Klen=2000, Rep=1000;

  seed = 12345;

  // obtain parameters at run time
  if (argc>1) { N    = atoi(argv[1]); }
  if (argc>2) { Klen = atoi(argv[2]); }
  if (argc>3) { Rep  = atoi(argv[3]); }

  printf("Challenge #1: Vector size is %d. Kernel size is %d. Repeat %d times\n", N, Klen, Rep);

  // Create Image & Kernel
  struct pixel *image= (struct pixel*) malloc( N*sizeof(struct pixel) );
  int *K = (int *) malloc( Klen*sizeof(int) );
 
  InitImage ( image, N );
  InitKernel( K, Klen );

  // Repeat
  for (i=0; i<Rep; i++)
  {
    TransfImage( image, N, K, Klen ); 
    int ii;
    ii = myRandom() % N;  image[ii].x = myRandom() % 31;
    ii = myRandom() % N;  image[ii].y = myRandom() % 31;
    ii = myRandom() % N;  image[ii].z = myRandom() % 31;
  }

  for (i=0, sumX=0; i<N; i++) sumX += image[i].x;
  for (i=0, sumY=0; i<N; i++) sumY += image[i].y;
  for (i=0, sumZ=0; i<N; i++) sumZ += image[i].z;
    
  printf("Result: sumX= %d, sumY= %d, sumZ= %d\n", sumX, sumY, sumZ);

  free(image); free(K);

  return 0;
}
