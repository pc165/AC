#include <stdio.h>
#include <stdlib.h>

// Define here as constant for easy change
#define REAL double

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Function to generate pseudo-random numbers
int myRandom() {
  seed = (214013*seed+2531011);
  return (seed>>3);
}

void printCheck ( REAL V[], int N )
{
  int x;

  REAL S=0;
  for (x=0; x<=N+1; x++)
    S = S + V[x];

  printf("\nCheckSum = %1.10e\n", S);

  //for (x=0; x<8; x++)
  //  printf("%1.10e\n", V[x]);
  for (x=0; x<=8; x++)
    printf("%1.10e  ", V[x*N/8]);
}


void SimulationStep ( REAL *In, REAL *Out, REAL L, int N )
{
  for (int x=1; x<=N; x++)
    Out[x] = In[x] - 2.0f*L*In[x] + L*In[x+1] + L*In[x-1];
}

REAL RandomFactor ( int K )
{
  int *List = malloc ( sizeof(*List)*K );
  for (int x=0; x<K; x++)
    List[x] = myRandom();

  int max= List[0];
  for (int x=1; x<K; x++)
    if (List[x]>max) max = List[x];

  int min= List[0];
  for (int x=1; x<K; x++)
    if (List[x]<min) min = List[x];

  free(List);
  //printf("max = %d, min = %d\n", max, min);

  return (REAL)(max-min)/(REAL)(max+min);
}

void CopyAndMultiply ( REAL *In, REAL *Out, int N, REAL F )
{
  for (int x=1; x<=N; x++)
    Out[x] = F*In[x];
}


int main(int argc, char **argv)
{
  int  x, t, N= 10000000, T=1000;
  REAL L= 0.123456789;
  REAL *U1, *U2;

  seed = 12345;

  if (argc>1) { T = atoi(argv[1]); }   // get  first command line parameter
  if (argc>2) { N = atoi(argv[2]); }   // get second command line parameter
  if (argc>3) { L = atof(argv[3]); }   // get  third command line parameter
  if (argc>4) { seed = atof(argv[4]);} // get fourth command line parameter
 
  if (N < 1 || L >= 0.5) {
    printf("arguments: T N L (T: steps, N: vector size, L < 0.5)\n");
    return 1;
  }

  U1 = malloc ( sizeof(*U1)*(N+2) );
  U2 = malloc ( sizeof(*U2)*(N+2) );
  if (!U1 || !U2) { printf("Cannot allocate vectors\n"); exit(1); }

  // initialize temperatures at time t=0  
  REAL t1, t2, t3;
  t1= 0.1234; t2 = -0.9456; t3 = 0.9789;
  for (x=0; x<=N+1; x++)
  {
    U1[x] = t1*t2*t3;
    t1= t1>1.3?   0.1234: t1+0.00567;
    t2= t2>1.8?  -0.9456: t2+0.00987;
    t3= t3<0.007? 0.9789: t3-0.00321;
  }
 
  // initialize fixed boundary conditions on U1
  {
    U1[0]  = 1.2345678901;
    U1[N+1]= 1.2345678901;
  }

  printf("Challenge #4: Simulate %d steps on 1-D vector of %d elements with L=%1.10e\n", T, N, L);

  for (t=1; t<=T; t++)
  {  // loop on time
    SimulationStep ( U1, U2, L, N );
    REAL Factor = RandomFactor(N/10);
    CopyAndMultiply( U2, U1, N, Factor );
    //printf("F = %1.10e\n", Factor);
  }

  printCheck(U1,N);
  free(U1); free(U2);
}
