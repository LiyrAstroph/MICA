#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "proto.h"
#include "allvars.h"

void test()
{
  double *A, *B, *C;
  int m=3, n=2, k=4;
  A=array_malloc(m*n);
  B=array_malloc(n*k);
  C=array_malloc(m*k);

  A[0*n+0] = 1.0;
  A[0*n+1] = 2.0;

  A[1*n+0] = 1.2;
  A[1*n+1] = 3.3;

  A[2*n+0] = 2.5;
  A[2*n+1] = 1.6;

  B[0*k + 0 ] = 0.2;
  B[0*k + 1 ] = 0.4;
  B[0*k + 2 ] = 0.5;
  B[0*k + 3 ] = 0.3;

  B[1*k + 0 ] = -0.2;
  B[1*k + 1 ] = -0.24;
  B[1*k + 2 ] = 0.5;
  B[1*k + 3 ] = -0.3;

  multiply_mat_MN(A, B, C, m, k, n);
  display_mat(C, m, k);
  printf("\n");
  multiply_mat_MN_transposeA(A, C, B, n, k, m);
  display_mat(B, n, k);
  printf("\n");
  multiply_mat_MN_transposeB(B, C, A, n, m, k);
  display_mat(A, n, m);  
}