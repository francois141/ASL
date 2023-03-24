#include "common.h"
#include <immintrin.h>
#include <math.h>
#include <iostream>

using namespace std;

void slow_performance1(double *x, double* h, double* y, int N, int M) {
  for (int i = 0; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += (i + k + 1) * h[k] * fabs(x[i + (M - 1) - k]);
    }
  }

}

inline __m256d abs(__m256d c) {
  return _mm256_andnot_pd(_mm256_set_pd(-0.0,-0.0,-0.0,-0.0),c); 
}

void maxperformance(double *x, double* h, double* y, int N, int M) {
  
  __m256d H0 = _mm256_mul_pd(_mm256_set_pd(4,3,2,1),_mm256_broadcast_sd(&h[0]));
  __m256d H1 = _mm256_mul_pd(_mm256_set_pd(5,4,3,2),_mm256_broadcast_sd(&h[1]));
  __m256d H2 = _mm256_mul_pd(_mm256_set_pd(6,5,4,3),_mm256_broadcast_sd(&h[2]));
  __m256d H3 = _mm256_mul_pd(_mm256_set_pd(7,6,5,4),_mm256_broadcast_sd(&h[3]));
  
  __m256d H0_inc =  _mm256_mul_pd(_mm256_set_pd(4,4,4,4),_mm256_broadcast_sd(&h[0]));
  __m256d H1_inc =  _mm256_mul_pd(_mm256_set_pd(4,4,4,4),_mm256_broadcast_sd(&h[1]));
  __m256d H2_inc =  _mm256_mul_pd(_mm256_set_pd(4,4,4,4),_mm256_broadcast_sd(&h[2]));
  __m256d H3_inc =  _mm256_mul_pd(_mm256_set_pd(4,4,4,4),_mm256_broadcast_sd(&h[3]));
  
  __m256d x0,x1,x2,x3;
  
  for(int i = 0; i < N - 3; i += 4) {
    // Load x
    x0 = abs(_mm256_load_pd(&x[i]));
    x1 = abs(_mm256_load_pd(&x[i+1]));
    x2 = abs(_mm256_load_pd(&x[i+2]));
    x3 = abs(_mm256_load_pd(&x[i+3]));
    
    // Apply fma
    __m256d tmp1 = _mm256_fmadd_pd(H1,x2,_mm256_mul_pd(H0,x3));
    __m256d tmp2 = _mm256_fmadd_pd(H3,x0,_mm256_mul_pd(H2,x1));
    
    // Store the value in y 
    _mm256_store_pd(&y[i], _mm256_add_pd(tmp1,tmp2));
    
    // Update the H vectors
    H0 = _mm256_add_pd(H0,H0_inc);
    H1 = _mm256_add_pd(H1,H1_inc);
    H2 = _mm256_add_pd(H2,H2_inc);
    H3 = _mm256_add_pd(H3,H3_inc); 
  }
  
  return;
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() 
{
  add_function(&slow_performance1, "slow_performance1",1);
  add_function(&maxperformance, "maxperformance",1);
}
