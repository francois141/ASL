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
  __m256d X_lower_vector,X_higher_vector;
  
  __m256d BASE_H = _mm256_set_pd(4,3,2,1);
  
  __m256d H0 = _mm256_broadcast_sd(&h[0]);
  __m256d H1 = _mm256_broadcast_sd(&h[1]);
  __m256d H2 = _mm256_broadcast_sd(&h[2]);
  __m256d H3 = _mm256_broadcast_sd(&h[3]);
  
  X_lower_vector = abs(_mm256_load_pd(&x[0]));
  for(int i = 0; i < N - 3; i += 4) {
    X_higher_vector = abs(_mm256_load_pd(&x[i+4]));
    __m256d csloc = BASE_H;

    __m256d shuffle_0 = X_lower_vector;
    __m256d shuffle_1 = _mm256_permute4x64_pd(_mm256_blend_pd(X_lower_vector, X_higher_vector, 0b0001), 0b00111001);
    __m256d shuffle_2 = _mm256_permute2f128_pd(X_lower_vector, X_higher_vector, 0b00100001);
    __m256d shuffle_3 = _mm256_permute4x64_pd(_mm256_blend_pd(X_lower_vector, X_higher_vector, 0b0111), 0b10010011);
    
    __m256d yval = _mm256_mul_pd(_mm256_mul_pd(csloc, H0), shuffle_3);
    yval = _mm256_fmadd_pd(_mm256_add_pd(csloc, _mm256_set_pd(1,1,1,1)), _mm256_mul_pd(shuffle_2,H1), yval);
    yval = _mm256_fmadd_pd(_mm256_add_pd(csloc, _mm256_set_pd(2,2,2,2)), _mm256_mul_pd(shuffle_1,H2), yval);
    yval = _mm256_fmadd_pd(_mm256_add_pd(csloc, _mm256_set_pd(3,3,3,3)), _mm256_mul_pd(shuffle_0,H3), yval);
    
    _mm256_store_pd(y + i, yval);
    
    BASE_H = _mm256_add_pd(BASE_H, _mm256_set_pd(4,4,4,4));
    
    X_lower_vector = X_higher_vector;
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
