#include <immintrin.h>
#include "common.h"
#include "complex.h"
#include <iostream>
using namespace std;
void slow_performance1(quaternion_t x[N], quaternion_t y[N],quaternion_t A[N][N]) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
            A[i][j] = mul(x[i], y[j]);
        }
    }
}

static quaternion_t mul2(quaternion_t a, quaternion_t b) {
    quaternion_t r;
    r.r = a.r*b.r - a.i*b.i - a.j*b.j - a.k*b.k;
    r.i = a.r*b.i + a.i*b.r + a.j*b.k - a.k*b.j;
    r.j = a.r*b.j - a.i*b.k + a.j*b.r + a.k*b.i;
    r.k = a.r*b.k + a.i*b.j - a.j*b.i + a.k*b.r;
    return r;
}

void maxperformance(quaternion_t x[N], quaternion_t y[N],quaternion_t A[N][N]) {
  /* Add your final implementation here */
  for (int j = 0; j < N; j++) {
      __m256d v2 = _mm256_load_pd((double*)&y[j]);
      __m256d v4 = _mm256_permute_pd(v2,0b0101);
      __m256d v7 = _mm256_permute4x64_pd(v2,0b01001110);
      __m256d v11 = _mm256_permute4x64_pd(v2,0b00011011);
      
      for (int i = 0; i < N; i++) {
            // Step 1
            __m256d v1 = _mm256_broadcast_sd(&x[i].r);
            __m256d v3 = _mm256_broadcast_sd(&x[i].i);
            __m256d v6 = _mm256_broadcast_sd(&x[i].j);
            __m256d v10 = _mm256_broadcast_sd(&x[i].k);
            
            // Step 2
            __m256d v5 = _mm256_mul_pd(v3,v4);
            
            // Step 3
            __m256d v8 = _mm256_mul_pd(v6,v7);
            __m256d v9 = _mm256_set_pd(-1,1,1,-1);
            v9 = _mm256_mul_pd(v9,v8);
            
            // Step 4
            __m256d v12 = _mm256_mul_pd(v10,v11);
            __m256d v13 = _mm256_set_pd(1,1,-1,-1);
            v13 = _mm256_mul_pd(v13,v12);
            
            // Merge together
            __m256d result = _mm256_mul_pd(v1,v2);
            result = _mm256_addsub_pd(result,v5);
            result = _mm256_add_pd(result,v9);
            result = _mm256_add_pd(result,v13);
            
            // Store result
            _mm256_store_pd((double*)&A[i][j],result);
        }
    }
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() {
  add_function(&slow_performance1, "slow_performance1",1);
  add_function(&maxperformance, "maxperformance",1);
} 
