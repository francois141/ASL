#include <immintrin.h>
#include "common.h"

#include <iostream>

using namespace std;

// Precondition: 0 <= y[i] < 1.0
void slow_performance1(complex_t *x, double *y, int n) {
    for (int i = 0; i < n; i++) {
        unsigned int k = floor(4.0*y[i]);
        switch (k) {
            case 0:  y[i] += fmin(re(sqr(x[i])), im(sqr(x[i]))); break;
            case 1:  y[i] += fmax(re(sqr(x[i])), im(sqr(x[i]))); break;
            default: y[i] += pow(abs(x[i]), 2.0); break;
        }
    }
}

void maxperformance(complex_t* x, double* y, int n) {

  __m256d comp_mask1 = _mm256_set_pd(0.5,0.5,0.5,0.5);
  __m256d comp_mask2 = _mm256_set_pd(0.25,0.25,0.25,0.25);

  for (int i = 0; i < n; i += 16) {
  __m256d first_x = _mm256_load_pd((double*)&x[i]); 
  __m256d second_x = _mm256_load_pd((double*)&x[i + 2]); 
  __m256d first_x_2 = _mm256_load_pd((double*)&x[i + 4]); 
  __m256d second_x_2 = _mm256_load_pd((double*)&x[i + 6]); 
  __m256d first_x_3 = _mm256_load_pd((double*)&x[i + 8]); 
  __m256d second_x_3 = _mm256_load_pd((double*)&x[i + 10]); 
  __m256d first_x_4 = _mm256_load_pd((double*)&x[i + 12]); 
  __m256d second_x_4 = _mm256_load_pd((double*)&x[i + 14]); 
  {
      __m256d y_base = _mm256_load_pd(y + i);
  
      __m256d real_part = _mm256_unpacklo_pd(first_x, second_x);
      __m256d imaginary_part = _mm256_unpackhi_pd(first_x, second_x);
  
      real_part = _mm256_permute4x64_pd(real_part, 0b11011000);
      imaginary_part = _mm256_permute4x64_pd(imaginary_part, 0b11011000);
      
      __m256d real_mult_imaginary = _mm256_mul_pd(real_part, imaginary_part);
      real_mult_imaginary = _mm256_add_pd(real_mult_imaginary, real_mult_imaginary);
  
      __m256d imaginary_square = _mm256_mul_pd(imaginary_part, imaginary_part);
      
      __m256d re_mult_sub = _mm256_fmsub_pd(real_part, real_part, imaginary_square);
      __m256d re_mult_add = _mm256_fmadd_pd(real_part, real_part, imaginary_square);
  
      __m256d min_values = _mm256_min_pd(re_mult_sub, real_mult_imaginary);
      __m256d max_values = _mm256_max_pd(re_mult_sub, real_mult_imaginary);
  
      __m256d mask1 = _mm256_cmp_pd(y_base, comp_mask1, _CMP_GE_OQ);
      __m256d mask2 = _mm256_cmp_pd(y_base, comp_mask2, _CMP_GE_OQ);
  
      __m256d y_to_add = _mm256_blendv_pd(min_values, max_values, mask2);
      y_to_add = _mm256_blendv_pd(y_to_add, re_mult_add, mask1);
      
      _mm256_store_pd(y + i, _mm256_add_pd(y_base, y_to_add));
  }
  {
      __m256d y_base = _mm256_load_pd(y + i + 4);
  
      __m256d real_part = _mm256_unpacklo_pd(first_x_2, second_x_2);
      __m256d imaginary_part = _mm256_unpackhi_pd(first_x_2, second_x_2);
  
      real_part = _mm256_permute4x64_pd(real_part, 0b11011000);
      imaginary_part = _mm256_permute4x64_pd(imaginary_part, 0b11011000);
      
      __m256d real_mult_imaginary = _mm256_mul_pd(real_part, imaginary_part);
      real_mult_imaginary = _mm256_add_pd(real_mult_imaginary, real_mult_imaginary);
  
      __m256d imaginary_square = _mm256_mul_pd(imaginary_part, imaginary_part);
      
      __m256d re_mult_sub = _mm256_fmsub_pd(real_part, real_part, imaginary_square);
      __m256d re_mult_add = _mm256_fmadd_pd(real_part, real_part, imaginary_square);
  
      __m256d min_values = _mm256_min_pd(re_mult_sub, real_mult_imaginary);
      __m256d max_values = _mm256_max_pd(re_mult_sub, real_mult_imaginary);
  
      __m256d mask1 = _mm256_cmp_pd(y_base, comp_mask1, _CMP_GE_OQ);
      __m256d mask2 = _mm256_cmp_pd(y_base, comp_mask2, _CMP_GE_OQ);
  
      __m256d y_to_add = _mm256_blendv_pd(min_values, max_values, mask2);
      y_to_add = _mm256_blendv_pd(y_to_add, re_mult_add, mask1);
      
      _mm256_store_pd(y + i + 4, _mm256_add_pd(y_base, y_to_add));
  }
  {
      __m256d y_base = _mm256_load_pd(y + i + 8);
  
      __m256d real_part = _mm256_unpacklo_pd(first_x_3, second_x_3);
      __m256d imaginary_part = _mm256_unpackhi_pd(first_x_3, second_x_3);
  
      real_part = _mm256_permute4x64_pd(real_part, 0b11011000);
      imaginary_part = _mm256_permute4x64_pd(imaginary_part, 0b11011000);
      
      __m256d real_mult_imaginary = _mm256_mul_pd(real_part, imaginary_part);
      real_mult_imaginary = _mm256_add_pd(real_mult_imaginary, real_mult_imaginary);
  
      __m256d imaginary_square = _mm256_mul_pd(imaginary_part, imaginary_part);
      
      __m256d re_mult_sub = _mm256_fmsub_pd(real_part, real_part, imaginary_square);
      __m256d re_mult_add = _mm256_fmadd_pd(real_part, real_part, imaginary_square);
  
      __m256d min_values = _mm256_min_pd(re_mult_sub, real_mult_imaginary);
      __m256d max_values = _mm256_max_pd(re_mult_sub, real_mult_imaginary);
  
      __m256d mask1 = _mm256_cmp_pd(y_base, comp_mask1, _CMP_GE_OQ);
      __m256d mask2 = _mm256_cmp_pd(y_base, comp_mask2, _CMP_GE_OQ);
  
      __m256d y_to_add = _mm256_blendv_pd(min_values, max_values, mask2);
      y_to_add = _mm256_blendv_pd(y_to_add, re_mult_add, mask1);
      
      _mm256_store_pd(y + i + 8, _mm256_add_pd(y_base, y_to_add));
  }
  {
      __m256d y_base = _mm256_load_pd(y + i + 12);
  
      __m256d real_part = _mm256_unpacklo_pd(first_x_4, second_x_4);
      __m256d imaginary_part = _mm256_unpackhi_pd(first_x_4, second_x_4);
  
      real_part = _mm256_permute4x64_pd(real_part, 0b11011000);
      imaginary_part = _mm256_permute4x64_pd(imaginary_part, 0b11011000);
      
      __m256d real_mult_imaginary = _mm256_mul_pd(real_part, imaginary_part);
      real_mult_imaginary = _mm256_add_pd(real_mult_imaginary, real_mult_imaginary);
  
      __m256d imaginary_square = _mm256_mul_pd(imaginary_part, imaginary_part);
      
      __m256d re_mult_sub = _mm256_fmsub_pd(real_part, real_part, imaginary_square);
      __m256d re_mult_add = _mm256_fmadd_pd(real_part, real_part, imaginary_square);
  
      __m256d min_values = _mm256_min_pd(re_mult_sub, real_mult_imaginary);
      __m256d max_values = _mm256_max_pd(re_mult_sub, real_mult_imaginary);
  
      __m256d mask1 = _mm256_cmp_pd(y_base, comp_mask1, _CMP_GE_OQ);
      __m256d mask2 = _mm256_cmp_pd(y_base, comp_mask2, _CMP_GE_OQ);
  
      __m256d y_to_add = _mm256_blendv_pd(min_values, max_values, mask2);
      y_to_add = _mm256_blendv_pd(y_to_add, re_mult_add, mask1);
      
      _mm256_store_pd(y + i + 12, _mm256_add_pd(y_base, y_to_add));
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
