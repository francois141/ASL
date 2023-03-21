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
  
  __m256d one_mask = _mm256_castsi256_pd(_mm256_set1_epi32(-1));

  for (int i = 0; i < n; i += 16) 
  {
    {      
      // Step 1 : Get values for k 
      __m256d y_curr = _mm256_load_pd((double*)&y[i]);
  
      // Step 2 : Get re and im
      __m256d buff1 = _mm256_load_pd((double*)&x[i]);
      __m256d buff2 = _mm256_load_pd((double*)&x[i+2]);
          
      __m256d re = _mm256_unpacklo_pd(buff1, buff2);
      __m256d im = _mm256_unpackhi_pd(buff1, buff2);
          
      re = _mm256_permute4x64_pd(re, 0b11011000);
      im = _mm256_permute4x64_pd(im, 0b11011000);
          
      __m256d im2 = _mm256_mul_pd(im,im);
      __m256d reim = _mm256_mul_pd(re,im);
          
      __m256d newRe = _mm256_fmsub_pd(re,re, im2);
      __m256d newIm = _mm256_add_pd(reim, reim); 
        
      // Step 3 : Get the three vectors
      __m256d v1 = _mm256_min_pd(newRe,newIm);
      __m256d v2 = _mm256_max_pd(newRe,newIm);
          
      // max vector computation - no need to compute the sqrt
      __m256d v3 = _mm256_fmadd_pd(re,re, im2);
          
      // Step 4 : Compute final vector with the three vectors and the values of k
      __m256d c1 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.25, 0.25, 0.25, 0.25), _CMP_LT_OQ);
      __m256d c2 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.5, 0.5, 0.5, 0.5), _CMP_LT_OQ);
  
      __m256d mask1 = c1;
      __m256d mask2 = _mm256_andnot_pd(c1,c2);
      __m256d mask3 = _mm256_xor_pd(c2,one_mask);
  
      // TODO : Improve bitwise operations
      __m256d final = _mm256_and_pd(v1,mask1);
      __m256d final2 = _mm256_an''d_pd(v2,mask2);
      __m256d final3 = _mm256_and_pd(v3,mask3);
      
      final = _mm256_or_pd(final,final2);
      final = _mm256_or_pd(final,final3);
      
      // Step 5 : Increment the current vector
      final = _mm256_add_pd(final,y_curr);
      _mm256_store_pd(&y[i],final);
    }
    {      
      // Step 1 : Get values for k 
      __m256d y_curr = _mm256_load_pd((double*)&y[i+4]);
  
      // Step 2 : Get re and im
      __m256d buff1 = _mm256_load_pd((double*)&x[i+4]);
      __m256d buff2 = _mm256_load_pd((double*)&x[i+6]);
          
      __m256d re = _mm256_unpacklo_pd(buff1, buff2);
      __m256d im = _mm256_unpackhi_pd(buff1, buff2);
          
      re = _mm256_permute4x64_pd(re, 0b11011000);
      im = _mm256_permute4x64_pd(im, 0b11011000);
          
      __m256d im2 = _mm256_mul_pd(im,im);
      __m256d reim = _mm256_mul_pd(re,im);
          
      __m256d newRe = _mm256_fmsub_pd(re,re, im2);
      __m256d newIm = _mm256_add_pd(reim, reim); 
        
      // Step 3 : Get the three vectors
      __m256d v1 = _mm256_min_pd(newRe,newIm);
      __m256d v2 = _mm256_max_pd(newRe,newIm);
          
      // max vector computation - no need to compute the sqrt
      __m256d v3 = _mm256_fmadd_pd(re,re, im2);
          
      // Step 4 : Compute final vector with the three vectors and the values of k
      __m256d c1 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.25, 0.25, 0.25, 0.25), _CMP_LT_OQ);
      __m256d c2 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.5, 0.5, 0.5, 0.5), _CMP_LT_OQ);
  
      __m256d mask1 = c1;
      __m256d mask2 = _mm256_andnot_pd(c1,c2);
      __m256d mask3 = _mm256_xor_pd(c2,one_mask);
  
      // TODO : Improve bitwise operations
      __m256d final = _mm256_and_pd(v1,mask1);
      __m256d final2 = _mm256_and_pd(v2,mask2);
      __m256d final3 = _mm256_and_pd(v3,mask3);
      
      final = _mm256_or_pd(final,final2);
      final = _mm256_or_pd(final,final3);
      
      // Step 5 : Increment the current vector
      final = _mm256_add_pd(final,y_curr);
      _mm256_store_pd(&y[i+4],final);
    }
    {      
      // Step 1 : Get values for k 
      __m256d y_curr = _mm256_load_pd((double*)&y[i+8]);
  
      // Step 2 : Get re and im
      __m256d buff1 = _mm256_load_pd((double*)&x[i+8]);
      __m256d buff2 = _mm256_load_pd((double*)&x[i+10]);
          
      __m256d re = _mm256_unpacklo_pd(buff1, buff2);
      __m256d im = _mm256_unpackhi_pd(buff1, buff2);
          
      re = _mm256_permute4x64_pd(re, 0b11011000);
      im = _mm256_permute4x64_pd(im, 0b11011000);
          
      __m256d im2 = _mm256_mul_pd(im,im);
      __m256d reim = _mm256_mul_pd(re,im);
          
      __m256d newRe = _mm256_fmsub_pd(re,re, im2);
      __m256d newIm = _mm256_add_pd(reim, reim); 
        
      // Step 3 : Get the three vectors
      __m256d v1 = _mm256_min_pd(newRe,newIm);
      __m256d v2 = _mm256_max_pd(newRe,newIm);
          
      // max vector computation - no need to compute the sqrt
      __m256d v3 = _mm256_fmadd_pd(re,re, im2);
          
      // Step 4 : Compute final vector with the three vectors and the values of k
      __m256d c1 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.25, 0.25, 0.25, 0.25), _CMP_LT_OQ);
      __m256d c2 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.5, 0.5, 0.5, 0.5), _CMP_LT_OQ);
  
      __m256d mask1 = c1;
      __m256d mask2 = _mm256_andnot_pd(c1,c2);
      __m256d mask3 = _mm256_xor_pd(c2,one_mask);
  
      // TODO : Improve bitwise operations
      __m256d final = _mm256_and_pd(v1,mask1);
      __m256d final2 = _mm256_and_pd(v2,mask2);
      __m256d final3 = _mm256_and_pd(v3,mask3);
      
      final = _mm256_or_pd(final,final2);
      final = _mm256_or_pd(final,final3);
      
      // Step 5 : Increment the current vector
      final = _mm256_add_pd(final,y_curr);
      _mm256_store_pd(&y[i+8],final);
    }
    {      
      // Step 1 : Get values for k 
      __m256d y_curr = _mm256_load_pd((double*)&y[i+12]);
  
      // Step 2 : Get re and im
      __m256d buff1 = _mm256_load_pd((double*)&x[i+12]);
      __m256d buff2 = _mm256_load_pd((double*)&x[i+14]);
          
      __m256d re = _mm256_unpacklo_pd(buff1, buff2);
      __m256d im = _mm256_unpackhi_pd(buff1, buff2);
          
      re = _mm256_permute4x64_pd(re, 0b11011000);
      im = _mm256_permute4x64_pd(im, 0b11011000);
          
      __m256d im2 = _mm256_mul_pd(im,im);
      __m256d reim = _mm256_mul_pd(re,im);
          
      __m256d newRe = _mm256_fmsub_pd(re,re, im2);
      __m256d newIm = _mm256_add_pd(reim, reim); 
        
      // Step 3 : Get the three vectors
      __m256d v1 = _mm256_min_pd(newRe,newIm);
      __m256d v2 = _mm256_max_pd(newRe,newIm);
          
      // max vector computation - no need to compute the sqrt
      __m256d v3 = _mm256_fmadd_pd(re,re, im2);
          
      // Step 4 : Compute final vector with the three vectors and the values of k
      __m256d c1 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.25, 0.25, 0.25, 0.25), _CMP_LT_OQ);
      __m256d c2 = _mm256_cmp_pd(y_curr, _mm256_set_pd(0.5, 0.5, 0.5, 0.5), _CMP_LT_OQ);
  
      __m256d mask1 = c1;
      __m256d mask2 = _mm256_andnot_pd(c1,c2);
      __m256d mask3 = _mm256_xor_pd(c2,one_mask);
  
      // TODO : Improve bitwise operations
      __m256d final = _mm256_and_pd(v1,mask1);
      __m256d final2 = _mm256_and_pd(v2,mask2);
      __m256d final3 = _mm256_and_pd(v3,mask3);
      
      final = _mm256_or_pd(final,final2);
      final = _mm256_or_pd(final,final3);
      
      // Step 5 : Increment the current vector
      final = _mm256_add_pd(final,y_curr);
      _mm256_store_pd(&y[i+12],final);
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
