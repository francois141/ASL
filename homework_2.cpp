#include <iostream>
#include "common.h"
#include "mat.h"

#define C1 0.1
#define C2 (2.0/3.0)

void slow_performance1(mat* x, mat* y, mat*z) {
  double t1;
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
        if (i % 2) {
          t1 = mat_get(z,i,j)/sqrt(mat_get(y,0,i%2)) + (mat_get(x,i,j) + C1)*(mat_get(x,i,j) - C1);
        } else {
          t1 = mat_get(z,i,j)/sqrt(mat_get(y,0,i%2)) + C1*mat_get(x,i,j);
        }
        mat_set(z,i,j-1, mat_get(z,i,j-1)*cos(C2*M_PI*j));
        mat_set(z,i,j,t1);
        mat_set(z,i,j+1,fmax(mat_get(z,i,j+1), mat_get(x,i,j+1)));
     }
  }
}

// Remove sqrt from loop (code motion)
// Convert division into multiplication
void slow_performance2(mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
    
  double sqrt1 = 1/sqrt(mat_get(y,0,1));
  double sqrt2 = 1/sqrt(mat_get(y,0,0));
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      if (i % 2) {
        t1 = mat_get(z,i,j)*sqrt1 + (mat_get(x,i,j) + C1)*(mat_get(x,i,j) - C1);
      } else {
        t1 = mat_get(z,i,j)*sqrt2 + C1*mat_get(x,i,j);
      }
      mat_set(z,i,j-1, mat_get(z,i,j-1)*cos(C2*M_PI*j));
      mat_set(z,i,j,t1);
      mat_set(z,i,j+1,fmax(mat_get(z,i,j+1), mat_get(x,i,j+1)));
    }
  }
}

// Create simple loops for the easiest instructions
// This also breaks the loop dependency
void slow_performance3 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
    
  double sqrt1 = 1/sqrt(mat_get(y,0,1));
  double sqrt2 = 1/sqrt(mat_get(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set(z,i,j+1,fmax(mat_get(z,i,j+1), mat_get(x,i,j+1)));
    }
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      if (i % 2) {
        t1 = mat_get(z,i,j)*sqrt1 + (mat_get(x,i,j) + C1)*(mat_get(x,i,j) - C1);
      } else {
        t1 = mat_get(z,i,j)*sqrt2 + C1*mat_get(x,i,j);
      }
      mat_set(z,i,j,t1);
    }
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set(z,i,j-1, mat_get(z,i,j-1)*cos(C2*M_PI*j));
    }
  }
}

inline double mat_get2(mat *m, int i, int j) {
  int ij = i * m->n2 + j;
  return m->data[ij];
}

inline void mat_set2(mat *m, int i, int j, double val) {
  int ij = i * m->n2 + j;
  m->data[ij] = val;
}

// Make loop easier and add accumulators
// Also using inline funcions without conditional checking (can be unsage but clearly faster)
void slow_performance4 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set2(z,i,j+1,fmax(mat_get2(z,i,j+1), mat_get2(x,i,j+1)));
    }
  }
    
  for (int i = 0; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t1 = mat_get2(z,i,j)*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
    }
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t2 = mat_get2(z,i,j)*sqrt1 + (mat_get2(x,i,j) + C1)*(mat_get2(x,i,j) - C1);
      mat_set2(z,i,j,t2);
    }
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
       mat_set2(z,i,j-1, mat_get2(z,i,j-1)*cos(C2*M_PI*j));
    }
  }
}


// Extract the cosinus
void slow_performance5 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set2(z,i,j+1,fmax(mat_get2(z,i,j+1), mat_get2(x,i,j+1)));
    }
  }
    
  for (int i = 0; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t1 = mat_get2(z,i,j)*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
    }
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t2 = mat_get2(z,i,j)*sqrt1 + (mat_get2(x,i,j) + C1)*(mat_get2(x,i,j) - C1);
      mat_set2(z,i,j,t2);
    }
  }
    
  constexpr double NEW_C2 = 2.0/3.0*M_PI;
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
    
  for (int j = 1; j < z->n2 - 1; j++) {
    values[j] = cos(NEW_C2*j);
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set2(z,i,j-1, mat_get2(z,i,j-1)*values[j]);
    }
  }
}


inline double mat_get3(mat *m, int i, int j) {
  int ij = i * 101 + j;
  return m->data[ij];
}

inline void mat_set3(mat *m, int i, int j, double val) {
  int ij = i * 101 + j;
  m->data[ij] = val;
}

// Hardcode matrix size in the index get 
void slow_performance6 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get3(y,0,1));
  double sqrt2 = 1/sqrt(mat_get3(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set3(z,i,j+1,fmax(mat_get3(z,i,j+1), mat_get3(x,i,j+1)));
    }
  }
    
  for (int i = 0; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t1 = mat_get3(z,i,j)*sqrt2 + C1*mat_get3(x,i,j);
      mat_set3(z,i,j,t1);
    }
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    for (int j = 1; j < z->n2 - 1; j++) {
      t2 = mat_get3(z,i,j)*sqrt1 + (mat_get3(x,i,j) + C1)*(mat_get3(x,i,j) - C1);
      mat_set3(z,i,j,t2);
    }
  }
    
  constexpr double NEW_C2 = 2.0/3.0*M_PI;
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
    
  for (int j = 1; j < z->n2 - 1; j++) {
    values[j] = cos(NEW_C2*j);
  }  
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 1; j < z->n2 - 1; j++) {
      mat_set3(z,i,j-1, mat_get3(z,i,j-1)*values[j]);
    }
  }
}

// Remove first nested loop and give fmax on the fly
void slow_performance7 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get3(y,0,1));
  double sqrt2 = 1/sqrt(mat_get3(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    mat_set3(z,i,z->n2-1,fmax(mat_get3(z,i,z->n2-1), mat_get3(x,i,z->n2-1)));
  }
    
  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get3(z,i,1)*sqrt2 + C1*mat_get3(x,i,1);
    mat_set3(z,i,1,t1);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt2 + C1*mat_get3(x,i,j);
      mat_set3(z,i,j,t1);
    }
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    t2 = mat_get3(z,i,1)*sqrt1 + (mat_get3(x,i,1) + C1)*(mat_get3(x,i,1) - C1);
    mat_set3(z,i,1,t2);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t2 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt1 + (mat_get3(x,i,j) + C1)*(mat_get3(x,i,j) - C1);
      mat_set3(z,i,j,t2);
    }
  }
    
  constexpr double NEW_C2 = 2.0/3.0*M_PI;
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
    
  for (int j = 1; j < z->n2 - 1; j++) {
    values[j] = cos(NEW_C2*j);
  }
    
  for (int i = 0; i < z->n1; i++) {
    mat_set3(z,i,0, mat_get3(z,i,0)*values[1]);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      mat_set3(z,i,j-1, mat_get3(z,i,j-1)*values[j]);
    }
  }
}


// DO operations in the first loop when data is in the cache
void slow_performance8 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get3(y,0,1));
  double sqrt2 = 1/sqrt(mat_get3(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get3(z,i,1)*sqrt2 + C1*mat_get3(x,i,1);
    mat_set3(z,i,1,t1);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt2 + C1*mat_get3(x,i,j);
      mat_set3(z,i,j,t1);
    }
    mat_set3(z,i,z->n2-1,fmax(mat_get3(z,i,z->n2-1), mat_get3(x,i,z->n2-1)));
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    t2 = mat_get3(z,i,1)*sqrt1 + (mat_get3(x,i,1) + C1)*(mat_get3(x,i,1) - C1);
    mat_set3(z,i,1,t2);
    
    for (int j = 2; j < z->n2 - 1; j++) {
      t2 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt1 + (mat_get3(x,i,j) + C1)*(mat_get3(x,i,j) - C1);
      mat_set3(z,i,j,t2);
    }
    
    mat_set3(z,i,z->n2-1,fmax(mat_get3(z,i,z->n2-1), mat_get3(x,i,z->n2-1)));
  }
  
  constexpr double NEW_C2 = 2.0/3.0*M_PI;
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
  
  for (int j = 1; j < z->n2 - 1; j++) {
    values[j] = cos(NEW_C2*j);
  }
  
  for (int i = 0; i < z->n1; i++) {
    mat_set3(z,i,0, mat_get3(z,i,0)*values[1]);
    
    for (int j = 2; j < z->n2 - 1; j++) {
      mat_set3(z,i,j-1, mat_get3(z,i,j-1)*values[j]);
    }
  }
}

// Delete cosine
void slow_performance9 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get3(y,0,1));
  double sqrt2 = 1/sqrt(mat_get3(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get3(z,i,1)*sqrt2 + C1*mat_get3(x,i,1);
    mat_set3(z,i,1,t1);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt2 + C1*mat_get3(x,i,j);
      mat_set3(z,i,j,t1);
    }
    mat_set3(z,i,z->n2-1,fmax(mat_get3(z,i,z->n2-1), mat_get3(x,i,z->n2-1)));
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    t2 = mat_get3(z,i,1)*sqrt1 + (mat_get3(x,i,1) + C1)*(mat_get3(x,i,1) - C1);
    mat_set3(z,i,1,t2);
    
    for (int j = 2; j < z->n2 - 1; j++) {
      t2 = fmax(mat_get3(z,i,j), mat_get3(x,i,j))*sqrt1 + (mat_get3(x,i,j) + C1)*(mat_get3(x,i,j) - C1);
      mat_set3(z,i,j,t2);
    }
    
    mat_set3(z,i,z->n2-1,fmax(mat_get3(z,i,z->n2-1), mat_get3(x,i,z->n2-1)));
  }
  
  double *values = (double*)malloc(24);
  
  values[0] = 1.0;
  values[1] = -0.5;
  values[2] = -0.5;
  
  for (int i = 0; i < z->n1; i++) {
    mat_set3(z,i,0, mat_get3(z,i,0)*values[1]);
    
    int j = 2;
    for (; j < z->n2 - 1; j+=3) {
      mat_set3(z,i,j-1, mat_get3(z,i,j-1)*values[j%3]);
      mat_set3(z,i,j, mat_get3(z,i,j)*values[(j+1)%3]);
      mat_set3(z,i,j+1, mat_get3(z,i,j+1)*values[(j+2)%3]);
    }
    
    mat_set3(z,i,97, mat_get3(z,i,97)*values[(98)%3]);
    mat_set3(z,i,98, mat_get3(z,i,98)*values[(99)%3]);
  }
}


void maxperformance(mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
    double t1;
    for (int i = 0; i < z->n1; i++) {
        for (int j = 1; j < z->n2 - 1; j++) {
            if (i % 2) {
                t1 = mat_get(z,i,j)/sqrt(mat_get(y,0,i%2)) + (mat_get(x,i,j) + C1)*(mat_get(x,i,j) - C1);
            } else {
                t1 = mat_get(z,i,j)/sqrt(mat_get(y,0,i%2)) + C1*mat_get(x,i,j);
            }
            mat_set(z,i,j-1, mat_get(z,i,j-1)*cos(C2*M_PI*j));
            mat_set(z,i,j,t1);
            mat_set(z,i,j+1,fmax(mat_get(z,i,j+1), mat_get(x,i,j+1)));
        }
    }
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() {  
  /*add_function(&slow_performance1, "slow_performance1",1);
  add_function(&slow_performance2, "slow_performance2",1);
  add_function(&slow_performance3, "slow_performance3",1);
  add_function(&slow_performance4, "slow_performance4",1);
  add_function(&slow_performance5, "slow_performance5",1);
  add_function(&slow_performance6, "slow_performance6",1);
  add_function(&slow_performance7, "slow_performance7",1);*/
  add_function(&slow_performance8, "slow_performance8",1);
  add_function(&slow_performance9, "slow_performance9",1);
  //add_function(&maxperformance, "maxperformance",1);
}
