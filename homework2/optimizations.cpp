#include "common.h"
#include "mat.h"

#define C1 0.1
#define C2 (2.0/3.0)
#define C12 0.01

#include <iostream>

using namespace std;

void base_function(mat* x, mat* y, mat*z) {
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
void optimization_1(mat* x, mat* y, mat* z) {
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
void optimization_2 (mat* x, mat* y, mat* z) {
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
void optimization_3 (mat* x, mat* y, mat* z) {
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
void optimization_4 (mat* x, mat* y, mat* z) {
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
    
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
    
  for (int j = 0; j < z->n2 - 1; j++) {
    values[j] = cos(C2*M_PI*(j+1));
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 0; j < z->n2 - 2; j++) {
      mat_set2(z,i,j, mat_get2(z,i,j)*values[j]);
    }
  }
}

// Remove first nested loop and give fmax on the fly
void optimization_5 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));
    
  for (int i = 0; i < z->n1; i++) {
    mat_set2(z,i,z->n2-1,fmax(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
  }
    
  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
    }
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    t2 = mat_get2(z,i,1)*sqrt1 + (mat_get2(x,i,1) + C1)*(mat_get2(x,i,1) - C1);
    mat_set2(z,i,1,t2);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t2 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt1 + (mat_get2(x,i,j) + C1)*(mat_get2(x,i,j) - C1);
      mat_set2(z,i,j,t2);
    }
  }
    
  double *values = (double*)malloc(sizeof(double)*z->n2 - 1);
    
  for (int j = 0; j < z->n2 - 1; j++) {
    values[j] = cos(C2*M_PI*(j+1));
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 0; j < z->n2 - 2; j++) {
      mat_set2(z,i,j, mat_get2(z,i,j)*values[j]);
    }
  }
}




// Delete cosine
void optimization_6 (mat* x, mat* y, mat* z) {
  /* Add your final implementation here */
  double t1;
  double t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
    }
    mat_set2(z,i,z->n2-1,fmax(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
  }
    
  for (int i = 1; i < z->n1; i+=2) {
    t2 = mat_get2(z,i,1)*sqrt1 + (mat_get2(x,i,1) + C1)*(mat_get2(x,i,1) - C1);
    mat_set2(z,i,1,t2);
    
    for (int j = 2; j < z->n2 - 1; j++) {
      t2 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt1 + (mat_get2(x,i,j) + C1)*(mat_get2(x,i,j) - C1);
      mat_set2(z,i,j,t2);
    }
    mat_set2(z,i,z->n2-1,fmax(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
  }
  
  for (int i = 0; i < z->n1; i++) {
    for (int j = 0; j < z->n2-2; j+=3) {
      mat_set2(z,i,j, mat_get2(z,i,j)*-0.5);
      mat_set2(z,i,j+1, mat_get2(z,i,j+1)*-0.5);
    }
  }
}

// Merge functions together
void optimization_7(mat* x, mat* y, mat* z) {
  double t1,t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
    
    t2 = mat_get2(z,i+1,1)*sqrt1 + (mat_get2(x,i+1,1) + C1)*(mat_get2(x,i+1,1) - C1);
    mat_set2(z,i+1,1,t2);
      
    for (int j = 2; j < z->n2 - 1; j++) {
      t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
      
      t2 = fmax(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
      mat_set2(z,i+1,j,t2);
    }
    
    mat_set2(z,i,z->n2-1,fmax(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
    mat_set2(z,i+1,z->n2-1,fmax(mat_get2(z,i+1,z->n2-1), mat_get2(x,i+1,z->n2-1)));
  }
    
  for (int i = 0; i < z->n1; i++) {
    for (int j = 0; j < z->n2-2; j+=3) {
      mat_set2(z,i,j, mat_get2(z,i,j)*-0.5);
      mat_set2(z,i,j+1, mat_get2(z,i,j+1)*-0.5);
    }
  }
}

// Unroll
void optimization_8(mat* x, mat* y, mat* z) {
  double t1,t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
    
    t2 = mat_get2(z,i+1,1)*sqrt1 + (mat_get2(x,i+1,1) + C1)*(mat_get2(x,i+1,1) - C1);
    mat_set2(z,i+1,1,t2);
    
    mat_set2(z,i,0, mat_get2(z,i,0)*-0.5);
    mat_set2(z,i,1, mat_get2(z,i,1)*-0.5);
    
    mat_set2(z,i+1,0, mat_get2(z,i+1,0)*-0.5);
    mat_set2(z,i+1,1, mat_get2(z,i+1,1)*-0.5);
      
    int j = 2;
    for (; j < z->n2 - 3; j+=3) {
      t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
      
      t2 = fmax(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
      mat_set2(z,i+1,j,t2);
      
      t1 = fmax(mat_get2(z,i,j+1), mat_get2(x,i,j+1))*sqrt2 + C1*mat_get2(x,i,j+1);
      mat_set2(z,i,j+1,-0.5*t1);
      
      t2 = fmax(mat_get2(z,i+1,j+1), mat_get2(x,i+1,j+1))*sqrt1 + (mat_get2(x,i+1,j+1) + C1)*(mat_get2(x,i+1,j+1) - C1);
      mat_set2(z,i+1,j+1,-0.5*t2);
      
      t1 = fmax(mat_get2(z,i,j+2), mat_get2(x,i,j+2))*sqrt2 + C1*mat_get2(x,i,j+2);
      mat_set2(z,i,j+2,-0.5*t1);
      
      t2 = fmax(mat_get2(z,i+1,j+2), mat_get2(x,i+1,j+2))*sqrt1 + (mat_get2(x,i+1,j+2) + C1)*(mat_get2(x,i+1,j+2) - C1);
      mat_set2(z,i+1,j+2,-0.5*t2);
    }
    
    t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
      
    t2 = fmax(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
    mat_set2(z,i+1,j,t2);
    
    j++;
      
    t1 = fmax(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
  
    t2 = fmax(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
    mat_set2(z,i+1,j,t2);
    
    mat_set2(z,i,z->n2-1,fmax(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
    mat_set2(z,i+1,z->n2-1,fmax(mat_get2(z,i+1,z->n2-1), mat_get2(x,i+1,z->n2-1)));
  }
}

// fmax --> max
void optimization_9(mat* x, mat* y, mat* z) {
  double t1,t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
    
    t2 = mat_get2(z,i+1,1)*sqrt1 + (mat_get2(x,i+1,1) + C1)*(mat_get2(x,i+1,1) - C1);
    mat_set2(z,i+1,1,t2);
    
    mat_set2(z,i,0, mat_get2(z,i,0)*-0.5);
    mat_set2(z,i,1, mat_get2(z,i,1)*-0.5);
    
    mat_set2(z,i+1,0, mat_get2(z,i+1,0)*-0.5);
    mat_set2(z,i+1,1, mat_get2(z,i+1,1)*-0.5);
      
    int j = 2;
    for (; j < z->n2 - 3; j+=3) {
      t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
      
      t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
      mat_set2(z,i+1,j,t2);
      
      t1 = max(mat_get2(z,i,j+1), mat_get2(x,i,j+1))*sqrt2 + C1*mat_get2(x,i,j+1);
      mat_set2(z,i,j+1,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+1), mat_get2(x,i+1,j+1))*sqrt1 + (mat_get2(x,i+1,j+1) + C1)*(mat_get2(x,i+1,j+1) - C1);
      mat_set2(z,i+1,j+1,-0.5*t2);
      
      t1 = max(mat_get2(z,i,j+2), mat_get2(x,i,j+2))*sqrt2 + C1*mat_get2(x,i,j+2);
      mat_set2(z,i,j+2,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+2), mat_get2(x,i+1,j+2))*sqrt1 + (mat_get2(x,i+1,j+2) + C1)*(mat_get2(x,i+1,j+2) - C1);
      mat_set2(z,i+1,j+2,-0.5*t2);
    }
    
    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
      
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
    mat_set2(z,i+1,j,t2);
    
    j++;
      
    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
  
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + (mat_get2(x,i+1,j) + C1)*(mat_get2(x,i+1,j) - C1);
    mat_set2(z,i+1,j,t2);
    
    mat_set2(z,i,z->n2-1,max(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
    mat_set2(z,i+1,z->n2-1,max(mat_get2(z,i+1,z->n2-1), mat_get2(x,i+1,z->n2-1)));
  }

}

// Transform algebraic expression
void optimization_10(mat* x, mat* y, mat* z) {
  double t1,t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));

  for (int i = 0; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
    
    t2 = mat_get2(z,i+1,1)*sqrt1 + mat_get2(x,i+1,1)*mat_get2(x,i+1,1) - C12;
    mat_set2(z,i+1,1,t2);
    
    mat_set2(z,i,0, mat_get2(z,i,0)*-0.5);
    mat_set2(z,i,1, mat_get2(z,i,1)*-0.5);
    
    mat_set2(z,i+1,0, mat_get2(z,i+1,0)*-0.5);
    mat_set2(z,i+1,1, mat_get2(z,i+1,1)*-0.5);
      
    int j = 2;
    for (; j < z->n2 - 3; j+=3) {
      t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
      
      t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
      mat_set2(z,i+1,j,t2);
      
      t1 = max(mat_get2(z,i,j+1), mat_get2(x,i,j+1))*sqrt2 + C1*mat_get2(x,i,j+1);
      mat_set2(z,i,j+1,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+1), mat_get2(x,i+1,j+1))*sqrt1 + mat_get2(x,i+1,j+1)*mat_get2(x,i+1,j+1) - C12;
      mat_set2(z,i+1,j+1,-0.5*t2);
      
      t1 = max(mat_get2(z,i,j+2), mat_get2(x,i,j+2))*sqrt2 + C1*mat_get2(x,i,j+2);
      mat_set2(z,i,j+2,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+2), mat_get2(x,i+1,j+2))*sqrt1 + mat_get2(x,i+1,j+2)*mat_get2(x,i+1,j+2) - C12;
      mat_set2(z,i+1,j+2,-0.5*t2);
    }
    
    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
      
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
    mat_set2(z,i+1,j,t2);
    
    j++;
      
    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
  
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
    mat_set2(z,i+1,j,t2);
    
    mat_set2(z,i,z->n2-1,max(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
    mat_set2(z,i+1,z->n2-1,max(mat_get2(z,i+1,z->n2-1), mat_get2(x,i+1,z->n2-1)));
  }
}

void maxperformance(mat* x, mat* y, mat* z) {
  double t1,t2;
    
  double sqrt1 = 1/sqrt(mat_get2(y,0,1));
  double sqrt2 = 1/sqrt(mat_get2(y,0,0));
  int i = 0;
  for (; i < z->n1; i+=2) {
    t1 = mat_get2(z,i,1)*sqrt2 + C1*mat_get2(x,i,1);
    mat_set2(z,i,1,t1);
    
    t2 = mat_get2(z,i+1,1)*sqrt1 + mat_get2(x,i+1,1)*mat_get2(x,i+1,1) - C12;
    mat_set2(z,i+1,1,t2);
    
    mat_set2(z,i,0, mat_get2(z,i,0)*-0.5);
    mat_set2(z,i,1, mat_get2(z,i,1)*-0.5);
    
    mat_set2(z,i+1,0, mat_get2(z,i+1,0)*-0.5);
    mat_set2(z,i+1,1, mat_get2(z,i+1,1)*-0.5);
      
    int j = 2;

    for (; j < z->n2 - 3; j+=3) {
      t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
      mat_set2(z,i,j,t1);
      
      t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
      mat_set2(z,i+1,j,t2);
      
      t1 = max(mat_get2(z,i,j+1), mat_get2(x,i,j+1))*sqrt2 + C1*mat_get2(x,i,j+1);
      mat_set2(z,i,j+1,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+1), mat_get2(x,i+1,j+1))*sqrt1 + mat_get2(x,i+1,j+1)*mat_get2(x,i+1,j+1) - C12;
      mat_set2(z,i+1,j+1,-0.5*t2);
      
      t1 = max(mat_get2(z,i,j+2), mat_get2(x,i,j+2))*sqrt2 + C1*mat_get2(x,i,j+2);
      mat_set2(z,i,j+2,-0.5*t1);
      
      t2 = max(mat_get2(z,i+1,j+2), mat_get2(x,i+1,j+2))*sqrt1 + mat_get2(x,i+1,j+2)*mat_get2(x,i+1,j+2) - C12;
      mat_set2(z,i+1,j+2,-0.5*t2);
    }
    

    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
      
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
    mat_set2(z,i+1,j,t2);
    
    j++;
      
    t1 = max(mat_get2(z,i,j), mat_get2(x,i,j))*sqrt2 + C1*mat_get2(x,i,j);
    mat_set2(z,i,j,t1);
  
    t2 = max(mat_get2(z,i+1,j), mat_get2(x,i+1,j))*sqrt1 + mat_get2(x,i+1,j)*mat_get2(x,i+1,j) - C12;
    mat_set2(z,i+1,j,t2);
    
    mat_set2(z,i,z->n2-1,max(mat_get2(z,i,z->n2-1), mat_get2(x,i,z->n2-1)));
    mat_set2(z,i+1,z->n2-1,max(mat_get2(z,i+1,z->n2-1), mat_get2(x,i+1,z->n2-1)));
  }

}



/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() {  
  add_function(&base_function, "base_function",1);
  add_function(&optimization_1, "optimization_1",1);
  add_function(&optimization_2, "optimization_2",1);
  add_function(&optimization_3, "optimization_3",1);
  add_function(&optimization_4, "optimization_4",1);
  add_function(&optimization_5, "optimization_5",1);
  add_function(&optimization_6, "optimization_6",1);
  add_function(&optimization_7, "optimization_7",1);
  add_function(&optimization_8, "optimization_8",1);
  add_function(&optimization_9, "optimization_9",1);
  add_function(&optimization_10, "optimization_10",1);
  add_function(&maxperformance, "maxperformance",1);
}
