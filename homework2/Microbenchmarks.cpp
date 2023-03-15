#include "include/microbenchmark.h"
#include "include/tsc_x86.h"
#include "include/foo.h"
#include <iostream>
#include <climits>

using namespace std;

#define DOUBLE_HIGH_LATENCY 467859236578264.33333458878685748889999878
#define DOUBLE_HIGH_LATENCY_2 46785923643525264.423543523
#define DOUBLE_LOW_LATENCY 1.0

#define FILL_ARRAY_HIGH_LATENCY for(int i = 0; i < n;i++) { arr[i] = DOUBLE_HIGH_LATENCY; }
#define FILL_ARRAY_HIGH_LATENCY_2 for(int i = 0; i < n;i++) { arr[i] = DOUBLE_HIGH_LATENCY_2; }
#define FILL_ARRAY_LOW_LATENCY  for(int i = 0; i < n;i++) { arr[i] = DOUBLE_LOW_LATENCY; }
#define LATENCY_ADD 4.0

myInt64 cycles;
myInt64 start;
  
const unsigned int n = 10000000;
    
double *arr;
  
double t;
double t1 = 1.0,t2 = 1.0,t3 = 1.0,t4 = 1.0,t5 = 1.0,t6 = 1.0,t7 = 1.0,t8 = 1.0,t9 = 1.0,t10 = 1.0,t11 = 1.0,t12 = 1.0,t13 = 1.0,t14 = 1.0,t15 = 1.0,t16 = 1.0;
double toAdd;

void initialize_microbenchmark_data (microbenchmark_mode_t mode) {
  srand((unsigned)time(0)); 
  
  arr = (double*)malloc(n*sizeof(double));
  /* You can use to initialize some data if needed */
    switch (mode) {
        case MAX_LAT:
          t = DOUBLE_HIGH_LATENCY;
          FILL_ARRAY_HIGH_LATENCY;
          break;
        case MAX_GAP:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_LOW_LATENCY;
          break;
        case DIV_LAT:
          t = DOUBLE_HIGH_LATENCY_2;
          FILL_ARRAY_HIGH_LATENCY;
          break;
        case DIV_GAP:
          t = DOUBLE_HIGH_LATENCY;
          FILL_ARRAY_HIGH_LATENCY_2;
          break;
        case DIV_LAT_MIN:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_LOW_LATENCY;
          break;
        case DIV_GAP_MIN:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_LOW_LATENCY;
          break;
        case FOO_LAT:
          FILL_ARRAY_HIGH_LATENCY;
          toAdd = DOUBLE_HIGH_LATENCY;
          break;
        case FOO_GAP:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_HIGH_LATENCY;
          break;
        case FOO_LAT_MIN:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_LOW_LATENCY;
          toAdd = 0;
          break;
        case FOO_GAP_MIN:
          t = DOUBLE_LOW_LATENCY;
          FILL_ARRAY_LOW_LATENCY;
          break;
        default: break;
    }
}


double microbenchmark_get_max_latency() {
  start = start_tsc();

  for(int i = 0; i < n;i++) {
    t = (t > arr[i]) ? t : arr[i];
  }

  cycles = stop_tsc(start);
  return (double)cycles/n;
}

double microbenchmark_get_max_gap() {
  start = start_tsc();

  for(int i = 0; i < n;i++) {
    t1 = (1.0 > t1) ? 1.0 : t1;
    t2 = (1.0 > t2) ? 1.0 : t2;
    t3 = (1.0 > t3) ? 1.0 : t3;
    t4 = (1.0 > t4) ? 1.0 : t4;
    t5 = (1.0 > t5) ? 1.0 : t5;
    t6 = (1.0 > t6) ? 1.0 : t6;
    t7 = (1.0 > t7) ? 1.0 : t7;
    t8 = (1.0 > t8) ? 1.0 : t8;
    t9 = (1.0 > t9) ? 1.0 : t9;
    t10 = (1.0 > t10) ? 1.0 : t10;
    t11 = (1.0 > t11) ? 1.0 : t11;
    t12 = (1.0 > t12) ? 1.0 : t12;
    t13 = (1.0 > t13) ? 1.0 : t13;
    t14 = (1.0 > t14) ? 1.0 : t14;
    t15 = (1.0 > t15) ? 1.0 : t15;
    t16 = (1.0 > t16) ? 1.0 : t16;
  }
  
  cycles = stop_tsc(start);
  return ((double)cycles/n/16);
}

double microbenchmark_get_div_latency() {
  start = start_tsc();
  
  for(int i = 0; i < n;i++) {
    t = arr[i] / t;
  }
  
  cycles = stop_tsc(start);
  return (double)cycles/n;
}

double microbenchmark_get_div_gap() {
  start = start_tsc();
  
  for(int i = 0; i < n;i++) {
    arr[i] = arr[i] / t1;
    arr[i] = arr[i] / t2;
    arr[i] = arr[i] / t3;
    arr[i] = arr[i] / t4;
    arr[i] = arr[i] / t5;
    arr[i] = arr[i] / t6;
    arr[i] = arr[i] / t7;
    arr[i] = arr[i] / t8;
    arr[i] = arr[i] / t9;
    arr[i] = arr[i] / t10;
    arr[i] = arr[i] / t11;
    arr[i] = arr[i] / t12;
    arr[i] = arr[i] / t13;
    arr[i] = arr[i] / t14;
    arr[i] = arr[i] / t15;
    arr[i] = arr[i] / t16;
  }
  
  cycles = stop_tsc(start);
  return (double)cycles/n/16;
}

double microbenchmark_get_foo_latency() {
  start = start_tsc();
  
  for(int i = 1; i < n;i++) {
    arr[i] = foo(arr[i-1]) + toAdd;
  }
  
  cycles = stop_tsc(start);
  // Remove the latency of an add operation
  return (double)cycles/(n - 1) - LATENCY_ADD;
}

double microbenchmark_get_foo_gap() {
  start = start_tsc();

  for(int i = 0; i < n;i++) {
    t1 = foo(arr[i]);
    t2 = foo(arr[i]);
    t3 = foo(arr[i]);
    t4 = foo(arr[i]);
    t5 = foo(arr[i]);
    t6 = foo(arr[i]);
    t7 = foo(arr[i]);
    t8 = foo(arr[i]);
    t9 = foo(arr[i]);
    t10 = foo(arr[i]);
    t11 = foo(arr[i]);
    t12 = foo(arr[i]);
    t13 = foo(arr[i]);
    t14 = foo(arr[i]);
    t15 = foo(arr[i]);
    t16 = foo(arr[i]);
  }
  
  cycles = stop_tsc(start);
  return (double)cycles/n/16;
}
