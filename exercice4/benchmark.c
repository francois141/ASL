//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
//#include <windows.h> // Include if under windows
#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define alpha 0.2

/*
 *	Initialize the input
 */
static inline double generate_double() {
    return ((double)rand()) / INT_MAX;
}

void fill_vector(double *x, int n) {
    for(int i=0; i < n; i++) {
        x[i] = generate_double();
    }
}

/*
 * Sort algorithm to sort the results
*/
void quick_sort(double *numbers, int left, int right)
{
    int l_hold = left;
    int r_hold = right;
    int pivot = numbers[left];

    while (left < right) {
        while ((numbers[right] >= pivot) && (left < right))
            right--;
        if (left != right) {
            numbers[left] = numbers[right];
            left++;
        }
        while ((numbers[left] <= pivot) && (left < right))
            left++;
        if (left != right) {
            numbers[right] = numbers[left];
            right--;
        }
    }

    numbers[left] = pivot;
    pivot = left;
    left = l_hold;
    right = r_hold;

    if (left < pivot)
        quick_sort(numbers, left, pivot-1);
    if (right > pivot)
        quick_sort(numbers, pivot+1, right);
}

/* 
 * execute operation
 */
void compute(double *x, double *y, int n) {
    double s = 0.0;
    for(int i = 0; i < n;i++) {
        s = (s + x[i]*x[i]) + y[i]*y[i]*y[i];
    }
    x[0] = s;
}

void compute_optimized(double *x, double *y, int n) {
    double acc1 = 0.0;
    double acc2 = 0.0;
    for(int i = 0; i < n;i++) {
       acc1 = acc1 + x[i]*x[i];
       acc2 = acc2 + y[i]*y[i]*y[i];
    }
    x[0] = acc1 + acc2;
}

/* 
 * Timing functions based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(double A[], double B[], int n) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /* 
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to 
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            compute(A, B, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(A, B, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

#ifdef __x86_64__
double rdtsc_optimized(double A[], double B[], int n) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /* 
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to 
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            compute_optimized(A, B, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute_optimized(A, B, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

int main(int argc, char **argv) {   
#ifdef __x86_64__
    srand(time(0));

    for(unsigned int nb_iterations = (1 << 4); nb_iterations <= (1 << 23); nb_iterations <<= 1) {
        const unsigned int n = nb_iterations;
        const unsigned int nb_tests = 30;

        double* A = (double *)malloc(n*sizeof(double));
        double* B = (double *)malloc(n*sizeof(double));

        double values[nb_tests];
        double values_optimized[nb_tests];

        for(int i = 0; i < nb_tests; i++) {
            fill_vector(A, n);
            fill_vector(B, n);
            double r = rdtsc(A, B, n);

            fill_vector(A, n);
            fill_vector(B, n);
            double r_optimized = rdtsc_optimized(A, B, n);

            values[i] = r;
            values_optimized[i] = r_optimized;
        }

        free(A);
        free(B);

        quick_sort(values,0,nb_tests);
        quick_sort(values_optimized,0,nb_tests);

        // N is always even ==> Median is average of arr[idx1] and arr[idx2]
        int idx1 = nb_tests/2 - 1;
        int idx2 = nb_tests/2;

        double cycles = (values[idx1] + values[idx2]) / 2.0;
        double cycles_optimized = (values_optimized[idx1] + values_optimized[idx2]) / 2.0;

        printf("=========================================\n");
        printf("N = %d \n",n);
        printf("=========================================\n");

        printf("[BASE VERSION]\n");
        printf("RDTSC instruction:\n %lf cycles measured => %lf seconds, assuming frequency is %lf MHz. (change in source file if different)\n\n", cycles, cycles/(FREQUENCY), (FREQUENCY)/1e6);

        printf("[OPTIMIZED VERSION]\n");
        printf("RDTSC instruction:\n %lf cycles measured => %lf seconds, assuming frequency is %lf MHz. (change in source file if different)\n\n", cycles_optimized, cycles_optimized/(FREQUENCY), (FREQUENCY)/1e6);
    }
#else
    printf("This benchmark works only for a x86_64 architecture");
#endif
    return 0;
}

