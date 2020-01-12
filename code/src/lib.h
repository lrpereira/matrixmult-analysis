#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "papi.h"

void init_app(int argc);
float** init_matrix(int size);
void fill_matrix(float** m, int size, int content);
void print_matrix(float** m, int size);
void transpose_matrix(float** m, int size);
int val_result_axb(float** c, int size);
int val_result_bxa(float** c, int size);
long long calc_miss_rate(long long tcm, long long tca);
void free_matrices(float** a, float** b, float** c1, int size);
