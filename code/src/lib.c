#include "lib.h"

void init_app(int argc)
{
    if (argc!=3) {
        printf("Usage: ./binz size_sqMatrix kernel_choice\n");
        exit(0);
    }
}

float** init_matrix(int size)
{
    float** aux = calloc(size, sizeof(float*));
    for (int i=0; i<size; i++) {
        aux[i] = calloc(size, sizeof(float));
    }
    return aux;
}

void fill_matrix(float** m, int size, int content)
{
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if (content == 1) {
                m[i][j] = 1;
            }
            else if (content == 9) {
                m[i][j] = ((float) rand()) / ((float) RAND_MAX);
            }
        }
    }
}

void print_matrix(float** m, int size)
{
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            printf("%f ", m[i][j]);
        }
        printf("\n");
    }
}

void transpose_matrix(float** m, int size)
{
    float tmp;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < i; ++j) {
            tmp = m[i][j];
            m[i][j] = m[j][i];
            m[j][i] = tmp;
        }
    }
}

int val_result_axb(float** c, int size)
{
    for (int i=0; i<size; i++)
        for (int j=0; j<size-1; j++)
            if (c[i][j]!=c[i][j+1])
                return -1;
    return 0;
}

int val_result_bxa(float** c, int size)
{
    for (int j=0; j<size; j++)
        for (int i=0; i<size-1; i++)
            if (c[i][j]!=c[i+1][j])
                return -1;
    return 0;
}

long long calc_miss_rate(long long tcm, long long tca)
{
    return (tcm/tca);
}

void free_matrices(float** a, float** b, float** c1, int size)
{
    for (int i=0; i<size; i++) {
        free(a[i]); free(b[i]); free(c1[i]);
    }
    free(a); free(b); free(c1);
}
