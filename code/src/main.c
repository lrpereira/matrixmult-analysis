#include "lib.h"

double clearcache [30000000];
#define NUM_EVENTS 3

void clearCache (void)
{
    for(int i = 0; i < 30000000; ++i)
        clearcache[i] = i;
}

void compute (float* x, float* y, float* r)
{
    (*r) += (*x) * (*y);
}

void multiplication_ijk(float** a, float** b, float** c, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                compute(&a[i][k], &b[k][j], &c[i][j]);
}

void multiplication_ijk_transpose(float** a, float** b, float** c, int size)
{
    transpose_matrix(b, size);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            for (int k = 0; k < size; k++)
                compute(&a[i][k], &b[j][k], &c[i][j]);
}

void multiplication_ikj(float** a, float** b, float** c, int size)
{
    for (int i = 0; i < size; i++)
        for (int k = 0; k < size; k++)
            for (int j = 0; j < size; j++)
                compute(&a[i][k], &b[k][j], &c[i][j]);
}

void multiplication_jki(float** a, float** b, float** c, int size)
{
    for (int j = 0; j < size; j++)
        for (int k = 0; k < size; k++)
            for (int i = 0; i < size; i++)
                compute(&a[i][k], &b[k][j], &c[i][j]);
}

void multiplication_jki_transpose(float** a, float** b, float** c, int size)
{
    transpose_matrix(a, size);
    transpose_matrix(b, size);

    for (int j = 0; j < size; ++j)
        for (int k = 0; k < size; ++k)
            for (int i = 0; i < size; ++i)
                compute(&a[k][i], &b[j][k], &c[j][i]);

    transpose_matrix(c, size);
}

int main(int argc, char *argv[])
{
    //srand(time(0));

    init_app(argc);
    int size = atoi(argv[1]);
    int option = atoi(argv[2]);

    double start, end;
    int Events[NUM_EVENTS] = {PAPI_L1_TCM, PAPI_LD_INS, PAPI_SR_INS};
                              //PAPI_L2_TCM, PAPI_L2_TCA, PAPI_L3_TCM, PAPI_L3_TCA};
    int EventSet = PAPI_NULL;
    long long papi[NUM_EVENTS];
    int retval = 0;

    float** a = init_matrix(size);
    float** b = init_matrix(size);
    float** c1 = init_matrix(size);
    //float** c2 = init_matrix(size);

    fill_matrix(a, size, 9);
    fill_matrix(b, size, 1);

    /* Initialize the PAPI Library */
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval > 0) {
        fprintf(stderr,"PAPI library version mismatch!\n");
        exit(1);
    }

    /* Allocate space for the new eventset and do setup */
    if((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK) {
        fprintf(stderr, "PAPI create event set error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    /* Add Flops and total cycles to the eventset */
    if((retval = PAPI_add_events(EventSet, Events, NUM_EVENTS)) != PAPI_OK) {
        fprintf(stderr, "PAPI add event set error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    //PAPI START
    if((retval = PAPI_start(EventSet)) != PAPI_OK) {
        fprintf(stderr, "PAPI start error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    start = omp_get_wtime();
    switch(option){
    case 1:
        multiplication_ijk(a, b, c1, size);
        break;
    case 2:
        multiplication_ikj(a, b, c1, size);
        break;
    case 3:
        multiplication_jki(a, b, c1, size);
        break;
    case 4:
        multiplication_ijk_transpose(a, b, c1, size);
        break;
    case 5:
        multiplication_jki_transpose(a, b, c1, size);
        break;
    default:
        printf("Option invalid.\n");
        break;
    }
    end = omp_get_wtime();

    //PAPI STOP
    if((retval = PAPI_stop(EventSet,papi)) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    double timeMs = (end-start)*1000;
    FILE* fp = fopen("resultados.csv", "a");

    float miss_rate_l1 = calc_miss_rate(papi[0], papi[1]+papi[2]);
    // float miss_rate_l2 = calc_miss_rate(papi[3], papi[4]);
    // float miss_rate_l3 = calc_miss_rate(papi[5], papi[6]);

    fprintf(fp, "%d\n%lf\n%f\n", option, timeMs, miss_rate_l1);

    fclose(fp);

    return 0;
}
