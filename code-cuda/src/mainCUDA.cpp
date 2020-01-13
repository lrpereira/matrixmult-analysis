#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "papi.h"
#include <sys/time.h>

double clearcache[30000000];
#define NUM_EVENTS 3
#define TIME_RESOLUTION 1000000
#define BLOCK_SIZE 32
long long unsigned initial_time;
long long unsigned execTime;
timeval t;


void validate_rows(float *C, int SIZE){
	int i, j, r = 1;

	for (i = 0; i < SIZE && r; i++){
		for (j = 0; j < SIZE && r; j++){
			if (C[i * SIZE + j] != C[i * SIZE + 0]){
				r = 0;
			}
		}
	}

	if(!r) {
		printf("ERRRO");
	}
}

void validate_columns(float *C, int SIZE){
	int i, j, r = 1;

	for (i = 0; i < SIZE && r; i++){
		for (j = 0; j < SIZE && r; j++){
			if (C[i * SIZE + j] != C[0 * SIZE + j]){
				r = 0;
			}
		}
	}

	if(!r) {
		printf("ERRRO");
	}
}

void fillMatrices (float *A,float *B,float *C, int SIZE) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      A[i * SIZE + j] = ((float) rand()) / ((float) RAND_MAX);
      B[i * SIZE + j] = 1;
      C[i * SIZE + j] = 0;
    }
  }
}

void start(){
	gettimeofday(&t, NULL);

	//Convert current time to milliseconds
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

void stop(){
	gettimeofday(&t, NULL);

	//Convert current time to milliseconds
	long long unsigned currentTime;
	currentTime = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
	execTime = currentTime - initial_time;
}

void clearCache(){
	for(int i = 0; i < 30000000; i++){
		clearcache[i] = i;
	}
}

__global__ void testKernel(float *a, float *b, float *c, int size){
}

void useCUDA(float *a, float *b, float *c, int size){
	float *aA, *bB, *cC;
	cudaMalloc((void**) &aA, size*size*sizeof(float));
	cudaMalloc((void**) &bB, size*size*sizeof(float));
	cudaMalloc((void**) &cC, size*size*sizeof(float));

	cudaMemcpy(aA, a, size*size*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(bB, b, size*size*sizeof(float), cudaMemcpyHostToDevice);

	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(size/BLOCK_SIZE + (size % BLOCK_SIZE > 0), size/BLOCK_SIZE + (size % BLOCK_SIZE > 0));
	cudaKernel<<<dimGrid, dimBlock>>>(aA, bB, cC, size);

	//Bloqueia ate tudo estar concluido
	cudaDeviceSynchronize();
	cudaMemcpy(cC, c, size*size*sizeof(float), cudaMemcpyHostToDevice);
	cudaFree(aA);
	cudaFree(bB);
	cudaFree(cC);
}

__global__ void cudaKernel(float *a, float *b, float *c, int size){
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;
	int k;

  	if (row>=size || col>=size){
    float res=0.0f;
    for(k=0; k<size; ++k){
      res += a[row*size + k] * b[k*size + col];
    }
    c[row*size + col]=res;
  }
}

void dotProductCUDA(float *a, float *b, float *c, int size){
  	float *aA, *bB, *cC;
  	cudaMalloc( (void**) &aA, size * size *sizeof(float) );
  	cudaMalloc( (void**) &bB, size * size *sizeof(float) );
  	cudaMalloc( (void**) &cC, size * size *sizeof(float) );

  	cudaMemcpy(aA, a, size * size * sizeof(float), cudaMemcpyHostToDevice);
  	cudaMemcpy(bB, b, size * size * sizeof(float), cudaMemcpyHostToDevice);

  	dim3 dimBlock (BLOCK_SIZE, BLOCK_SIZE);
  	dim3 dimGrid (size / BLOCK_SIZE + (size % BLOCK_SIZE>0), size / BLOCK_SIZE + (size % BLOCK_SIZE>0) );
  	cudaKernel<<<dimGrid, dimBlock>>>(aA, bB, cC, size);

  	cudaDeviceSynchronize();
  	cudaMemcpy(c, cC, size * size * sizeof(float), cudaMemcpyDeviceToHost);
  	cudaFree(aA);
  	cudaFree(bB);
  	cudaFree(cC);
}


__global__ void cudaBlockKernel(float *a, float *b, float *c, int size){
	float value = 0;

	int row = blockId.y * BLOCK_SIZE + threadIdx.y;
	int col = blockId.x * BLOCK_SIZE + threadIdx.x;

	__shared__ float As[BLOCK_SIZE][BLOCK_SIZE];
  	__shared__ float Bs[BLOCK_SIZE][BLOCK_SIZE];

  	for (int k = 0; k < (BLOCK_SIZE + size - 1)/BLOCK_SIZE; k++) {
    	if (k*BLOCK_SIZE + threadIdx.x < size && row < size){
      		As[threadIdx.y][threadIdx.x] = a[row*size + k*BLOCK_SIZE + threadIdx.x];
    	}else{
      		As[threadIdx.y][threadIdx.x] = 0.0;
    	}

    	if (k*BLOCK_SIZE + threadIdx.y < size && col < size){
      		Bs[threadIdx.y][threadIdx.x] = b[(k*BLOCK_SIZE + threadIdx.y)*size + col];
    	}else{
      		Bs[threadIdx.y][threadIdx.x] = 0.0;
    	}

    	__syncthreads();

    	for (int n = 0; n < BLOCK_SIZE; ++n){
      		value += As[threadIdx.y][n] * Bs[n][threadIdx.x];
    	}

    	__syncthreads();
  	}

	if (row < size && col < size){
    	c[((blockIdx.y * blockDim.y + threadIdx.y)*size) + (blockIdx.x * blockDim.x)+ threadIdx.x] = value;
  	}
}

void dotProductBlockCUDA(float *a, float *b, float *c, int size){
  float *aA, *bB, *cC;

  size_t size = size * size * sizeof(float);

  cudaMalloc(&aA, size);
  cudaMalloc(&bB, size);

  cudaMemcpy(aA, a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(bB, b, size, cudaMemcpyHostToDevice);

  cudaMalloc(&cC, size);

  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid((size + dimBlock.x -1) / dimBlock.x, (size + dimBlock.y -1 )/ dimBlock.y);

  cudaBlockKernel<<<dimGrid, dimBlock>>>(aA, bB, cC, size);

  cudaDeviceSynchronize();

  cudaMemcpy(c, cC, size, cudaMemcpyDeviceToHost);

  cudaFree(aA);
  cudaFree(bB);
  cudaFree(cC);
}

void run_dotproduct(float *A, float *B, float *C, int SIZE, void (*f) (float*,float*,float*, int)){

	//fill matrix with random numbers
	fillMatrices(A, B, C, SIZE);

  	start();
  	f(A,  B,  C, SIZE);
  	stop();

  	validate_rows(C, SIZE);

  	printf("Execution time taked: %llu milliseconds", execTime);
}

int main (int argc, char *argv[]) {
    
    float *A, *B, *C;
    int size = atoi(argv[2]);

    A = (float *) malloc(size * size * sizeof(float));
    B = (float *) malloc(size * size * sizeof(float));
    C = (float *) malloc(size * size * sizeof(float));

    int Events[NUM_EVENTS] = {PAPI_L1_TCM, PAPI_LD_INS, PAPI_SR_INS};
                              //PAPI_L2_TCM, PAPI_L2_TCA, PAPI_L3_TCM, PAPI_L3_TCA};
    int EventSet = PAPI_NULL;
    long long papi[NUM_EVENTS];
    int retval = 0;


    void (*implentation) (float *,float *,float *,int)


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


    //define function
    if(!strcmp("dotProductCUDA", argv[1])){
    	implentation = dotProductCUDA;
    }else if(!strcmp("dotProductBlockCUDA", argv[1])){
    	implentation = dotProductBlockCUDA;
    }else if(!strcmp("useCUDA", argv[1])){
    	implentation = useCUDA;
    }

    //execute functions
    if(!strcmp("time", argv[3])){
    	run_dotproduct(a, b, c1, size, implentation);
    }


    //PAPI STOP
    if((retval = PAPI_stop(EventSet,papi)) != PAPI_OK) {
        fprintf(stderr, "PAPI stop error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    FILE* fp = fopen("resultadosCUDA.csv", "a");

    float miss_rate_l1 = calc_miss_rate(papi[0], papi[1]+papi[2]);
    // float miss_rate_l2 = calc_miss_rate(papi[3], papi[4]);
    // float miss_rate_l3 = calc_miss_rate(papi[5], papi[6]);

    fprintf(fp, "%llu\n%f\n", timeMs, miss_rate_l1);

    fclose(fp);

    return 0;
}