#include "libCUDA.h"


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