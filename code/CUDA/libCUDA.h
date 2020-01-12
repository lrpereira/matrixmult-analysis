#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "papi.h"
#include <time.h>

void validate_rows(float *C, int SIZE);
void validate_columns(float *C, int SIZE);
void fillMatrices (float *A,float *B,float *C, int SIZE);