#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>

#define N 10
#define gamma 1
#define dim (int) ((gamma + 1)*(N + 1)*(N + 2)/2)

typedef struct {
	int n[3];
	int nph;
} tvec;

typedef struct {
	tvec vectors[dim];
	int nvec;
} tbasis;

typedef struct {
	float ** data;
	int nrow;
	int ncol;
} tmatrix;

tvec * initVec();

tvec copy(tvec *vecin);

tbasis * initBasis();

void push(tbasis *basis, tvec *vec);

tmatrix * initMatrix(int size_row, int size_col);

tbasis * basisVectors();

float Hact(tvec *vi, tvec *vj, float wf, float wn[3], float D, float f[3], float E);

tmatrix * H(tbasis *basis, float wf, float wn[3], float D, float f[3], float E);

float * eigenvalues(float *HMatrix);
