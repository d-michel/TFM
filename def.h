#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>

typedef struct {
	int n[3];
	int nph;
} tvec;

typedef struct {
	tvec * data;
	int nvec;
} tbasis;

typedef struct {
	float ** data;
	int nrow;
	int ncol;
} tmatrix;

tvec * initVec();

tvec copy(tvec *vecin);

tbasis * initBasis(int size_basis);

void push(tbasis *basis, tvec *vec, int max_size_basis);

tmatrix * initMatrix(int size_row, int size_col);

tbasis * basisVectors(int N, int gamma, int dim);

float Hact(tvec *vi, tvec *vj, int N, float wf, float wn[3], float D, float f[3], float E);

tmatrix * H(tbasis *basis, int N, int dim, float wf, float wn[3], float D, float f[3], float E);

float * eigenvalues(float *HMatrix);
