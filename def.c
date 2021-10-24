#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "lapacke.h"
#include "include/def.h"

/*int basisDimension(int N, int gamma) {
	int dim;
	dim = (int) ((gamma + 1)*(N + 1)*(N + 2)/2);
	return dim;
}*/

tvec * initVec() {
	tvec * vec = malloc(sizeof(tvec));
	for (int i = 0; i < 3; i++) {
		vec->n[i] = 0;
	}
	vec->nph = 0;
	return vec;
}

tvec copy(tvec *vecin) {
	tvec vecout;
	for (int i = 0; i < 3; i++) {
		vecout.n[i] = vecin->n[i];
	}
	vecout.nph = vecin->nph;
	return vecout;
}

tbasis * initBasis() {
	tbasis * basis = malloc(sizeof(tbasis));
	basis->nvec = 0;
	return basis;
}

void push(tbasis *basis, tvec *vec) {
	if (basis->nvec == dim) {
		printf("\n Full tbasis \n");
	} else {
		basis->vectors[basis->nvec] = copy(vec);
		basis->nvec++;
	}
}

tmatrix * initMatrix(int size_row, int size_col) {
	tmatrix * matrix = malloc(sizeof(tmatrix *));
	
	matrix->nrow = size_row;
	matrix->ncol = size_col;
	matrix->data = malloc(matrix->nrow*sizeof(float *));
	
	for (int i = 0; i < matrix->ncol; i++) {
		matrix->data[i] = malloc(sizeof(float));
	}
	
	return matrix;
}

tbasis * basisVectors() {
		
	tbasis * basis = initBasis();
	tvec * vec = initVec();
	
	for (int i = 0; i <= N; i++) {
		int N0, N1, N2;
		N0 = N - i;
		N1 = 0;
		N2 = i;
		while (N0 >= 0) {
			initVec(&vec);
			vec->n[0] = N0;
			vec->n[1] = N1;
			vec->n[2] = N2;
			for (int j = 0; j <= gamma; j++) {
				vec->nph = j;
				push(basis, vec);
			}
			N0--;
			N1++;
		}
	}
	//printf("\n %d \n", basis.nvec);
	return basis;
}

float Hact(tvec *vi, tvec *vj, float wf, float wn[3], float D, float f[3], float E) {
	float H;
	H = 0;
	
	if (vi->n[0] == vj->n[0] && vi->n[1] == vj->n[1] && vi->n[2] == vj->n[2]) {
		if (vi->nph == vj->nph) {
			
			H += wf*vi->nph + D*(2*vi->nph + 1);
				
			for (int alpha = 0; alpha < 3; alpha++) {
				H += wn[alpha]*vi->n[alpha];
			}
			
		} else if ((vi->nph + 1) == vj->nph) {
			
			H += E*sqrt((vi->nph + 1)/N);
			
		} else if ((vi->nph - 1) == vj->nph) {
			
			H += E*sqrt(vi->nph/N);
			
		} else if ((vi->nph + 2) == vj->nph) {
			
			H += D*sqrt((vi->nph + 1)*(vi->nph + 2));
			
		} else if ((vi->nph - 2) == vj->nph) {
			
			H += D*sqrt(vi->nph*(vi->nph - 1));
			
		}
	}
	
	int c[3] = {2, 1, 0};
	float Om[3][3];
	for (int b = 1; b <= 2; b++) {
		for (int a = 0; a < b; a++) {
			Om[a][b] = sqrt((wn[b] - wn[a])*f[a+b-1]*D);
			if (vi->n[c[a+b-1]] == vj->n[c[a+b-1]]) {
				if ((vi->nph - 1) == vj->nph) {
					if (((vi->n[a] + 1) == vj->n[a]) && ((vi->n[b] - 1) == vj->n[b])) {
						H += Om[a][b]*sqrt((vi->n[a] + 1)*vi->n[b]*vi->nph/N);
					} else if (((vi->n[a] - 1) == vj->n[a]) && ((vi->n[b] + 1) == vj->n[b])) {
						H += Om[a][b]*sqrt((vi->n[b] + 1)*vi->n[a]*vi->nph/N);
					}
				} else if ((vi->nph + 1) == vj->nph) {
					if (((vi->n[a] + 1) == vj->n[a]) && ((vi->n[b] - 1) == vj->n[b])) {
						H += Om[a][b]*sqrt((vi->n[a] + 1)*vi->n[b]*(vi->nph + 1)/N);
					} else if (((vi->n[a] - 1) == vj->n[a]) && ((vi->n[b] + 1) == vj->n[b])) {
						H += Om[a][b]*sqrt((vi->n[b] + 1)*vi->n[a]*(vi->nph + 1)/N);
					}
				}
			}
		}
	}
	return H;
}

tmatrix * H(tbasis *basis, float wf, float wn[3], float D, float f[3], float E) {
	tmatrix * Hmatrix = initMatrix(dim,dim);
	for (int i = 0; i < Hmatrix->nrow; i++) {
		for (int j = 0; j < Hmatrix->ncol; j++) {
			Hmatrix->data[i][j] = Hact(&basis->vectors[i], &basis->vectors[j], wf, wn, D, f, E);
		}
	}
	return Hmatrix;
}

/*float * eigenvalues(float *HMatrix) {
	float * eigvals = malloc(dim*sizeof(float));
	int info = 1;
	// LAPACKE_ssyev(int matrix_layout, char jobz, char uplo, lapack_int n, float* a, lapack_int lda, float* w);
//	info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'N', 'U', (lapack_int) dim, HMatrix, (lapack_int) dim, eigvals);
	if (info > 0) {
		printf("the algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
	return eigvals;
}*/