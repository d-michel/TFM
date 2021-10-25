#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>
#include "include/def.h"

int main(int argc, char **argv) {
	int N = 12;
	int gamma = 1;
	int dim = basisDimension(N, gamma);
	float wf = 1;
	float wn[3] = {1, 1.1, 2.1};
	float D = 3;
	float f[3] = {6, 6, 6};
	float E = 0;
	
	tbasis * basis = initBasis(dim);
	basis = basisVectors(N, gamma, dim);
	/* for (int i = 0; i < basis->nvec; i++) {
		printf("%d, %d, %d, %d\n", basis->vectors[i].n[0], basis->vectors[i].n[1], basis->vectors[i].n[2], basis->vectors[i].nph);
	} */
	
	tmatrix * Hmatrix = initMatrix(dim, dim);
	Hmatrix = H(basis, N, dim, wf, wn, D, f, E);
	//Hmatrix->data[181][181] = Hact(&basis->vectors[181], &basis->vectors[181], wf, wn, D, f, E);
	//printf("%d, %d\n", Hmatrix->nrow, Hmatrix->ncol);
	//for (int i = 0; i < Hmatrix->nrow; i++) {
	//	for (int j = 0; j < Hmatrix->ncol; j++) {
	//		printf("%f", Hmatrix->data[i][j]);
	//	}
	//	printf("\n");
	//}
	
	return 0;
}
