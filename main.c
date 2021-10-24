#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>
#include "include/def.h"

int main(int argc, char **argv) {
	
	float wf = 1;
	float wn[3] = {1, 1.1, 2.1};
	float D = 3;
	float f[3] = {6, 6, 6};
	float E = 0;
	
	tbasis * basis = initBasis();
	basis = basisVectors();
	/* for (int i = 0; i < basis->nvec; i++) {
		printf("%d, %d, %d, %d\n", basis->vectors[i].n[0], basis->vectors[i].n[1], basis->vectors[i].n[2], basis->vectors[i].nph);
	} */
	
	tmatrix * Hmatrix = initMatrix(dim, dim);
	Hmatrix = H(basis, wf, wn, D, f, E);

	/*for (int i = 0; i < Hmatrix->nrow; i++) {
		for (int j = 0; j < Hmatrix->ncol; j++) {
	printf("%.2f, ", Hmatrix->data[0][0]);
		}
		printf("\n");
	}*/
	
	return 0;
}
