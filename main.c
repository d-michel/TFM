#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>

/*int basisDimension(int N, int gamma) {
	int dim;
	dim = (int) ((gamma + 1)*(N + 1)*(N + 2)/2);
	return dim;
}*/

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

void initVec(tvec *vec) {
	for (int i = 0; i < 3; i++) {
		vec->n[i] = 0;
	}
	vec->nph = 0;
}

tvec copy(tvec vecin) {
	tvec vecout;
	for (int i = 0; i < 3; i++) {
		vecout.n[i] = vecin.n[i];
	}
	vecout.nph = vecin.nph;
	return vecout;
}

void initBasis(tbasis *basis) {
	basis->nvec = 0;
}

void push(tbasis *basis, tvec vec) {
	if (basis->nvec == dim) {
		printf("\n Full tbasis \n");
	} else {
		basis->vectors[basis->nvec] = copy(vec);
		basis->nvec++;
	}
}

tbasis basisVectors() {
		
	tbasis basis;
	initBasis(&basis);
	
	tvec vec;
	for (int i = 0; i <= N; i++) {
		int N0, N1, N2;
		N0 = N - i;
		N1 = 0;
		N2 = i;
		while (N0 >= 0) {
			initVec(&vec);
			vec.n[0] = N0;
			vec.n[1] = N1;
			vec.n[2] = N2;
			for (int j = 0; j <= gamma; j++) {
				vec.nph = j;
				push(&basis, vec);
			}
			N0--;
			N1++;
		}
	}
	//printf("\n %d \n", basis.nvec);
	return basis;
}

float Hact(tvec vi, tvec vj, float wf, float wn[3], float D, float f[3], float E) {
	float H;
	H = 0;
	
	if (vi.n[0] == vj.n[0] && vi.n[1] == vj.n[1] && vi.n[2] == vj.n[2]) {
		if (vi.nph == vj.nph) {
			
			H += wf*vi.nph + D*(2*vi.nph + 1);
				
			for (int alpha = 0; alpha < 3; alpha++) {
				H += wn[alpha]*vi.n[alpha];
			}
			
		} else if ((vi.nph + 1) == vj.nph) {
			
			H += E*sqrt((vi.nph + 1)/N);
			
		} else if ((vi.nph - 1) == vj.nph) {
			
			H += E*sqrt(vi.nph/N);
			
		} else if ((vi.nph + 2) == vj.nph) {
			
			H += D*sqrt((vi.nph + 1)*(vi.nph + 2));
			
		} else if ((vi.nph - 2) == vj.nph) {
			
			H += D*sqrt(vi.nph*(vi.nph - 1));
			
		}
	}
	
	int c[3] = {2, 1, 0};
	float Om[3][3];
	for (int b = 1; b <= 2; b++) {
		for (int a = 0; a < b; a++) {
			Om[a][b] = sqrt((wn[b] - wn[a])*f[a+b-1]*D);
			if (vi.n[c[a+b-1]] == vj.n[c[a+b-1]]) {
				if ((vi.nph - 1) == vj.nph) {
					if (((vi.n[a] + 1) == vj.n[a]) && ((vi.n[b] - 1) == vj.n[b])) {
						H += Om[a][b]*sqrt((vi.n[a] + 1)*vi.n[b]*vi.nph/N);
					} else if (((vi.n[a] - 1) == vj.n[a]) && ((vi.n[b] + 1) == vj.n[b])) {
						H += Om[a][b]*sqrt((vi.n[b] + 1)*vi.n[a]*vi.nph/N);
					}
				} else if ((vi.nph + 1) == vj.nph) {
					if (((vi.n[a] + 1) == vj.n[a]) && ((vi.n[b] - 1) == vj.n[b])) {
						H += Om[a][b]*sqrt((vi.n[a] + 1)*vi.n[b]*(vi.nph + 1)/N);
					} else if (((vi.n[a] - 1) == vj.n[a]) && ((vi.n[b] + 1) == vj.n[b])) {
						H += Om[a][b]*sqrt((vi.n[b] + 1)*vi.n[a]*(vi.nph + 1)/N);
					}
				}
			}
		}
	}
	return H;
}

float H(tbasis basis, float wf, float wn[3], float D, float f[3], float E) {
	float Hmatrix[dim][dim];
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			Hmatrix[i][j] = Hact(basis.vectors[i], basis.vectors[j], wf, wn, D, f, E);
		}
	}
	return Hmatrix;
}

float eigenvalues(HMatrix) {
	float eigvals [dim];
	info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'N', 'U', dim, &HMatrix, dim, &eigvals);
	if (info > 0) {
		printf("the algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
	return eigvals;
}

int main(int argc, char **argv)
{
	tbasis basis;
	basis = basisVectors();
	for (int i = 0; i < basis.nvec; i++) {
		printf("%d, %d, %d, %d\n", basis.vectors[i].n[0], basis.vectors[i].n[1], basis.vectors[i].n[2], basis.vectors[i].nph);
	}
	return 0;
}
