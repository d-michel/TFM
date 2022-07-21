#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>
#include <lapacke_utils.h>
//#include <cblas.h>
//#include <omp.h>
#include "include/def.h"

//// CLASS OPERATIONS ////

/* STATE
typedef struct {
	float * n;		// (3-array) number of atoms in each energy level
	float nph;		// number of photons
} tstate;
*/
// INITIALIZE STATE
tstate * init_state() {
	tstate * state = malloc(sizeof(tstate));	// (state-malloc) create malloc of tstate
	state->n = malloc(3*sizeof(float));			// (3*float-malloc) create malloc of 3-array (number of atoms)
	
	return state;
}

/* MIXED STATE
typedef struct {
	tstate ** data;
	float * width;
	int nstates;
} tmixstate;
*/

// INITIALIZE MIXED STATE
tmixstate * init_mix_state(int n_states) {
	tmixstate * mixstate = malloc(sizeof(tmixstate));
	mixstate->nstates = n_states;
	mixstate->width = malloc(mixstate->nstates*sizeof(float));
	mixstate->data = malloc(mixstate->nstates*sizeof(tstate*));
	for (int i = 0; i < mixstate->nstates; i++) {
		mixstate->data[i] = init_state();
	}
	
	return mixstate;
}

/* BASIS
typedef struct {
	tstate ** data;		// (2D-array) basis states
	int nstates;		// number of states
} tbasis;
*/
// INITIALIZE BASIS
tbasis * init_basis(int size_basis) {
	tbasis * basis = malloc(sizeof(tbasis));					// (tbasis-malloc) create malloc of tbasis
	
	basis->nstates = size_basis;								// number of states equal to basis size
	basis->data = malloc(basis->nstates*sizeof(tstate *));		// (nstates*tstate-malloc) create malloc of states data
	
	for (int i = 0; i < basis->nstates; i++) {
		basis->data[i] = init_state();							// initialize each state of data
	}
	
	return basis;
}

/* FLOAT VECTOR
typedef struct {
	float * data;		// (1D-array) vector data
	int length;			// vector length
} tvector;
*/

// INITIALIZE FLOAT VECTOR
tvector * init_vector(int size_vector) {
	tvector * vector = malloc(sizeof(tvector));				// (tvector-malloc) create malloc of tvector
	
	vector->length = size_vector;							// vector length equal to vector size
	vector->data = malloc(vector->length*sizeof(float));	// (length*float-malloc) create malloc of vector data

	return vector;
}

float mean_vector(tvector *vector, int start, int finish) {
	float mean = 0;
	
	for (int i = start; i < finish; i++) {
		mean += vector->data[i];
	}
	mean /= (float) (finish - start);
	
	return mean;
}

/* INTEGER VECTOR
typedef struct {
	int * data;			// (1D-array) vector data
	int length;			// vector length
} tvector_int;
*/

// INITIALIZE INTEGER VECTOR
tvector_int * init_vector_int(int size_vector) {
	tvector_int * vector = malloc(sizeof(tvector_int));				// (tvector_int-malloc) create malloc of tvector_int
	
	vector->length = size_vector;							// vector length equal to vector size
	vector->data = malloc(vector->length*sizeof(int));	// (length*int-malloc) create malloc of vector data

	return vector;
}

/* SQUARE MATRIX
typedef struct {
	float ** data;		// (2D-array) matrix data
	int nrow;			// number of rows
	int ncol;			// number of columns
} tsmatrix;
*/

// INITIALIZE SQUARE MATRIX
tsmatrix * init_square_matrix(int size_row, int size_col) {
	tsmatrix * smatrix = malloc(sizeof(tsmatrix));					// (tsmatrix-malloc) create malloc of s-matrix
	
	smatrix->nrow = size_row;										// number of rows equal to row size
	smatrix->ncol = size_col;										// number of columns equal to column size
	smatrix->data = malloc(smatrix->nrow*sizeof(float *));			// (nrow*float_pointer-malloc) cr mall of row data
	
	for (int i = 0; i < smatrix->nrow; i++) {						// for each row
		smatrix->data[i] = malloc(smatrix->ncol*sizeof(float));		// (ncol*float-malloc) create malloc of columns data
	}
	
	return smatrix;
}

// COPY COLUMN
void copy_column(tsmatrix *out_matrix, tsmatrix *in_matrix, int id_out, int id_in) {
	for (int i = 0; i < in_matrix->nrow; i++) {
		out_matrix->data[i][id_out] = in_matrix->data[i][id_in];
	}
}

// COPY ROW
void copy_row(tsmatrix *out_matrix, tsmatrix *in_matrix, int id_out, int id_in) {
	for (int i = 0; i < in_matrix->ncol; i++) {
		out_matrix->data[id_out][i] = in_matrix->data[id_in][i];
	}
}

// MATRIX - MATRIX PRODUCT
tsmatrix * matrix_matrix_product(tsmatrix *matrix1, tsmatrix *matrix2) {
	if (matrix1->ncol != matrix2->nrow) {
		printf("Matrices size must agree (matrix_matrix_product error).\n");
		exit(1);
	}
	
	tsmatrix * matrix_result = init_square_matrix(matrix1->nrow, matrix2->ncol);
	
	//#pragma omp parallel private(i,j,k) shared(matrix_result,matrix1,matrix2)
	//{
	//#pragma omp for
	for (int i = 0; i < matrix_result->nrow; i++) {
		for (int k = 0; k < matrix1->ncol; k++) {
		//for (int j = 0; j < matrix_result->ncol; j++) {
			//matrix_result->data[i][j] = 0;
			for (int j = 0; j < matrix_result->ncol; j++) {
			//for (int k = 0; k < matrix1->ncol; k++) {
				matrix_result->data[i][j] += matrix1->data[i][k]*matrix2->data[k][j];
			}
		}
	}
	//}
	//sgemm(transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
	//cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, matrix_result->nrow, matrix_result->ncol, matrix1->ncol, 1, matrix1->data, matrix1->nrow, matrix2->data, matrix2->nrow, 0, matrix_result->data, matrix_result->nrow);
	
	return matrix_result;
}

tsmatrix * diagonal_matrix_matrix_product(tvector *matrix1, tsmatrix *matrix2) {
/*	if (matrix1->ncol != matrix2->nrow) {
		printf("Matrices size must agree (matrix_matrix_product error).\n");
		exit(1);
	}
*/
	tsmatrix * matrix_result = init_square_matrix(matrix1->length, matrix2->ncol);
	
	for (int i = 0; i < matrix_result->nrow; i++) {
		for (int j = 0; j < matrix_result->ncol; j++) {
			matrix_result->data[i][j] = matrix1->data[i]*matrix2->data[i][j];
		}
	}
	
	return matrix_result;
}

tsmatrix * change_basis_eigvectors_to_identity(tvector *Ematrix, tspectrum *spectrum) {
/*	if (Ematrix->nrow != Ematrix->ncol) {
		printf("Input matrix must be square (change_basis_eigvectors_to_identity error).\n");
		exit(1);
	}
*/
	tsmatrix * Imatrix = init_square_matrix(Ematrix->length, Ematrix->length);
	tsmatrix * P = init_square_matrix(Ematrix->length, Ematrix->length);
	tsmatrix * P_inverse = init_square_matrix(Ematrix->length, Ematrix->length);
	
	for (int i = 0; i < P->nrow; i++) {
		for (int j = 0; j < P->ncol; j++) {
			P->data[i][j] = spectrum->eigvectors->data[j][i];
		}
	}
	
	tvector * P_array = init_vector(P->nrow*P->ncol);	// (ncol*nrow_Hmatrix vector) needed conversion
																			// from tsmatrix to float-array
	for (int i = 0; i < P->nrow; i++) {								// for each row
		for (int j = 0; j < P->ncol; j++) {							// and for each column
			P_array->data[i*P->nrow+j] = P->data[i][j];	// insert matrix elements into the float-array matrix_array->data
		}
	}
	
	tvector_int * ipvt = init_vector_int(P->nrow);
	int INFO = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, P->nrow, P->ncol, P_array->data, P->nrow, ipvt->data);
	if (INFO == 0) {
		INFO = LAPACKE_sgetri(LAPACK_ROW_MAJOR, P->nrow, P_array->data, P->nrow, ipvt->data);
		
		if (INFO > 0) {
			printf("the algorithm failed to compute inverse matrix (change_basis_eigvectors_to_identity error).\n");
			exit(1);
		}
	} else {
		printf("the algorithm failed to compute factorization of matrix (change_basis_eigvectors_to_identity error).\n");
		exit(1);
	}
	
	for (int i = 0; i < P_inverse->nrow; i++) {
		for (int j = 0; j < P_inverse->ncol; j++) {
			P_inverse->data[i][j] = P_array->data[i*P_inverse->nrow+j];
		}
	}
	
	tsmatrix * matrix_aux = init_square_matrix(Imatrix->nrow, Imatrix->ncol);
	matrix_aux = diagonal_matrix_matrix_product(Ematrix, P_inverse);
	Imatrix = matrix_matrix_product(P, matrix_aux);
	
	return Imatrix;
}

// MATRIX - ARRAY PRODUCT
float * matrix_array_product(tsmatrix *matrix, float *array) {
	float * result = malloc(matrix->nrow*sizeof(float));
	for (int i = 0; i < matrix->nrow; i++) {
		for (int j = 0; j < matrix->ncol; j++) {
			result[i] += matrix->data[i][j]*array[j];
		}
	}
	
	return result;
}

// SCALAR PRODUCT
float scalar_product(float *array_1, float *array_2, int length) {
	float result = 0;
	for (int i = 0; i < length; i++) {
		result += array_1[i]*array_2[i];
	}
	
	return result;
}

/* NOT SQUARE MATRIX
typedef struct {
	float ** data;		// (2D-array) matrix data
	int nrow;			// number of rows
	int * ncol;			// (2D-array) number of columns for each row
} tnsmatrix;
*/

// INITIALIZE NOT SQUARE MATRIX
tnsmatrix * init_not_square_matrix(int size_row, int *size_col_vector) {
	tnsmatrix * nsmatrix = malloc(sizeof(tnsmatrix));					// (tnsmatrix-malloc) create malloc of ns-matrix
	
	nsmatrix->nrow = size_row;											// number of rows equal to row size
	nsmatrix->ncol = malloc(nsmatrix->nrow*sizeof(float));				// (nrow*float-malloc) create malloc of columns
	nsmatrix->data = malloc(nsmatrix->nrow*sizeof(float *));			// (nrow*float_pointer-malloc) cr mll f row data
	
	for (int i = 0; i < nsmatrix->nrow; i++) {							// for each row
		nsmatrix->data[i] = malloc(size_col_vector[i]*sizeof(float));	// (i_ncol*float-malloc) cr mall of i-col data
		nsmatrix->ncol[i] = size_col_vector[i];							// number of i-col equal to i-col size
	}
	
	return nsmatrix;
}

/* SPECTRUM
typedef struct {
	float * eigvals;			// (1D-array) eigenvalues
	tsmatrix * eigvectors;		// (tsmatrix-malloc) eigenvectors
	int length;					// (int) number of eigenvalues
} tspectrum;
*/

// INITIALIZE SPECTRUM
tspectrum * init_spectrum(int n_eigvals, int dim) {
	tspectrum * spectrum = malloc(sizeof(tspectrum));				// malloc of tspectrum
	
	spectrum->length = n_eigvals;									// length = number of eigenvalues
	spectrum->eigvals = malloc(spectrum->length*sizeof(float));		// (length*float) malloc of eigvals array
	spectrum->eigvectors = init_square_matrix(spectrum->length, dim);	// smatrix of eigvectors (sz =lngth)
	spectrum->positions = malloc(spectrum->length*sizeof(int));
	
	return spectrum;
}

/* SPECTRA
typedef struct {
	tspectrum ** data;	// list of spectra
	int nspectra;		// number of spectra
} tspectra;
*/

// INITIALIZE SPECTRA
tspectra * init_spectra(int n_spectra, int *length_spectra, int *dim) {
	tspectra * spectra = malloc(sizeof(tspectra));
	
	spectra->nspectra = n_spectra;
	spectra->data = malloc(spectra->nspectra*sizeof(tspectrum *));
	
	for (int i = 0; i < spectra->nspectra; i++) {
		spectra->data[i] = init_spectrum(length_spectra[i], dim[i]);
	}
	
	return spectra;
}

//// SYSTEM PREPARATION ////

// BASIS DIMENSION
int basis_dimension(float N, float gamma) {
	int dim;										// (int) basis dimension
	dim = (int) ((gamma + 1)*(N + 1)*(N + 2)/2);	// basis dimension equal to [...]
	return dim;
}

// BASIS STATES
void basis_states(tbasis *basis, float N, float gamma) {
		
	int N0, N1, N2;
	
	int j = 0;
	for (int i = 0; i <= N; i++) {
		N0 = (int) N - i;
		N1 = 0;
		N2 = i;
		while (N0 >= 0) {
			for (int k = 0; k <= gamma; k++) {
				basis->data[j]->n[0] = (float) N0;
				basis->data[j]->n[1] = (float) N1;
				basis->data[j]->n[2] = (float) N2;
				basis->data[j]->nph = (float) k;
				j++;
			}
			N0--;
			N1++;
		}
	}
}

// HAMILTONIAN ACTION ON ONE STATE
float Hact(tstate *vi, tstate *vj, float N, float wf, float *wn, float D, float *f) {
	float H;
	H = 0;
	
	if (vi->n[0] == vj->n[0] && vi->n[1] == vj->n[1] && vi->n[2] == vj->n[2]) {
		if (vi->nph == vj->nph) {
			H += wf*vi->nph + D*(2*vi->nph + 1);
			for (int alpha = 0; alpha < 3; alpha++) {
				H += wn[alpha]*vi->n[alpha];
			}
		//} else if ((vi->nph + 1) == vj->nph) {
		//	H += E*sqrt((vi->nph + 1)/N);
		//} else if ((vi->nph - 1) == vj->nph) {
		//	H += E*sqrt(vi->nph/N);
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

// HAMILTONIAN ACTION ON THE ENTIRE BASIS
tsmatrix * H(tbasis *basis, float N, int dim, float wf, float *wn, float D, float *f) {
	tsmatrix * Hmatrix = init_square_matrix(basis->nstates, basis->nstates);				// (tsmatrix) H matrix
	for (int i = 0; i < Hmatrix->nrow; i++) {												// for each row
		for (int j = 0; j < Hmatrix->ncol; j++) {														// and for each column
			Hmatrix->data[i][j] = Hact(basis->data[i], basis->data[j], N, wf, wn, D, f);	// Hact on ch pair of states
		}
	}
	
	return Hmatrix;
}

//// EIGENVALUES AND EIGENVECTORS ////

// EIGENVALUES AND VECTORS OF ONE MATRIX
tspectrum * eigenvalues(tsmatrix *Hmatrix) {
	int dim = Hmatrix->nrow;
	tspectrum * spectrum = init_spectrum(dim, dim);		// (tspectrum) initialize spectrum
	int INFO = 1;			// (int) out-parameter of LAPACKE_ssyev function
								// = 0: successful exit
								// < 0: if INFO = -i, the i-th argument had an illegal value
								// > 0: if INFO = i, the algorithm failed to converge; i off-diagonal elements of an
									// intermediate tridiagonal form did not converge to zero
	tvector * matrix_array = init_vector(Hmatrix->ncol*Hmatrix->nrow);	// (ncol*nrow_Hmatrix vector) needed conversion
																			// from tsmatrix to float-array
	
	for (int i = 0; i < Hmatrix->nrow; i++) {								// for each row
		for (int j = 0; j < Hmatrix->ncol; j++) {										// and for each column
			matrix_array->data[i*Hmatrix->nrow+j] = Hmatrix->data[i][j];	// insert matrix elements into the
			//printf("%f	", matrix_array->data[i*Hmatrix->nrow+j]);					// float-array matrix_array->data
		}
	}
	
/*	for (int i = 0; i < Hmatrix->nrow*Hmatrix->ncol; i++) {
		printf("%f	", matrix_array->data[i]);
	}
	printf("\n");
*/	
	//LAPACKE_ssyev(int matrix_layout, char jobz, char uplo, lapack_int n, float* a, lapack_int lda, float* w);
		// (int) matrix_layout: indicates whether the input and output matrices are stored in row mayor order or column
			// major order, where:
				// = LAPACK_ROW_MAJOR: the matrices are stored in row major order (C/C++ default)
				// = LAPACK_COL_MAJOR: the matrices are stored in column major order (Fortran default)
		// (char) jobz: indicates the type of computation to be performed, where:
				// = 'N': eigenvalues only are computed
				// = 'V': eigenvalues and eigenvectors are computed
		// (char) uplo: indicates whether the upper or lower part of matrix A is referenced, where:
				// = 'U': the upper triangular part is referenced
				// = 'L': the lower triangular part is referenced
		// (lapack_int) n: is the order of matrix A used in computation (n >= 0)
		// (float_pointer) a: is the real symmetric or complex Hermitian matrix A of order n
		// (lapack_int) lda: is the leading dimension of the array specified for A (lda >= n)
		// (float_pointer) w: is the vector w of length n, which contains the computed eigenvalues in acending order
	INFO = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'V', 'L', dim, matrix_array->data, dim, spectrum->eigvals);
	if (INFO > 0) {
		printf("the algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
	
/*	for (int i = 0; i < Hmatrix->nrow*Hmatrix->ncol; i++) {
		printf("%f	", matrix_array->data[i]);
	}
	printf("\n");
*/	
	// copy eigenvectors from matrix_array->data[] to eigvector_matrix->data[][]
	for (int i = 0; i < Hmatrix->nrow; i++) {
		for (int j = 0; j < Hmatrix->ncol; j++) {
			spectrum->eigvectors->data[i][j] = matrix_array->data[i*Hmatrix->nrow+j];
		}
	}
	
	return spectrum;
}

// EIGENVECTORS PROBLEM
tsmatrix * eigvector_problem(tspectrum *spectrum, tsmatrix *Hmatrix) {
	tsmatrix * A_lambda_I = init_square_matrix(Hmatrix->nrow, Hmatrix->ncol);
	tsmatrix * result_matrix = init_square_matrix(spectrum->length, A_lambda_I->nrow);
	for (int k = 0; k < spectrum->length; k++) {
		for (int i = 0; i < Hmatrix->nrow; i++) {
			for (int j = 0; j < Hmatrix->ncol; j++) {
				A_lambda_I->data[i][j] = Hmatrix->data[i][j];
			}
			A_lambda_I->data[i][i] -= spectrum->eigvals[k];
		}
	//tsmatrix * result_matrix = init_square_matrix(spectrum->length, A_lambda_I->nrow);
	//printf("%d %d\n", spectrum->eigvectors->nrow, result_matrix->nrow);
	//for (int i = 0; i < spectrum->length; i++) {
		result_matrix->data[k] = matrix_array_product(A_lambda_I, spectrum->eigvectors->data[k]);
	}
	
	return result_matrix;
}

// EIGENVECTORS ORTHONORMALITY
tsmatrix * eigvector_orthonormality(tspectrum *spectrum) {
	tsmatrix * result = init_square_matrix(spectrum->length, spectrum->length);
	
	for (int i = 0; i < spectrum->length; i++) {
		for (int j = 0; j < spectrum->length; j++) {
			result->data[i][j] = scalar_product(spectrum->eigvectors->data[i], spectrum->eigvectors->data[j], spectrum->eigvectors->nrow);
		}
	}
	
	return result;
}

// EIGENVALUES SERIES OF N SAMPLES
tspectra * eigval_series(float initial_gamma, int nsamples, float N, float wf, float *wn, float D, float *f) {
	float * g = malloc(nsamples*sizeof(float));		// (nsamples*float) malloc of n-photons array
	int * d = malloc(nsamples*sizeof(int));			// (nsamples*int) malloc of n-dimensions array
	
	g[0] = initial_gamma;							// initial value of number of photons
	d[0] = basis_dimension(N, g[0]);				// initial value of basis dimension
	for (int i = 0; i < (nsamples-1); i++) {		// for each sample
		g[i+1] = round(g[i]*1.1);					// next number of photons = round(1.1*previous)
		d[i+1] = basis_dimension(N, g[i+1]);		// compute dimension for new number of photons
	}
	
	tbasis * basis;														// (tbasis) malloc of system basis
	tsmatrix * Hmatrix;													// (tsmatrix) malloc of hamiltonian action
	tspectra * spectra = init_spectra(nsamples, d, d);				// (tspectra) initialize list of spectra (d-dim)
/*	tsmatrix ** eigvector_problem_result = malloc(nsamples*sizeof(tsmatrix*));
	tsmatrix ** eigvector_orthonormality_result = malloc(nsamples*sizeof(tsmatrix*));
	
	FILE *fp_p, *fp_o;
	fp_p = fopen("data/eigvector_problem.txt", "w");
	fp_o = fopen("data/eigvector_ortho.txt", "w");
*/
	for (int i = 0; i < nsamples; i++) {					// for each sample
		basis = init_basis(d[i]);							// initialize basis system (d[i]-dim)
		basis_states(basis, N, g[i]);						// compute basis states
		Hmatrix = init_square_matrix(d[i], d[i]);			// initialize hamiltonian actiom
		Hmatrix = H(basis, N, d[i], wf, wn, D, f);		// compute hamiltonian action
		spectra->data[i] = eigenvalues(Hmatrix);			// compute eigvalues and vectors
	}
/*		printf("MATRIX:	%d, %d\n", Hmatrix->nrow, Hmatrix->ncol);
		for (int k = 0; k < Hmatrix->nrow; k++) {
			for (int j = 0; j < Hmatrix->ncol; j++) {
				printf("%f ", Hmatrix->data[k][j]);
			}
			printf("\n");
		}
		printf("================================\n");
	
		
		eigvector_problem_result[i] = eigvector_problem(spectra->data[i], Hmatrix);
	
		//printf("%d, %d\n", eigvector_problem_result[i]->nrow, eigvector_problem_result[i]->ncol);
		for (int k = 0; k < eigvector_problem_result[i]->nrow; k++) {
			for (int j = 0; j < eigvector_problem_result[i]->ncol; j++) {
				//printf("%f ", eigvector_problem_result[i]->data[k][j]);
				fprintf(fp_p, "%f	", eigvector_problem_result[i]->data[k][j]);
			}
			//printf("\n");
			fprintf(fp_p, "\n");
		}
	}
		//printf("================================\n");
		
		eigvector_orthonormality_result[i] = eigvector_orthonormality(spectra->data[i]);
		
		//printf("%d, %d\n", eigvector_orthonormality_result[i]->nrow, eigvector_orthonormality_result[i]->ncol);
		for (int k = 0; k < eigvector_orthonormality_result[i]->nrow; k++) {
			for (int j = 0; j < eigvector_orthonormality_result[i]->ncol; j++) {
				//printf("%f ", eigvector_orthonormality_result[i]->data[k][j]);
				fprintf(fp_o, "%f	", eigvector_orthonormality_result[i]->data[k][j]);
			}
			//printf("\n");
			fprintf(fp_o, "\n");
		}
		fprintf(fp_p, "\n");
		fprintf(fp_o, "\n");	
	}
	
	fclose(fp_p);
	fclose(fp_o);
*/
	return spectra;
}

tspectra * fast_eigval_series(float initial_gamma, float N, float wf, float *wn, float D, float *f) {
	float * g = malloc(2*sizeof(float));		// (nsamples*float) malloc of n-photons array
	g[0] = initial_gamma;
	g[1] = round(initial_gamma*1.1);
	int * d = malloc(2*sizeof(int));
	d[0] = basis_dimension(N, g[0]);
	d[1] = basis_dimension(N, g[1]);
	tspectra * spectra = init_spectra(2, d, d);
	
	tbasis * basis;
	tsmatrix * Hmatrix;
	
	for (int i = 0; i < 2; i++) {					// for each sample
		basis = init_basis(d[i]);							// initialize basis system (d[i]-dim)
		basis_states(basis, N, g[i]);						// compute basis states
		Hmatrix = H(basis, N, d[i], wf, wn, D, f);		// compute hamiltonian action
		spectra->data[i] = eigenvalues(Hmatrix);			// compute eigvalues and vectors
		free(basis);
		free(Hmatrix);
	}
	
	return spectra;
}

// CONVERGENCE ANALYSIS
tvector_int * convergence_analysis(tspectra *spectra, float tol) {
	tvector_int * converg_analysis_array = init_vector_int(spectra->nspectra-1);
	int count;		// (int) counter of converged eigenvalues
	
	for (int i = 0; i < (spectra->nspectra-1); i++) {	// for each sample
		count = 0;											// initialize counter to zero
		for (int j = 0; j < spectra->data[i]->length; j++) {	// for each eigenvalue
			if (fabs(spectra->data[i+1]->eigvals[j] - spectra->data[i]->eigvals[j]) < tol) {
				count++;	// if the difference of eigvals of adjacent samples is less than a tolerance, increase count
			}
		}
		converg_analysis_array->data[i] = count;	// number of eigvals converged
	}
	
	return converg_analysis_array;
}

// CONVERGED EIGENVALUES
	// (int) convergence_type: type of computation for converged eigenvalues, where:
		// = 0: final configuration of the system
		// = 1: every configuration of the system
tspectra * converged_eigvals(tspectra *spectra, tvector_int *converg_analysis_array, float tol, int convergence_type) {
	tspectra * converged_spectra;			// (tspectra) converged eigenvalues and eigenvectors
	if (convergence_type == 0) {			// initialize conv spectra: 1 spectrum, length = number of conv eigvals
		converged_spectra = init_spectra(1, &(converg_analysis_array->data[converg_analysis_array->length-1]), &(spectra->data[spectra->nspectra-1]->length));
		int j = 0;	// id of last converged eigenvalue
		for (int i = 0; i <	 spectra->data[spectra->nspectra-2]->length; i++) {	// for each eigval of first config
			if (fabs(spectra->data[spectra->nspectra-1]->eigvals[i] - spectra->data[spectra->nspectra-2]->eigvals[i]) < tol) {
				converged_spectra->data[0]->eigvals[j] = spectra->data[spectra->nspectra-1]->eigvals[i];
				copy_row(converged_spectra->data[0]->eigvectors, spectra->data[spectra->nspectra-1]->eigvectors, j, i);
				converged_spectra->data[0]->positions[j] = i;
				j++;		// next id
			}
		}
	} else if (convergence_type == 1) {
		tvector_int * dim_spectra = init_vector_int(spectra->nspectra-1);
		for (int i = 0; i < spectra->nspectra-1; i++) {
			dim_spectra->data[i] = spectra->data[i]->length;
		}
		converged_spectra = init_spectra(spectra->nspectra-1, converg_analysis_array->data, dim_spectra->data);
		int j;
		for (int k = 0; k < converged_spectra->nspectra; k++) {	// for each configuration
			j = 0;	// id of last converged eigenvalue
			for (int i = 0; i < converged_spectra->data[k]->length; i++) {
				if (fabs(spectra->data[k+1]->eigvals[i] - spectra->data[k]->eigvals[i]) < tol) {
				converged_spectra->data[k]->eigvals[j] = spectra->data[k+1]->eigvals[i];	// add next conv eigval
				copy_row(converged_spectra->data[k]->eigvectors, spectra->data[k+1]->eigvectors, j, i);
				j++;		// next id
				}
			}
		}
	} else {
		printf("Convergence type error. Values admitted: 0 for final configuration and 1 for total configuration.\n");
		exit(1);
	}
	
	return converged_spectra;
}

tspectrum * fast_converged_spectrum(tspectra *spectra, float tol) {
	int j = 0;
	while ((fabs(spectra->data[1]->eigvals[j] - spectra->data[0]->eigvals[j]) < tol) && (j < spectra->data[0]->length)) {
		j++;
	}
	
	tspectrum * converged_spectrum = init_spectrum(j, spectra->data[1]->length);
	for (int i = 0; i < j; i++) {
		converged_spectrum->eigvals[i] = spectra->data[1]->eigvals[i];
		copy_row(converged_spectrum->eigvectors, spectra->data[1]->eigvectors, i, i);
	}
	
	return converged_spectrum;
}

//// STATISTICS ANALYSIS ////

// EIGENVALUES SPACINGS
tvector * eigval_spacing(tspectrum *spectrum) {
	tvector * spacing = init_vector(spectrum->length-1);					// (tvector) initialize spacings array
	for (int i = 0; i < spacing->length; i++) {								// for each pair of eigenvalues
		spacing->data[i] = spectrum->eigvals[i+1] - spectrum->eigvals[i];	// compt difference
	}
	
	return spacing;
}

// EIGENVALUES RATIOS
tvector * eigval_ratio(tvector *spacing) {
	tvector * ratio = init_vector(spacing->length-1);			// (tvector) initialize ratios array
	for (int i = 0; i < ratio->length; i++) {					// for each pair of spacings
		ratio->data[i] = spacing->data[i+1]/spacing->data[i];	// compute division
	}
	
	return ratio;
}

// EIGENVALUES MINIMAL RATIOS
tvector * eigval_min_ratio(tvector *ratio) {
	tvector * min_ratio = init_vector(ratio->length);			// (tvector) initialize min-ratios array
	tvector * inverse_ratio = init_vector(min_ratio->length);	// (tvector) initialize inverse ratios array
	
	for (int i = 0; i < min_ratio->length; i++) {				// for each ratio
		inverse_ratio->data[i] = 1./ratio->data[i];				// calculate inverse
		if (ratio->data[i] <= inverse_ratio->data[i]) {			// choose minimum between ratio and its inverse
			min_ratio->data[i] = ratio->data[i];
		} else {
			min_ratio->data[i] = inverse_ratio->data[i];
		}
	}
	
	return min_ratio;
}

// ENERGY WINDOW ANALYSIS
tsmatrix * window_analysis(tspectrum *spectrum, tvector *min_ratio, float num_windows, float overlap) {
	tsmatrix * end_points_matrix = init_square_matrix(num_windows, 2);
	int j;
	float E_0 = spectrum->eigvals[0];				// initial energy
	float length_window = (E_0 - (overlap*(num_windows - 1)))/num_windows;
	for (int i = 0; i < num_windows; i++) {
		j = 0;
		while (spectrum->eigvals[j] < (E_0 + i*(length_window + overlap))) {
			j++;
		}
		end_points_matrix->data[i][0] = j;
		while (spectrum->eigvals[j] < (E_0 + (i + 1)*length_window + (i + 2)*overlap)) {
			j++;
		}
		end_points_matrix->data[i][1] = j;
	}
	tsmatrix * mean_window_min_ratio = init_square_matrix(num_windows, 3);	// initialize matrix of w-min-ratios
																			// n-rows = number of windows
																			// n-cols = 3
																				// 1: mean of min-ratios in each window
																				// 2: starting energy of the window
																				// 3: last energy of the window
	
	for (int i = 0; i < num_windows; i++) {
		mean_window_min_ratio->data[i][0] = mean_vector(min_ratio, end_points_matrix->data[i][0], end_points_matrix->data[i][1]);
		mean_window_min_ratio->data[i][1] = spectrum->eigvals[(int)end_points_matrix->data[i][0]];
		mean_window_min_ratio->data[i][2] = spectrum->eigvals[(int)end_points_matrix->data[i][1]];
	}
	
	return mean_window_min_ratio;
}

// INVERSE PARTICIPATION RATIO (IPR)
float * inverse_participation_ratio(tspectrum *spectrum) {
	float * IPR = malloc(spectrum->length*sizeof(float));			// (eigvectors*float) malloc of IPR
	float sum;														// sum part of IPR (denominator)
	for (int i = 0; i < spectrum->length; i++) {					// for each eigenvector
		sum = 0;													// reload sum
		for (int j = 0; j < spectrum->eigvectors->ncol; j++) {				// for each element of each eigenvector
			sum += pow(spectrum->eigvectors->data[i][j], 4);		// increase sum with (element)^4
		}
		IPR[i] = 1./sum;												// IPR = inverse of sum
	}
	return IPR;
}

tmixstate * array_to_state(tbasis *basis, float *vector) {	
	tmixstate * mix_state = init_mix_state(basis->nstates);
	for (int i = 0; i < basis->nstates; i++) {
		mix_state->width[i] = vector[i];
		for (int j = 0; j < 3; j++) {
			mix_state->data[i]->n[j] = basis->data[i]->n[j];
		}
		mix_state->data[i]->nph = basis->data[i]->nph;
	}
	
	return mix_state;
}

float * op_sigma(tmixstate *mix_state) {
	float * sigma = malloc(3*sizeof(float));
	//sigma[0] = 0;
	//sigma[1] = 0;
	//sigma[2] = 0;
	for (int i = 0; i < mix_state->nstates; i++) {
		sigma[0] += mix_state->width[i]*mix_state->width[i]*mix_state->data[i]->n[0];
		sigma[1] += mix_state->width[i]*mix_state->width[i]*mix_state->data[i]->n[1];
		sigma[2] += mix_state->width[i]*mix_state->width[i]*mix_state->data[i]->n[2];
	}
	
	return sigma;
}

float op_ada(tmixstate *mix_state) {
	float ada = 0;
	for (int i = 0; i < mix_state->nstates; i++) {
		ada += pow(mix_state->width[i], 2)*mix_state->data[i]->nph;
	}
	
	return ada;
}

// PERES LATTICES
void peres_lattice(tsmatrix *On_sigma, tvector *On_ada, tbasis *basis, tspectrum *spectrum) {
	//tsmatrix * On_sigma = init_square_matrix(spectrum->length, 3);
	//tvector * On_ada = init_vector(spectrum->length);
	tmixstate ** mix_states = malloc(spectrum->length*sizeof(tmixstate*));

	for (int i = 0; i < spectrum->length; i++) {
		mix_states[i] = array_to_state(basis, spectrum->eigvectors->data[i]);
		On_sigma->data[i] = op_sigma(mix_states[i]);
		On_ada->data[i] = op_ada(mix_states[i]);
	}
}

float Hact_a_plus_ad(tstate *vi, tstate *vj) {
	float a_plus_ad = 0;
	
	if (vi->n[0] == vj->n[0] && vi->n[1] == vj->n[1] && vi->n[2] == vj->n[2]) {
		if ((vi->nph - 1) == vj->nph) {
			a_plus_ad += sqrt(vi->nph);
		} else if ((vi->nph + 1) == vj->nph) {
			a_plus_ad += sqrt(vi->nph + 1);
		}
	}
	return a_plus_ad;
}

tsmatrix * H_a_plus_ad(tbasis *basis) {
	tsmatrix * Hmatrix_a_plus_ad = init_square_matrix(basis->nstates, basis->nstates);		// (tsmatrix) H matrix
	for (int i = 0; i < Hmatrix_a_plus_ad->nrow; i++) {										// for each row
		for (int j = 0; j < Hmatrix_a_plus_ad->ncol; j++) {									// and for each column
			Hmatrix_a_plus_ad->data[i][j] = Hact_a_plus_ad(basis->data[i], basis->data[j]);	// Hact on ch pair of states
		}
	}
	
	return Hmatrix_a_plus_ad;
}

float sign(float number) {
	float sgn = 0;
	if (number > 0) {
		sgn = 1;
	} else if (number < 0) {
		sgn = -1;
	}
	
	return sgn;
}

tvector * peres_lattice_sign_a_plus_ad(tbasis *basis, tspectrum *spectrum) {
	tsmatrix * Hmatrix_a_plus_ad = H_a_plus_ad(basis);					// 1. M = a + ad
	tspectrum * spectrum_a_plus_ad = eigenvalues(Hmatrix_a_plus_ad);	// 2. Diagonalize M
	tvector * Ematrix_sign_a_plus_ad = init_vector(Hmatrix_a_plus_ad->nrow);
	
	for (int i = 0; i < spectrum_a_plus_ad->length; i++) {
		Ematrix_sign_a_plus_ad->data[i] = sign(spectrum_a_plus_ad->eigvals[i]);	// 3. C in diagonal basis of M
	}
	
	tsmatrix * Imatrix_sign_a_plus_ad = change_basis_eigvectors_to_identity(Ematrix_sign_a_plus_ad, spectrum_a_plus_ad);	// 4. Return to original basis
	
	tvector * C_expected_value = init_vector(spectrum->length);
	tvector * vector_aux = init_vector(Imatrix_sign_a_plus_ad->nrow);
	for (int i = 0; i < spectrum->length; i++) {
		vector_aux->data = matrix_array_product(Imatrix_sign_a_plus_ad, spectrum->eigvectors->data[i]);
		C_expected_value->data[i] = scalar_product(vector_aux->data, spectrum->eigvectors->data[i], vector_aux->length);
	}
	
	return C_expected_value;
}

tsmatrix * critical_energy(tspectrum *spectrum, tvector *C_expected_value, float tol) {
	tsmatrix * Ec = init_square_matrix(3,2);
	int j = 0;
	while ((fabs(1. - fabs(C_expected_value->data[j])) < tol) && (j < C_expected_value->length)) {
		j++;
	}
	
	for (int i = 0; i < Ec->nrow; i++) {
		Ec->data[i][0] = spectrum->eigvals[j-2+i];
		Ec->data[i][1] = C_expected_value->data[j-2+i];
	}
	
	return Ec;
}

tvector * number_of_photons(tbasis *basis, tspectrum *spectrum) {
	tvector * nphotons = init_vector(spectrum->length);
	
	for (int i = 0; i < spectrum->length; i++) {
		tmixstate * aux_mixstate = array_to_state(basis, spectrum->eigvectors->data[i]);
		for (int j = 0; j < aux_mixstate->nstates; j++) {
			nphotons->data[i] += pow(aux_mixstate->width[i], 2) * aux_mixstate->data[i]->nph;
		}
	}
	
	return nphotons;
}