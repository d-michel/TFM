#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <lapacke.h>

//// CLASS DEFINITIONS ////
typedef struct {
	float * n;
	float nph;
} tstate;

typedef struct {
	tstate ** data;
	float * width;
	int nstates;
} tmixstate;

typedef struct {
	tstate ** data;
	int nstates;
} tbasis;

typedef struct {
	float * data;
	int length;
} tvector;

typedef struct {
	int * data;
	int length;
} tvector_int;

typedef struct {
	float ** data;
	int nrow;
	int ncol;
} tsmatrix;

typedef struct {
	float ** data;
	int nrow;
	int * ncol;
} tnsmatrix;

typedef struct {
	float * eigvals;
	tsmatrix * eigvectors;
	int length;
	int * positions;
} tspectrum;

typedef struct {
	tspectrum ** data;
	int nspectra;
} tspectra;

//// CLASS OPERATIONS ////

tstate * init_state();

tmixstate * init_mix_state(int n_states);

tbasis * init_basis(int size_basis);

tvector * init_vector(int size_vector);

float mean_vector(tvector *vector, int start, int finish);

tvector_int * init_vector_int(int size_vector);

tsmatrix * init_square_matrix(int size_row, int size_col);

void copy_column(tsmatrix *out_matrix, tsmatrix *in_matrix, int id_out, int id_in);

void copy_row(tsmatrix *out_matrix, tsmatrix *in_matrix, int id_out, int id_in);

tsmatrix * matrix_matrix_product(tsmatrix *matrix1, tsmatrix *matrix2);

tsmatrix * diagonal_matrix_matrix_product(tvector *matrix1, tsmatrix *matrix2);

tsmatrix * change_basis_eigvectors_to_identity(tvector *Ematrix, tspectrum *spectrum);

float * matrix_array_product(tsmatrix *matrix, float *array);

float scalar_product(float *array_1, float *array_2, int length);

tnsmatrix * init_not_square_matrix(int size_row, int *size_col_vector);

tspectrum * init_spectrum(int n_eigvals, int dim);

tspectra * init_spectra(int n_spectra, int *n_eigvals, int *dim);

//// SYSTEM PREPARATION ////

int basis_dimension(float N, float gamma);

void basis_states(tbasis *basis, float N, float gamma);

float Hact(tstate *vi, tstate *vj, float N, float wf, float *wn, float D, float *f);	// without E parameter (26 feb)

tsmatrix * H(tbasis *basis, float N, int dim, float wf, float *wn, float D, float *f);	// without E parameter (26 feb)

//// EIGENVALUES AND EIGENVECTORS ////

tspectrum * eigenvalues(tsmatrix *Hmatrix);

tsmatrix * eigvector_problem(tspectrum *spectrun, tsmatrix *Hmatrix);

tsmatrix * eigvector_orthonormality(tspectrum *spectrum);

tspectra * eigval_series(float initial_gamma, int nsamples, float N, float wf, float *wn, float D, float *f); // without E parameter (26 feb)

tspectra * fast_eigval_series(float initial_gamma, float N, float wf, float *wn, float D, float *f);

tvector_int * convergence_analysis(tspectra *spectra, float tol);

tspectra * converged_eigvals(tspectra *spectra, tvector_int *converg_analysis_array, float tol, int convergence_type);

tspectrum * fast_converged_spectrum(tspectra *spectra, float tol);

//// STATISTICS ANALYSIS ////

tvector * eigval_spacing(tspectrum *spectrum);

tvector * eigval_ratio(tvector *spacing);

tvector * eigval_min_ratio(tvector *ratio);

tsmatrix * window_analysis(tspectrum *spectrum, tvector *min_ratio, float num_windows, float overlap);

float * inverse_participation_ratio(tspectrum *spectrum);

tmixstate * array_to_state(tbasis *basis, float *vector);

float * op_sigma(tmixstate *mix_state);

float op_ada(tmixstate *mix_state);

void peres_lattice(tsmatrix *On_sigma, tvector *On_ada, tbasis *basis, tspectrum *spectrum);

float Hact_a_plus_ad(tstate *vi, tstate *vj);

tsmatrix * H_a_plus_ad(tbasis *basis);

float sign(float number);

tvector * peres_lattice_sign_a_plus_ad(tbasis *basis, tspectrum *spectrum);

tsmatrix * critical_energy(tspectrum *spectrum, tvector *C_expected_value, float tol);

tvector * number_of_photons(tbasis *basis, tspectrum *spectrum);