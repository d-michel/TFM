#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include "include/def.h"
#include <time.h>

int main(int argc, char **argv) {
	
	time_t t_0, t_1;//, t_2, t_3;
	//struct tm *local_time;
	t_0 = time(NULL);
	//t_1 = time(NULL);
	
	float N = 10;
	float gamma = 150;
	float wf = 1;
	float * wn = malloc(3*sizeof(float));
	wn[0] = 1;
	wn[1] = 1.1;
	wn[2] = 2.1;
	float D = 3;
	float * f = malloc(3*sizeof(float));
	f[0] = 0;
	f[1] = 6;
	f[2] = 0;
	
	float tol = 0.1;
	tspectra * spectra;
	tspectrum * converged_spectrum;
	spectra = fast_eigval_series(gamma, N, wf, wn, D, f);
	converged_spectrum = fast_converged_spectrum(spectra, tol);
	
/*	int dim = basis_dimension(N, round(1.1*gamma));
	tbasis * basis = init_basis(dim);
	basis_states(basis, N, round(1.1*gamma));
	tsmatrix * On_sigma = init_square_matrix(converged_spectrum->length, 3);
	tvector * On_ada = init_vector(converged_spectrum->length);
	peres_lattice(On_sigma, On_ada, basis, converged_spectrum);
	tvector * On_sign_a_plus_ad = peres_lattice_sign_a_plus_ad(basis, converged_spectrum);
	FILE *fp_peres;
	fp_peres = fopen("peres_lattices_150_10_1_6_sign.txt", "w");
	for (int j = 0; j < converged_spectrum->length; j++) {
				//	  eig	ipr	00	11	22	ada
		fprintf(fp_peres, "%f	%f	%f	%f	%f	%f\n", converged_spectrum->eigvals[j],
														On_sigma->data[j][0],
														On_sigma->data[j][1],
														On_sigma->data[j][2],
														On_ada->data[j],
														On_sign_a_plus_ad->data[j]);
	}
	fclose(fp_peres);
*/
/*	tvector * IPR = init_vector(converged_spectrum->length);
	IPR->data = inverse_participation_ratio(converged_spectrum);
	FILE *fp_ipr;
	fp_ipr = fopen("data/peres_lattices/IPR_150_10_1_6.txt", "w");
	for (int i = 0; i < IPR->length; i++) {
		fprintf(fp_ipr, "%f	%f\n", converged_spectrum->eigvals[i], IPR->data[i]);
	}
	fclose(fp_ipr);
*/
/*	tsmatrix * density = density_of_states(converged_spectrum, 15);
	FILE *fp_density;
	fp_density = fopen("data/density/150_6_03_6_15.txt", "w");
	for (int i = 0; i < density->nrow; i++) {
			fprintf(fp_density, "%f	%f\n", density->data[i][0], density->data[i][1]);
	}
	fclose(fp_density);
*/
	int dim = basis_dimension(N, round(gamma*1.1));
	tbasis * basis = init_basis(dim);
	basis_states(basis, N, round(gamma*1.1));
	tvector * On_a_plus_ad = peres_lattice_a_plus_ad(basis, converged_spectrum);
	//tvector * C_expected_value = peres_lattice_sign_a_plus_ad(basis, converged_spectrum);
	
	FILE *fp_Ec;
	fp_Ec = fopen("symmetries_150_10_1_6_f0f2.txt", "w");
	for (int j = 0; j < converged_spectrum->length; j++) {
		fprintf(fp_Ec, "%f	%f\n", converged_spectrum->eigvals[j], On_a_plus_ad->data[j]);
	}
	fclose(fp_Ec);

/*	tvector * spacing = eigval_spacing(converged_spectrum);
	tvector * ratio = eigval_ratio(spacing);
	tvector * min_ratio = eigval_min_ratio(ratio);
    FILE *fp_statistics;
	fp_statistics = fopen("statistics_150_10_1_6.txt", "w");
	for (int i = 0; i < ratio->length; i++) {
		fprintf(fp_statistics, "%f	%f	%f\n", converged_spectrum->eigvals[i], ratio->data[i], min_ratio->data[i]);
	}
	fclose(fp_statistics);
*/
/*	
	int dim = basis_dimension(N, round(1.1*gamma));
	tbasis * basis = init_basis(dim);
	basis_states(basis, N, round(1.1*gamma));
	
	tvector * IPR = init_vector(converged_spectrum->length);
	
	IPR->data = inverse_participation_ratio(converged_spectrum);
	
	tsmatrix * On_sigma = init_square_matrix(converged_spectrum->length, 3);
	tvector * On_ada = init_vector(converged_spectrum->length);
	
	peres_lattice(On_sigma, On_ada, basis, converged_spectrum);
	
	FILE *fp;
	fp = fopen("data/peres_lattices/peres_lattices_150_5_wf05.txt", "w");
	for (int j = 0; j < converged_spectrum->length; j++) {
				//	  eig	ipr	00	11	22	ada
		fprintf(fp, "%f	%f	%f	%f	%f	%f\n", converged_spectrum->eigvals[j],
														IPR->data[j],
														On_sigma->data[j][0],
														On_sigma->data[j][1],
														On_sigma->data[j][2],
														On_ada->data[j]);
	}
	fclose(fp);
*/
	t_1 = time(NULL);
	printf("Timelapse: %ld seconds\n", (t_1 - t_0));
	
	return 0;
}

//	tsmatrix * On_sigma = init_square_matrix(converged_spectra->data[0]->length, 3);
//	tvector * On_ada = init_vector(converged_spectra->data[0]->length);
//	peres_lattice(On_sigma, On_ada, basis, converged_spectra->data[0]);
//	tvector * nphotons = number_of_photons(basis, converged_spectra->data[0]);

	// ANÁLISIS DE CONVERGENCIA
//	tvector_int * converg_analysis_array = convergence_analysis(spectra, tol);
	
	// CÁLCULO DE AUTOVALORES CONVERGIDOS
//	tspectra * converged_spectra = converged_eigvals(spectra, converg_analysis_array, tol, 0);
	
/*	// ANÁLISIS ESPECTRAL
	tvector * spacing = eigval_spacing(converged_spectra->data[0]);
	tvector * ratio = eigval_ratio(spacing);
	tvector * min_ratio = eigval_min_ratio(ratio);
	float num_windows = 10;
	float overlap = 10;
	tsmatrix * window_mean_min_ratio = window_analysis(converged_spectra->data[0], min_ratio, num_windows, overlap);
	
	for (int i = 0; i < ratio->length; i++) {
		fprintf(fp_ratio, "%f	%f	%f\n", converged_spectra->data[0]->eigvals[i], ratio->data[i], min_ratio->data[i]);
	}
	for (int i = 0; i < window_mean_min_ratio->nrow; i++) {
		fprintf(fp_window, "%f	%f	%f\n", window_mean_min_ratio->data[i][0], window_mean_min_ratio->data[i][1], window_mean_min_ratio->data[i][2]);
	}
	fclose(fp_ratio);
	fclose(fp_window);
*/
	// CÁLCULO DE INVERSE PARTICIPATION RATIO
//	float * IPR = inverse_participation_ratio(converged_spectra->data[0]);
	
/*	// GUARDADO DE DATOS
	fprintf(fp0, "N	gamma	wf	wn	D	f\n");
	fprintf(fp0, "%.0f	%.0f	%.1f	%.1f	%.1f	%.1f\n", N, gamma, wf, wn[0], D, f[0]);
	fprintf(fp0, "			%.1f		%.1f\n", wn[1], f[1]);
	fprintf(fp0, "			%.1f		%.1f\n", wn[2], f[2]);
	for (int j = 0; j < converged_spectra->data[0]->length; j++) {
				//	 pos	eig	ipr	00	11	22	ada
		fprintf(fp1, "%d	%f	%f	%f	%f	%f	%f\n", converged_spectra->data[0]->positions[j],
														converged_spectra->data[0]->eigvals[j],
														IPR[j],
														On_sigma->data[j][0],
														On_sigma->data[j][1],
														On_sigma->data[j][2],
														On_ada->data[j]);
	}
	for (int i = 0; i < window_mean_min_ratio->nrow; i++) {
		fprintf(fp2, "%f	%f	%f\n", window_mean_min_ratio->data[i][0], window_mean_min_ratio->data[i][1], window_mean_min_ratio->data[i][2]);
	}
	for (int i = 0; i < ratio->length; i++) {
		fprintf(fp3, "%f	%f	%f\n", spacing->data[i], ratio->data[i], min_ratio->data[i]);
	}
	fprintf(fp3, "%f", spacing->data[spacing->length-1]);
*/	

	/* IMPRESIÓN DE OBJETOS
	// Print basis
	for (int i = 0; i < basis->nvec; i++) {
		printf("%.0f, %.0f, %.0f, %.0f\n", basis->data[i]->n[0], basis->data[i]->n[1], basis->data[i]->n[2], basis->data[i]->nph);
	}
	
	// Print matrix
	printf("%d, %d\n", Hmatrix->nrow, Hmatrix->ncol);
	for (int i = 0; i < Hmatrix->nrow; i++) {
		for (int j = 0; j < Hmatrix->ncol; j++) {
			printf("%f ", Hmatrix->data[i][j]);
		}
		printf("\n");
	}
	
	// Print eigenvalues
	for (int i = 0; i < dim; i++) {
		printf("%f\n", eigvals[i]);
	}
	*/
	
	//////////////////////////////////////////////////////////////
	
	/* ANÁLISIS ESPECTRAL PARA DIFERENTES VALORES DE LOS PARÁMETROS f_ij
	int nsamples = 2;
	float tol = 0.1;
	int length_window = 700;
	int overlap = 50;
	
	int nfsamples = 10;
	
	tvector * mean_min_ratio_fseries = init_vector(nfsamples);
	
	for (int i = 0; i < nfsamples; i++) {
		// CÁLCULO DE AUTOVALORES DE LA SERIE DE CONFIGURACIONES
		tnsmatrix * eigvals_matrix = eigval_series(gamma, nsamples, N, wf, wn, D, f, E);
		
		// ANÁLISIS DE CONVERGENCIA
		tsmatrix * converg_analysis_matrix = convergence_analysis(eigvals_matrix, tol);
		
		// CÁLCULO DE AUTOVALORES CONVERGIDOS
		tnsmatrix * eigvals = converged_eigvals(eigvals_matrix, converg_analysis_matrix, tol, 0);
		
		// ANÁLISIS DE LOS ESTADÍSTICOS ESPECTRALES
		tvector * spacing = init_vector(eigvals->ncol[0]-1);						// spacing := vector de espaciamientos
		tvector * ratio = init_vector(spacing->length-1);						// ratio := vector de ratios
		tvector * min_ratio = init_vector(spacing->length-1);					// min_ratio := vector de minratios
		//float mean_min_ratio;
		
		eigval_spacing(eigvals, spacing);							// calculamos espaciamientos
		eigval_ratio(spacing, ratio);								// calculamos ratios
		eigval_min_ratio(ratio, min_ratio);							// calculamos minratios
		mean_min_ratio_fseries->data[i] = mean_vector(min_ratio, 0, min_ratio->length);
		printf("%f %f\n", mean_min_ratio_fseries->data[i], f[0]);
		for (int j = 0; j < 3; j++) {
			f[j] = f[j]*1.01;
		}
	}
	*/
	
	//////////////////////////////////////////////////////////////
	
	/* ANÁLISIS ESPECTRAL POR VENTANAS DE ENERGÍA
	printf("%f\n", mean_min_ratio);
	
	// ANÁLISIS DE VENTANA
	tsmatrix * mean_window_min_ratio = window_analysis(min_ratio, length_window, overlap);
	
	
	for (int i = 0; i < converg_analysis_matrix->nrow; i++) {			// se imprimen los contadores
		for (int j = 0; j < converg_analysis_matrix->ncol; j++) {
			printf("%f ", converg_analysis_matrix->data[i][j]);
		}
		printf("\n");
	}
	
	printf("%d\n", eigvals->ncol[0]);							// se imprime el número de autovalores convergidos de la última configuración
	
	printf("%d, %d\n", mean_window_min_ratio->nrow, mean_window_min_ratio->ncol);
	for (int i = 0; i < mean_window_min_ratio->nrow; i++) {
		for (int j = 0; j < mean_window_min_ratio->ncol; j++) {
			printf("%f ", mean_window_min_ratio->data[i][j]);
		}
		printf("\n");
	}
	*/
	
	//////////////////////////////////////////////////////////////
	
	/* GUARDADO DE DATOS
	// GUARDADO DE DATOS
	FILE *fp;													// guardamos los resultados en un fichero txt
	fp = fopen("data_f20.txt", "w");
	//fprintf(fp, "eigenvalues	spacings	ratios	min_ratios\n");
	for (int i = 0; i < ratio->length; i++) {
		fprintf(fp, "%f	%f	%f	%f\n", eigvals->data[0][i], spacing->data[i], ratio->data[i], min_ratio->data[i]);
	}
	fprintf(fp, "%f	%f\n", eigvals->data[0][eigvals->ncol[0]-2], spacing->data[eigvals->ncol[0]-1]);
	fprintf(fp, "%f\n", eigvals->data[0][eigvals->ncol[0]-1]);
	fclose(fp);
	
	fp = fopen("window_700_f20.txt", "w");
	for (int i = 0; i < mean_window_min_ratio->nrow; i++) {
		fprintf(fp, "%f	%f	%f\n", mean_window_min_ratio->data[i][0], mean_window_min_ratio->data[i][1], mean_window_min_ratio->data[i][2]);
	}
	fclose(fp);
	*/
	
	// PRUEBA DE AUTOVALORES CON LAPACKE_SSYEV()
/*	tsmatrix * matrix = init_square_matrix(3,3);
	matrix->data[0][0] = 1;
	matrix->data[0][1] = -2;
	matrix->data[0][2] = 3;
	matrix->data[1][0] = -2;
	matrix->data[1][1] = 2;
	matrix->data[1][2] = 4;
	matrix->data[2][0] = 3;
	matrix->data[2][1] = 4;
	matrix->data[2][2] = 3;
	
	
	tspectrum * spectrum = eigenvalues(matrix);
	tsmatrix * ev_problem = eigvector_problem(spectrum, matrix);
	tsmatrix * ev_orthonormality = eigvector_orthonormality(spectrum);
	
	printf("======== MATRIX ========\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f	", matrix->data[i][j]);
		}
		printf("\n");
	}
	
	printf("======== EIGENVALUES ========\n");
	for (int i = 0; i < 3; i++) {
		printf("%f	", spectrum->eigvals[i]);
	}
	printf("\n");
	
	printf("======== EIGENVECTORS ========\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f	", spectrum->eigvectors->data[i][j]);
		}
		printf("\n");
	}
	
	printf("======== EV PROBLEM ========\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f	", ev_problem->data[i][j]);
		}
		printf("\n");
	}
	
	printf("======== EV ORTHONORMALITY ========\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%f	", ev_orthonormality->data[i][j]);
		}
		printf("\n");
	}
*/

/*	// COMPROBACIÓN DE AUTOVALORES
	dim = basis_dimension(N, round(gamma*1.1));
	basis = init_basis(dim);
	tsmatrix * Hmatrix = H(basis, N, dim, wf, wn, D, f);
	tsmatrix * eigvector_problem_result = eigvector_problem(spectra->data[1], Hmatrix);
	
	printf("%d, %d\n", eigvector_problem_result->nrow, eigvector_problem_result->ncol);
	for (int k = 0; k < eigvector_problem_result->nrow; k++) {
		for (int j = 0; j < eigvector_problem_result->ncol; j++) {
			printf("%f ", eigvector_problem_result->data[k][j]);
		}
		printf("\n");
	}
	printf("================================\n");
	
	tsmatrix * eigvector_orthonormality_result = eigvector_orthonormality(spectra->data[1]);
	
	printf("%d, %d\n", eigvector_orthonormality_result->nrow, eigvector_orthonormality_result->ncol);
	for (int k = 0; k < eigvector_orthonormality_result->nrow; k++) {
		for (int j = 0; j < eigvector_orthonormality_result->ncol; j++) {
			printf("%f ", eigvector_orthonormality_result->data[k][j]);
		}
		printf("\n");
	}
*/



/*	// CANTIDAD CONSERVADA
	tspectra * spectra;
	tspectrum * converged_spectrum;
	int dim;
	tbasis * basis;
	tvector * C_expected_value;
	float Ec;
	int nNsamples = 1;
	float tol = 0.1;
	float tolN = 0.05;
	//float window_width = 2;
	//tsmatrix * density;

	for (int i = 0; i < nNsamples; i++) {
		//t_2 = time(NULL);
		// CÁLCULO DE AUTOVALORES DE LA SERIE DE CONFIGURACIONES
		spectra = fast_eigval_series(gamma, N, wf, wn, D, f);
		converged_spectrum = fast_converged_spectrum(spectra, tol);
		
		timebreak(&t_0, &t_1, "Converged spectrum calculated!");
		
		// CÁLCULO DE DENSIDAD DE NIVELES
		//density = density_of_states(converged_spectrum, window_width);
		//for (int j = 0; j < density->nrow; j++) {
			//fprintf(fp_density, "%f	%f\n", density->data[j][0], density->data[j][1]);
		//}
	
		// CÁLCULO DE PERES LATTICES
		dim = basis_dimension(N, round(gamma*1.1));
		basis = init_basis(dim);
		basis_states(basis, N, round(gamma*1.1));
		
		timebreak(&t_0, &t_1, "Basis states calculated!");
		
		C_expected_value = peres_lattice_sign_a_plus_ad(basis, converged_spectrum);
		
		timebreak(&t_0, &t_1, "C operator calculated!");
		
		Ec = critical_energy(converged_spectrum, C_expected_value, tolN);
		
		
		timebreak(&t_0, &t_1, "Critical energy calculated!");
		
		fprintf(fp_Ec, "\n%f	%f	%f	%f\n\n", N, Ec, converged_spectrum->eigvals[0], converged_spectrum->eigvals[converged_spectrum->length-1]);
		for (int j = 0; j < converged_spectrum->length; j++) {
			fprintf(fp_Ec, "%f	%f\n", converged_spectrum->eigvals[j], C_expected_value->data[j]);
		}
		//t_3 = time(NULL);
		//fprintf(fp_time, "%f	%ld\n", gamma, t_3-t_2);
		
		N += 1;
		free(spectra);
		free(converged_spectrum);
		free(basis);
		free(C_expected_value);
		
		printf("========== End of simulation %d of %d ==========\n\n", i+1, nNsamples);
	}
*/