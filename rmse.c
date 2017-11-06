#include <stdio.h>
#include <stdlib.h>
#include <mat.h>
#include <mkl.h>
#include <math.h>

double rmse(const mxArray *mx_truth, const double *estimated, const int m, const int n, const int k)
{
	const double *truth = (const double *)mxGetPr(mx_truth);
	
	double *copy_truth = (double *)mkl_malloc(m*n*k * sizeof(double), 32);
	if (copy_truth == NULL) {
		printf("Error: can not allocate memory for copy_truth. Aborting ......\n");
		return -1;
	}
	double max_truth = 0;
	for (int i = 0; i < m*n*k; i++) {
		copy_truth[i] = truth[i];
		if (truth[i] > max_truth) {
			max_truth = truth[i];
		}
	}
	for (int i = 0; i < m*n*k; i++) {
		copy_truth[i] *= (255 / max_truth);
	}

	double max_estimated = 0;
	double *copy_estimated = (double *)mkl_malloc(m*n*k * sizeof(double), 32);
	for (int i = 0; i < m*n*k; i++) {
		if (estimated[i] > 1) {
			copy_estimated[i] = 1;
		}
		else if (estimated[i] < 0) {
			copy_estimated[i] = 0;
		}
		else {
			copy_estimated[i] = estimated[i];
		}
		
		if (copy_estimated[i] > max_estimated) {
			max_estimated = copy_estimated[i];
		}
	}
	for (int i = 0; i < m*n*k; i++) {
		copy_estimated[i] *= 255;
	}


	double err_sum = 0;
	for (int i = 0; i < m*n*k; i++) {
		err_sum += (copy_estimated[i] - copy_truth[i]) * (copy_estimated[i] - copy_truth[i]);
	}

	mkl_free(copy_truth);
	mkl_free(copy_estimated);

	return sqrt(err_sum / (m*n*k));
}
