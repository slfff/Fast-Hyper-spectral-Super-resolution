/* Hyperspectral Super-Resolution using Linear Regression 
 * C version
 * Lingfei Song 2017.8.16 */

#include <stdio.h>
#include <stdlib.h>
#include <mat.h>
#include <mkl.h>
#include <time.h>
#include <string.h>
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>

struct thread_arg {
	double *D_multi;
	double *D_hyper;
	double *multi;
	double *hyper;
	int i;
	int j;
	mwSize size_multi[2];
};

void *sub_thread(void *arg);

void hss_lr(const mxArray *mx_D_hyper, const mxArray *mx_multi, const mxArray *mx_subsample, double *estimated_hyper)
{
	const mwSize *size_multi = mxGetDimensions(mx_multi);
	const mwSize *size_subsample = mxGetDimensions(mx_subsample);

	int i, j;

	//printf("The dimension of D_hyper is %d %d %d.\n", size_D_hyper[0], size_D_hyper[1], size_D_hyper[2]);
	//printf("The dimension of multi is %d %d %d.\n", size_multi[0], size_multi[1], size_multi[2]);
	//printf("The dimension of subsample is %d %d.\n", size_subsample[0], size_subsample[1]);

	double *D_hyper = (double *)mxGetPr(mx_D_hyper);
	double *multi = (double *)mxGetPr(mx_multi);

	/* sparse matrix, compressed sparse column(csc) format. 
	 * pr is a double-precision array of length nzmax contains the nonzero values of the matrix;
	 * ir points to an integer array also of length nzmax contatining the row indices of the 
	 * corresponding elements in pr;
	 * jc points to an integer array of length n+1, where n is the number of columns in the sparse
	 * matrix. The jc array contains column index information. If the j-th column of the sparse matrix
	 * has any nonzero elements, jc[j] is the index of the first nonzero element in the j-th column, 
	 * and jc[j+1]-1 is the index of the last nonzero element in that column. */
	double *subsample_pr = (double *)mxGetPr(mx_subsample);
	mwIndex *subsample_ir = mxGetIr(mx_subsample);
	mwIndex *subsample_jc = mxGetJc(mx_subsample);
	mwSize subsample_nzmax = mxGetNzmax(mx_subsample);

	double *D_multi = (double *)mkl_malloc((size_multi[0] / 32)*(size_multi[1] / 32)*size_multi[2] * sizeof(double), 64);
	if (D_multi == NULL) {
		printf("Error: Can not allocate memory for D_multi. Aborting ... ...\n");
		return;
	}

	/* Compute C=alpha*op(A)*B+beta*C, A is an m-by-k sparse matrix in csc formate. 
	 * MKL: The csc format is specified by four arrays: values, columns, pointerB, pointerE;
	 * values: A real or complex array that contains the nonzero elements of the sparse matrix;
	 * rows: Element i of the array rows is the number of the row in the sparse matrix that
	 * contains the i-th value in the values array;
	 * pointerB: Element j of this integer array gives the index of the element in the values array 
	 * that is first nonzero element in a column j of A;
	 * pointerE: An integer array that contains column indices, such that pointerE[j] - pointerB[0] 
	 * is the index of the element in the values array that is last nonzero element in a column j of A. */
	const char T = 't';
	const int m = size_subsample[0];		// number of rows of the matrix A
	const int n = size_multi[2];				// number of columns of the matrix C
	const int k = size_subsample[1];		// number of columns of the matrix A
	const double alpha = 1.0;
	const char matdescra[6] = { 'G', 'L', 'N', 'F','*','*' };	// Specifies properties of the matrix used for operation;
																						// Only first four array elements are used;
																						// General; One-based
	int *indx = (int *)mkl_malloc(subsample_nzmax * sizeof(int), 32);
	int *pntrb = (int *)mkl_malloc(size_subsample[1] * sizeof(int), 32);
	int *pntre = (int *)mkl_malloc(size_subsample[1] * sizeof(int), 32);
	if (indx == NULL || pntrb == NULL || pntre == NULL) {
		printf("Error: can not allocate memory for indx, pntrb, or pntre. Aborting ......\n");
		return;
	}
	for (i = 0; i < subsample_nzmax; i++) {
		indx[i] = subsample_ir[i];
	}
	for (i = 0; i < size_subsample[1]; i++) {
		pntrb[i] = subsample_jc[i];
		pntre[i] = subsample_jc[i + 1];
	}
	
	const int ldb = size_subsample[0];		// Specifies the leading dimension of B for one-based indexing, 
																	// and the second dimension of B for zero-based indexing. 
	const double beta = 0.0;
	const int ldc = size_subsample[1];		// Spedifies the leading dimension of C for one-based indexing,
																	// and the second dimension of C for zero-based indexing. 
	mkl_dcscmm(&T, &m, &n, &k, &alpha, matdescra, (const double *)subsample_pr, (const int *)indx,
		(const int *)pntrb, (const int *)pntre, (const double *)multi, &ldb, &beta, D_multi, &ldc);
	/**************************************************************************** 
	 * CAUTION: ?cscmm interface for 0-based csc matrix supposes that the dense matrices are stored in 
	 * row-major order (C style); ?cscmm interface for 1-based csc matrix supposes that the matrices are
	 * stored in column-major order (Fortran style) 
	 ****************************************************************************/

	mkl_free(indx);
	mkl_free(pntrb);
	mkl_free(pntre);


	pthread_t *id = (pthread_t *)mkl_malloc(size_multi[0]/128 * size_multi[1]/128 * sizeof(pthread_t), 8);
	if (id == NULL) {
		printf("Error: can not allocate memory for id!\n");
		return;
	}
	struct thread_arg *arg= (struct thread_arg *)mkl_malloc(size_multi[0]/128 * size_multi[1]/128 * sizeof(struct thread_arg), 8);
	if (arg == NULL) {
		printf("Error: can not allocate memory for arg!\n");
		return;
	}

	for (i = 0; i < size_multi[0] / (128); i++) {		// 128 = 4*32
		for (j = 0; j < size_multi[1] / (128); j++) {	// 128 = 4*32
			(arg + i*size_multi[1]/128 + j)->D_multi = D_multi;
			(arg + i*size_multi[1]/128 + j)->D_hyper = D_hyper;
			(arg + i*size_multi[1]/128 + j)->multi = multi;
			(arg + i*size_multi[1]/128 + j)->hyper = estimated_hyper;
			(arg + i*size_multi[1]/128 + j)->i = i;
			(arg + i*size_multi[1]/128 + j)->j = j;
			(arg + i*size_multi[1]/128 + j)->size_multi[0] = size_multi[0];
			(arg + i*size_multi[1]/128 + j)->size_multi[1] = size_multi[1];

			pthread_create(id + (i*size_multi[1] / (128) + j), NULL, sub_thread, (void *)(arg + i*size_multi[1]/128 + j));			
		}
	}

	for (i = 0; i < size_multi[0] * size_multi[1] / 16384; i++) {	// 16384 = 128*128
		pthread_join(id[i], NULL);
	}

	mkl_free(id);
	mkl_free(arg);
	mkl_free(D_multi);
}

void *sub_thread(void *arg)
{
	struct thread_arg *p = (struct thread_arg *)arg;
	int k, ii, jj;

	double Y[4 * 4 * 31];
	double X[4 * 4 * 3];
	double A[3 * 31];
	//double YY[4 * 32 * 4 * 32 * 31] = { 0 };	TO LARGE!!! Run out of stack space.	
	double XX[4 * 32 * 4 * 32 * 3];		
	double *YY = (double *)mkl_malloc(4 * 32 * 4 * 32 * 31 * sizeof(double), 32);	// 507904 = 4*32*4*32*31
	//double *XX = (double *)mkl_malloc(4 * 32 * 4 * 32 * 3 * sizeof(double), 32);		
	if (YY == NULL || XX == NULL) {
		printf("Error: can not allocate memory for XX or YY!\n");
		return NULL;
	}

	double mean_multi[3];
	double mean_hyper[31];
	double sum;

	for (k = 0; k < 3; k++) {
		for (jj = 0; jj < 4; jj++) {
			//for (int ii = 0; ii < 4; ii++) {
			//	X[ii + jj * 4 + k * 16] = p->D_multi[p->i * 4 + ii + p->j * p->size_multi[0] / 32 * 4 + jj * p->size_multi[0] / 32 + k * p->size_multi[0] / 32 * p->size_multi[1] / 32];
			//}
			memcpy((void *)(X + (jj * 4 + k * 16)), 
				(void *)(p->D_multi + (p->i * 4 + p->j * p->size_multi[0] / 8 + jj * p->size_multi[0] / 32 + k * p->size_multi[0] / 32 * p->size_multi[1] / 32)), 
				4*sizeof(double));
		}
		// get mean value and subtract mean value.
		sum = 0;
		for (ii = 0; ii < 16; ii++) {
			sum += X[ii + k * 16];
		}
		mean_multi[k] = sum / 16.0;

		for (ii = 0; ii < 16; ii++) {
			X[ii + k * 16] -= mean_multi[k];
		}
	}

	for (k = 0; k < 31; k++) {
		for (jj = 0; jj < 4; jj++) {
			//for (int ii = 0; ii < 4; ii++) {
			//	Y[ii + jj * 4 + k * 16] = p->D_hyper[p->i * 4 + ii + p->j * p->size_multi[0] / 32 * 4 + jj * p->size_multi[0] / 32 + k * p->size_multi[0] / 32 * p->size_multi[1] / 32];
			//}
			memcpy((void *)(Y + (jj * 4 + k * 16)),
				(void *)(p->D_hyper + (p->i * 4 + p->j * p->size_multi[0] / 8 + jj * p->size_multi[0] / 32 + k * p->size_multi[0] / 32 * p->size_multi[1] / 32)),
				4 * sizeof(double));
		}
		// get mean value and substract mean value.
		sum = 0;
		for (ii = 0; ii < 16; ii++) {
			sum += Y[ii + k * 16];
		}
		mean_hyper[k] = sum / 16.0;

		for (ii = 0; ii < 16; ii++) {
			Y[ii + k * 16] -= mean_hyper[k];
		}
	}

	// least square estimation
	// A = (X^T*X)^-1*X^T*Y, where X is n-by-3, Y is n-by-31
	// \hat{Y} = X*A + epsilon
	lapack_int info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 16, 3, 31, X, 16, Y, 16);
	if (info == 0) {
		//printf("The execution is successful.\n");
	}
	else if (info < 0) {
		printf("illegal parameter!\n");
	}
	else {
		printf("Warning: the algprithm for computing the SVD failed to converge!\n");
	}
	/* Y is overwritten by the n-by-nrhs solution matrix A.
	* rows 1 to n of Y contain the least squares solution vectors. */

	for (k = 0; k < 3; k++) {
		for (jj = 0; jj < 128; jj++) {		// 128 = 4*32
			for (int ii = 0; ii < 4 * 32; ii++) {
				XX[ii + jj * 128 + k * 128 * 128] = 
					p->multi[p->i * 128 + ii + p->j * p->size_multi[0] * 128 + jj * p->size_multi[0] + k * p->size_multi[0] * p->size_multi[1]] - mean_multi[k];
			}
		}
	}

	// Compute hyperspectra using learned linear relation between hyperspectra and multispectra.
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16384, 31, 3, 1.0, XX, 16384, Y, 16, 0.0, YY, 16384);		// 16384 = 128*128

	// feed YY into estimated_hyper
	for (k = 0; k < 31; k++) {
		for (jj = 0; jj < 128; jj++) {		// 128 = 4*32
			for (int ii = 0; ii < 4 * 32; ii++) {
				p->hyper[p->i * 128 + ii + p->j * p->size_multi[0] * 128 + jj * p->size_multi[0] + k * p->size_multi[0] * p->size_multi[1]] = 
					YY[ii + jj * 128 + k * 128 * 128] + mean_hyper[k];
			}
		}
	}

	//mkl_free(XX);
	mkl_free(YY);

	return NULL;
}
