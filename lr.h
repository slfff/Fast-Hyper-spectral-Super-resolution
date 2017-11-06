#ifndef LR_H
#define LR_H

#include <mat.h>
#include <mkl.h>

double rmse(const mxArray *mx_truth, const double *estimated, const int m, const int n, const int k);
void hss_lr(const mxArray *mx_D_hyper, const mxArray *mx_multi, const mxArray *mx_subsample, double *estimated_hyper);

#endif
