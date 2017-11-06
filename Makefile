LR_C: LR_C.c hss_lr.c rmse.c
	gcc LR_C.c hss_lr.c rmse.c -o LR_C -I/usr/local/MATLAB/R2014a/extern/include -I/opt/intel/mkl/include -L/usr/local/MATLAB/R2014a/bin/glnxa64 -L/opt/intel/mkl/lib/intel64 -lpthread -lm -lmat -lmx -lmkl_core -lmkl_rt -Wl,-rpath,/usr/local/MATLAB/R2014a/bin/glnxa64 -O2 -std=c99
