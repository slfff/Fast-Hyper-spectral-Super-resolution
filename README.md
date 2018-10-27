# Fast Hyper-spectral Image Super-resolution using Linear Regression

# C language implementation
  This version of C language implementation is developed on Linux. Intel Math Kernel Library (MKL) and Matlab R2014a extern library are
  internally used. Please make sure that these two dependencies are properly installed on your computer. Use Makefile to compile the 
  source code on your computer, but it is desired to modify the MKL and Matlab path according to your installation (our MKL is located 
  at /usr/opt, Matlab is located /usr/local.).
  
  The implementation of our method is contained in hss_lr.c, which explicitly requires low spatial resolution hyperspectral image 
  (LR-HSI), high spatial resolution multispectral image (HR-MSI), and degradation operator (D). The output is high spatial resolution 
  hyperspectral image (HR-HSI). In our demo (LR_C.c), we test our method on CAVE, Havard, and ICVL datasets respectively. 
  
# matlab implementation
  LR_latest*.m contains the implementation of our method for different patch size, i.e., 2*2, 4*4, 8*8, 16*16, 32*32. For proper run of 
  the code, please move the *.m into dictionary Dependency. 
