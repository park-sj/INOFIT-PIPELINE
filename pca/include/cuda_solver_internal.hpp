#include <iostream>
#include <cassert>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cusparse.h>
#include <cublas_v2.h>

class CUDASolverInternal {
public:
    CUDASolverInternal();
    ~CUDASolverInternal();

    void printDeviceInfo(int i);

    void solveCholesky(int n, const double *A, const double *b, double *x);
    void solveLU(int n, const double *A, const double *b, double *x);

    void solveConjugateGradient(int n, int nz, const int *I, const int *J, const double *val, const double *b, double *x);

private:
    cusolverDnHandle_t cusolver_handle;
    cudaStream_t cuda_stream;
    cublasHandle_t cublas_handle;
    cusparseHandle_t cusparse_handle;

};