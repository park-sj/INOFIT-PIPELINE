#include "cuda_solver_internal.hpp"

// cuSOLVER handle and CUDA stream are initialized on construction
CUDASolverInternal::CUDASolverInternal() {
    cusolverStatus_t cusolver_status;
    cudaError_t cuda_error;
    cublasStatus_t cublas_status;
    cusparseStatus_t cusparse_status;
    
    cusolver_status = cusolverDnCreate(&cusolver_handle);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaStreamCreate(&cuda_stream);
    assert(cuda_error == cudaSuccess);

    cublas_status = cublasCreate(&cublas_handle);
    assert(cublas_status == CUBLAS_STATUS_SUCCESS);

    cusparse_status = cusparseCreate(&cusparse_handle);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
}

// cuSOLVER handle and CUDA stream are destoryed on descruction
CUDASolverInternal::~CUDASolverInternal() {
    cusolverStatus_t cusolver_status;
    cudaError_t cuda_error;
    cublasStatus_t cublas_status;
    cusparseStatus_t cusparse_status;

    cusolver_status = cusolverDnDestroy(cusolver_handle);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaStreamDestroy(cuda_stream);
    assert(cuda_error == cudaSuccess);

    cublas_status = cublasDestroy(cublas_handle);
    assert(cublas_status == CUBLAS_STATUS_SUCCESS);

    cusparse_status = cusparseDestroy(cusparse_handle);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
}

// Prints information of the i-th CUDA device
void CUDASolverInternal::printDeviceInfo(int i) {
    int device_count;
    cudaGetDeviceCount(&device_count);

    assert(i >= 0 && i < device_count);

    cudaDeviceProp device_prop;
    cudaGetDeviceProperties(&device_prop, i);

    std::cout << "Name: " << device_prop.name << std::endl;
    std::cout << "Compute capability: " << device_prop.major << "." << device_prop.minor
        << " (should be 1.3 or above to support double precision)" << std::endl;
    std::cout << "Clock rate: " << device_prop.clockRate << " kHz" << std::endl;
    std::cout << "Device copy overlap: ";
    if(device_prop.deviceOverlap) std::cout << "Enabled" << std::endl;
    else std::cout << "Disabled" << std::endl;
    std::cout << "Kernel execution timeout: ";
    if(device_prop.kernelExecTimeoutEnabled) std::cout << "Enabled" << std::endl;
    else std::cout << "Disabled" << std::endl;
    std::cout << "Total global mem: " << device_prop.totalGlobalMem << " B" << std::endl;
    std::cout << "Total constant mem: " << device_prop.totalConstMem << " B" << std::endl;
    std::cout << "Multiprocessor count: " << device_prop.multiProcessorCount << std::endl;
    std::cout << "Shared mem per mp: " << device_prop.sharedMemPerBlock << " B" << std::endl;
    std::cout << "Registers per mp: " << device_prop.regsPerMultiprocessor << std::endl;
    std::cout << "Threads in warp: " << device_prop.warpSize << std::endl;
    std::cout << "Max threads per block: " << device_prop.maxThreadsPerBlock << std::endl;
    std::cout << "Max thread dimensions: " << device_prop.maxThreadsDim[0] << " " << device_prop.maxThreadsDim[1] << " " << device_prop.maxThreadsDim[2] << std::endl;
    std::cout << "Max grid dimensions: " << device_prop.maxGridSize[0] << " " << device_prop.maxGridSize[1] << " " << device_prop.maxGridSize[2] << std::endl;

    return;
}

// Solves Ax = b via Cholesky decomposition
// Assumes that A is an n x n positive definite matrix
// Implemented with cuSOLVER library
void CUDASolverInternal::solveCholesky(int n, const double *A, const double *b, double *x) {
    cusolverStatus_t cusolver_status;
    cudaError_t cuda_error;

    // Assuming that the leading dimension is equal to n
    int lda = n;

    // Preparing device memories
    double *d_A = NULL;
    double *d_b = NULL;
    double *d_x = NULL;

    cuda_error = cudaMalloc((void **)&d_A, sizeof(double)*lda*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_b, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);

    // Copying A and b from host to device
    cuda_error = cudaMemcpy(d_A, A, sizeof(double)*lda*n, cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemcpy(d_b, b, sizeof(double)*n, cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);

    // Preparing cuSOLVER-related memories
    int buffer_size = 0;
    double *d_buffer = NULL;
    int *d_info = NULL;
    int info = 0;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    
    // Calculating the buffer size needed
    cusolver_status = cusolverDnDpotrf_bufferSize(cusolver_handle, uplo, n, (double*)A, n, &buffer_size);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaMalloc(&d_info, sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc(&d_buffer, sizeof(double) * buffer_size);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemset(d_info, 0, sizeof(int));
    assert(cuda_error == cudaSuccess);

    // Computes Cholesky decomposition
    cusolver_status = cusolverDnDpotrf(cusolver_handle, uplo, n, d_A, lda, d_buffer, buffer_size, d_info);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cuda_error == cudaSuccess);
    assert(info == 0);

    cuda_error = cudaMemcpy(d_x, d_b, sizeof(double)*n, cudaMemcpyDeviceToDevice);
    assert(cuda_error == cudaSuccess);

    // Solving Ax = b
    cusolver_status = cusolverDnDpotrs(cusolver_handle, uplo, n, 1, d_A, lda, d_x, n, d_info);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cuda_error = cudaDeviceSynchronize();
    assert(cuda_error == cudaSuccess);

    // Copying the result x from device to host
    cuda_error = cudaMemcpy(x, d_x, sizeof(double)*n, cudaMemcpyDeviceToHost);
    assert(cuda_error == cudaSuccess);

    // Freeing the device memories
    cuda_error = cudaFree(d_info);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_buffer);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_A);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_b);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_x);
    assert(cuda_error == cudaSuccess);
    return;
}

// Solves Ax = b via LU decomposition
// Assumes that A is a square matrix
// Implemented with cuSOLVER library
void CUDASolverInternal::solveLU(int n, const double *A, const double *b, double *x) {
    cusolverStatus_t cusolver_status;
    cudaError_t cuda_error;

    // Assuming that the leading dimension is equal to n
    int lda = n;

    // Preparing device memories
    double *d_A = NULL;
    double *d_b = NULL;
    double *d_x = NULL;

    cuda_error = cudaMalloc((void **)&d_A, sizeof(double)*lda*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_b, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);

    // Copying A and b from host to device
    cuda_error = cudaMemcpy(d_A, A, sizeof(double)*lda*n, cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemcpy(d_b, b, sizeof(double)*n, cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);

    // Preparing cuSOLVER-related memories
    int buffer_size = 0;
    double *d_buffer = NULL;
    int *d_info = NULL;
    int info = 0;
    int *d_ipiv = NULL;  // Pivoting sequence
    
    // Calculating the buffer size needed
    cusolver_status = cusolverDnDgetrf_bufferSize(cusolver_handle, n, n, (double*)A, lda, &buffer_size);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaMalloc(&d_info, sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc(&d_buffer, sizeof(double) * buffer_size);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc(&d_ipiv, sizeof(int)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemset(d_info, 0, sizeof(int));
    assert(cuda_error == cudaSuccess);

    // Computes LU decomposition
    cusolver_status = cusolverDnDgetrf(cusolver_handle, n, n, d_A, lda, d_buffer, d_ipiv, d_info);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_error = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cuda_error == cudaSuccess);
    assert(info == 0);

    cuda_error = cudaMemcpy(d_x, d_b, sizeof(double)*n, cudaMemcpyDeviceToDevice);
    assert(cuda_error == cudaSuccess);

    // Solving Ax = b
    cusolver_status = cusolverDnDgetrs(cusolver_handle, CUBLAS_OP_N, n, 1, d_A, lda, d_ipiv, d_x, n, d_info);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cuda_error = cudaDeviceSynchronize();
    assert(cuda_error == cudaSuccess);

    // Copying the result x from device to host
    cuda_error = cudaMemcpy(x, d_x, sizeof(double)*n, cudaMemcpyDeviceToHost);
    assert(cuda_error == cudaSuccess);

    // Freeing the device memories
    cuda_error = cudaFree(d_info);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_buffer);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_A);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_b);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaFree(d_x);
    assert(cuda_error == cudaSuccess);
    return;
}

// Solves Ax = b via conjugate gradient
// Assumes that A is a n x n symmetric positive definite sparse matrix
// Implemented with cuSPARSE library
void CUDASolverInternal::solveConjugateGradient(int n, int nz, const int *I, const int *J, const double *val, const double *rhs, double *x) {
    cusparseStatus_t cusparse_status;
    cublasStatus_t cublas_status;
    cudaError_t cuda_error;
    
    cusparseMatDescr_t descr = 0;
    cusparse_status = cusparseCreateMatDescr(&descr);
    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    int *d_col = NULL;
    int *d_row = NULL;
    double *d_val = NULL;
    double *d_x = NULL;
    double *d_r = NULL;
    double *d_p = NULL;
    double *d_Ax = NULL;

    cuda_error = cudaMalloc((void **)&d_col, sizeof(int)*nz);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_row, sizeof(int)*(n+1));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_val, sizeof(double)*nz);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_r, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_p, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_Ax, sizeof(double)*n);
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, I, (n+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, n*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, rhs, n*sizeof(double), cudaMemcpyHostToDevice);

    double alpha = 1.0, alpham1 = -1.0, beta = 0.0, r0 = 0.0, r1;

    cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasDaxpy(cublas_handle, n, &alpham1, d_Ax, 1, d_r, 1);
    cublas_status = cublasDdot(cublas_handle, n, d_r, 1, d_r, 1, &r1);
    assert(cublas_status == CUBLAS_STATUS_SUCCESS);

    int k = 1;
    const double tol = 0.00001;
    const int max_iter = 10000;
    double a, b, na, dot;

    while (r1 > tol*tol && k <= max_iter) {
        if (k > 1) {
            b = r1 / r0;
            cublas_status = cublasDscal(cublas_handle, n, &b, d_p, 1);
            assert(cublas_status == CUBLAS_STATUS_SUCCESS);
            cublas_status = cublasDaxpy(cublas_handle, n, &alpha, d_r, 1, d_p, 1);
            assert(cublas_status == CUBLAS_STATUS_SUCCESS);
        } else {
            cublas_status = cublasDcopy(cublas_handle, n, d_r, 1, d_p, 1);
            assert(cublas_status == CUBLAS_STATUS_SUCCESS);
        }

        cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublas_status = cublasDdot(cublas_handle, n, d_p, 1, d_Ax, 1, &dot);
        assert(cublas_status == CUBLAS_STATUS_SUCCESS);
        a = r1 / dot;

        cublas_status = cublasDaxpy(cublas_handle, n, &a, d_p, 1, d_x, 1);
        assert(cublas_status == CUBLAS_STATUS_SUCCESS);
        na = -a;
        cublas_status = cublasDaxpy(cublas_handle, n, &na, d_Ax, 1, d_r, 1);
        assert(cublas_status == CUBLAS_STATUS_SUCCESS);

        r0 = r1;
        cublas_status = cublasDdot(cublas_handle, n, d_r, 1, d_r, 1, &r1);
        assert(cublas_status == CUBLAS_STATUS_SUCCESS);
        cudaDeviceSynchronize();
        //printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }

    cudaMemcpy(x, d_x, n*sizeof(double), cudaMemcpyDeviceToHost);

    double rsum, diff, err = 0.0;

    for (int i = 0; i < n; i++) {
        rsum = 0.0;

        for (int j = I[i]; j < I[i+1]; j++) {
            rsum += val[j]*x[J[j]];
        }

        diff = fabs(rsum - rhs[i]);

        if (diff > err) {
            err = diff;
        }
    }

    //printf("Test Summary:  Error amount = %lf\n", err);

    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ax);

    return;
}