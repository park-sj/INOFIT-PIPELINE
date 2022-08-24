#include "cuda_solver.hpp"

CUDASolver::CUDASolver() {}

// Prints information of the first CUDA device
void CUDASolver::printDeviceInfo() {
    internal.printDeviceInfo(0);
    return;
}

// Solves Ax = b via Cholesky decomposition
// Assumes that A is an n x n positive definite matrix
Eigen::VectorXd CUDASolver::solveCholesky(Eigen::MatrixXd A, Eigen::VectorXd b) {
    assert(A.cols() == A.rows());
    Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
    internal.solveCholesky(b.size(), A.data(), b.data(), x.data());
    return x;
}

// Solves Ax = b via LU decomposition
// Assumes that A is a square matrix
Eigen::VectorXd CUDASolver::solveLU(Eigen::MatrixXd A, Eigen::VectorXd b) {
    assert(A.cols() == A.rows());
    Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
    internal.solveLU(b.size(), A.data(), b.data(), x.data());
    return x;
}

// Solves Ax = b via conjugate gradient
// Assumes that A is a n x n symmetric positive definite sparse matrix
Eigen::VectorXd CUDASolver::solveConjugateGradient(Eigen::SparseMatrix<double> A, Eigen::VectorXd b) {
    assert(A.cols() == A.rows());
    Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());

    A.makeCompressed();

    // Not checking if A is row major or col major since A is assumed to be symmetric
    internal.solveConjugateGradient(b.size(), A.nonZeros(), A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr(), b.data(), x.data());

    return x;
}