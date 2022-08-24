#include <cassert>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "cuda_solver_internal.hpp"

class CUDASolver {
public:
    CUDASolver();
    
    void printDeviceInfo();
    
    Eigen::VectorXd solveCholesky(Eigen::MatrixXd A, Eigen::VectorXd b);
    Eigen::VectorXd solveLU(Eigen::MatrixXd A, Eigen::VectorXd b);

    Eigen::VectorXd solveConjugateGradient(Eigen::SparseMatrix<double> A, Eigen::VectorXd b);

private:
    CUDASolverInternal internal;

};