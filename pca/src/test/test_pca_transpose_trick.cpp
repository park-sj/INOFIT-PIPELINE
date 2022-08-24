#include <iostream>
#include <string>
#include <cstdlib>
#include <Eigen/Dense>

int main(int argc, char* argv[]) {
    int row;
    int col;

    if (argc != 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [row size] [column size]" << std::endl;
        return EXIT_FAILURE;
    }

    row = atoi(argv[1]);
    col = atoi(argv[2]);

    // Generating row x col random matrix A
    std::cout << "Randomly generating a " << row << "x" << col << " matrix A (it will be centered)" << std::endl;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(row, col);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(row);
    for(int i = 0; i < col; i++) {
        sum += A.col(i);
    }
    Eigen::VectorXd mean;

    // Comparing the eigenvalues of A^T*A and A*A^T when mean is divided by (col)
    mean = sum / col;
    Eigen::MatrixXd A1 = A;
    A1.colwise() -= mean;

    Eigen::MatrixXd C11 = A1 * A1.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver11(C11);
    if(solver11.info() != Eigen::Success) abort();
    Eigen::MatrixXd S11 = solver11.eigenvalues();
    std::cout << "Eigenvalues of A*A^T (when mean is divided by (col)): " << std::endl;
    std::cout << S11 << std::endl;

    Eigen::MatrixXd C12 = A1.transpose() * A1;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver12(C12);
    if(solver12.info() != Eigen::Success) abort();
    Eigen::MatrixXd S12 = solver12.eigenvalues();
    std::cout << "Eigenvalues of A^T*A (when mean is divided by (col)): " << std::endl;
    std::cout << S12 << std::endl;

    // Comparing the eigenvalues of A^T*A and A*A^T when mean is divided by (col-1)
    mean = sum / (col-1);
    Eigen::MatrixXd A2 = A;
    A2.colwise() -= mean;

    Eigen::MatrixXd C21 = A2 * A2.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver21(C21);
    if(solver21.info() != Eigen::Success) abort();
    Eigen::MatrixXd S21 = solver21.eigenvalues();
    std::cout << "Eigenvalues of A*A^T (when mean is divided by (col-1)): " << std::endl;
    std::cout << S21 << std::endl;

    Eigen::MatrixXd C22 = A2.transpose() * A2;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver22(C22);
    if(solver22.info() != Eigen::Success) abort();
    Eigen::MatrixXd S22 = solver22.eigenvalues();
    std::cout << "Eigenvalues of A^T*A (when mean is divided by (col-1)): " << std::endl;
    std::cout << S22 << std::endl;

    return EXIT_SUCCESS;
}