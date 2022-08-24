#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "cuda_solver.hpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    SimpleClock sc;
    int r;
    int c;

    if (argc != 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [number of rows] [number of columns]" << std::endl;
        return EXIT_FAILURE;
    }

    r = std::stoi(std::string(argv[1]));
    c = std::stoi(std::string(argv[2]));

    CUDASolver solver;
    std::cout << "Current device information" << std::endl;
    solver.printDeviceInfo();

    std::cout << "** Using double precision **" << std::endl;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(r, c);
    Eigen::VectorXd b = Eigen::VectorXd::Random(r);

    {
        std::cout << "(1) Solving Ax=b via Eigen bdcSvd()" << std::endl;
        sc.start();
        Eigen::VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(2) Solving Ax=b via Eigen fullPivHouseholderQr()" << std::endl;
        sc.start();
        Eigen::VectorXd x = A.fullPivHouseholderQr().solve(b);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(3) Solving Ax=b via converting to A^T*Ax=A^T*b with Eigen fullPivLu()" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = ATA.fullPivLu().solve(ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(4) Solving Ax=b via converting to A^T*Ax=A^T*b with Eigen llt() (Cholesky decomposition)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = ATA.llt().solve(ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(5) Solving Ax=b via converting to A^T*Ax=A^T*b with cuSOLVER (Cholesky decomposition)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = solver.solveCholesky(ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(5\') Solving Ax=b via converting to A^T*Ax=A^T*b with cuSOLVER (Cholesky decomposition) (repeated)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = solver.solveCholesky(ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(6) Solving Ax=b via converting to A^T*Ax=A^T*b with cuSOLVER (LU decomposition)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = solver.solveLU(ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(6\') Solving Ax=b via converting to A^T*Ax=A^T*b with cuSOLVER (LU decomposition) (repeated)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::VectorXd x = solver.solveLU(ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(7) Solving Ax=b via converting to A^T*Ax=A^T*b with cuSPARSE (conjugate gradient)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::SparseMatrix<double> sparse_ATA = ATA.sparseView();
        Eigen::VectorXd x = solver.solveConjugateGradient(sparse_ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    {
        std::cout << "(7\') Solving Ax=b via converting to A^T*Ax=A^T*b with cuSPARSE (conjugate gradient) (repeated)" << std::endl;
        sc.start();
        Eigen::MatrixXd ATA = A.transpose() * A;
        Eigen::VectorXd ATb = A.transpose() * b;
        Eigen::SparseMatrix<double> sparse_ATA = ATA.sparseView();
        Eigen::VectorXd x = solver.solveConjugateGradient(sparse_ATA, ATb);
        sc.stopAndPrint("Time");
        //std::cout << x << std::endl;
        std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
    }

    return 0;
}