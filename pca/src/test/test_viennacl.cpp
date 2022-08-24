#include <iostream>
#include <string>

#ifdef __APPLE__
    #include "TargetConditionals.h"
    #ifdef TARGET_OS_MAC
        #include <OpenCL/opencl.h>
    #endif
#else
    #include <CL/cl.h>
#endif

#include <viennacl/matrix.hpp>
#include <viennacl/linalg/lu.hpp>
#include <viennacl/linalg/cg.hpp>
#include <Eigen/Dense>
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

    std::cout << "Current device information:" << std::endl;
    std::cout << viennacl::ocl::current_device().info() << std::endl;

    if(viennacl::ocl::current_device().double_support()) {
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

        viennacl::matrix<double> vcl_ATA(c, c);
        viennacl::vector<double> vcl_ATb(c);
        viennacl::vector<double> vcl_x(c);

        {
            std::cout << "(5) Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL lu_factorize() and lu_substitute()" << std::endl;
            sc.start();
            Eigen::MatrixXd ATA = A.transpose() * A;
            Eigen::VectorXd ATb = A.transpose() * b;
            Eigen::VectorXd x = Eigen::VectorXd::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            viennacl::linalg::lu_factorize(vcl_ATA);
            viennacl::linalg::lu_substitute(vcl_ATA, vcl_ATb);
            viennacl::copy(vcl_ATb, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(5\') Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL lu_factorize() and lu_substitute() (repeated)" << std::endl;
            sc.start();
            Eigen::MatrixXd ATA = A.transpose() * A;
            Eigen::VectorXd ATb = A.transpose() * b;
            Eigen::VectorXd x = Eigen::VectorXd::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            viennacl::linalg::lu_factorize(vcl_ATA);
            viennacl::linalg::lu_substitute(vcl_ATA, vcl_ATb);
            viennacl::copy(vcl_ATb, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(6) Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL solve(cg_tag())" << std::endl;
            sc.start();
            Eigen::MatrixXd ATA = A.transpose() * A;
            Eigen::VectorXd ATb = A.transpose() * b;
            Eigen::VectorXd x = Eigen::VectorXd::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            vcl_x = viennacl::linalg::solve(vcl_ATA, vcl_ATb, viennacl::linalg::cg_tag());
            viennacl::copy(vcl_x, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(6\') Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL solve(cg_tag()) (repeated)" << std::endl;
            sc.start();
            Eigen::MatrixXd ATA = A.transpose() * A;
            Eigen::VectorXd ATb = A.transpose() * b;
            Eigen::VectorXd x = Eigen::VectorXd::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            vcl_x = viennacl::linalg::solve(vcl_ATA, vcl_ATb, viennacl::linalg::cg_tag());
            viennacl::copy(vcl_x, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }
    } else {
        std::cout << "** Using single precision **" << std::endl;

        Eigen::MatrixXf A = Eigen::MatrixXf::Random(r, c);
        Eigen::VectorXf b = Eigen::VectorXf::Random(r);

        {
            std::cout << "(1) Solving Ax=b via Eigen bdcSvd()" << std::endl;
            sc.start();
            Eigen::VectorXf x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(2) Solving Ax=b via Eigen fullPivHouseholderQr()" << std::endl;
            sc.start();
            Eigen::VectorXf x = A.fullPivHouseholderQr().solve(b);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(3) Solving Ax=b via converting to A^T*Ax=A^T*b with Eigen fullPivLu()" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = ATA.fullPivLu().solve(ATb);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(4) Solving Ax=b via converting to A^T*Ax=A^T*b with Eigen llt() (Cholesky decomposition)" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = ATA.llt().solve(ATb);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        viennacl::matrix<float> vcl_ATA(c, c);
        viennacl::vector<float> vcl_ATb(c);
        viennacl::vector<float> vcl_x(c);

        {
            std::cout << "(5) Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL lu_factorize() and lu_substitute()" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = Eigen::VectorXf::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            viennacl::linalg::lu_factorize(vcl_ATA);
            viennacl::linalg::lu_substitute(vcl_ATA, vcl_ATb);
            viennacl::copy(vcl_ATb, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(5\') Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL lu_factorize() and lu_substitute() (repeated)" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = Eigen::VectorXf::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            viennacl::linalg::lu_factorize(vcl_ATA);
            viennacl::linalg::lu_substitute(vcl_ATA, vcl_ATb);
            viennacl::copy(vcl_ATb, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(6) Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL solve(cg_tag())" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = Eigen::VectorXf::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            vcl_x = viennacl::linalg::solve(vcl_ATA, vcl_ATb, viennacl::linalg::cg_tag());
            viennacl::copy(vcl_x, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }

        {
            std::cout << "(6\') Solving Ax=b via converting to A^T*Ax=A^T*b with ViennaCL solve(cg_tag()) (repeated)" << std::endl;
            sc.start();
            Eigen::MatrixXf ATA = A.transpose() * A;
            Eigen::VectorXf ATb = A.transpose() * b;
            Eigen::VectorXf x = Eigen::VectorXf::Zero(c);
            viennacl::copy(ATA, vcl_ATA);
            viennacl::copy(ATb, vcl_ATb);
            vcl_x = viennacl::linalg::solve(vcl_ATA, vcl_ATb, viennacl::linalg::cg_tag());
            viennacl::copy(vcl_x, x);
            sc.stopAndPrint("Time");
            //std::cout << x << std::endl;
            std::cout << "Residual (l2 norm): " << (b - A*x).norm() << std::endl;
        }
    }

    return 0;
}