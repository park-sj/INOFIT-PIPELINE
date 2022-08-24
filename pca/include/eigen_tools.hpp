#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>

#ifndef EIGEN_TOOLS_HPP
#define EIGEN_TOOLS_HPP

class EigenTools {
public:
    EigenTools();

    bool hasNaN(Eigen::MatrixXd mat);
    bool hasInf(Eigen::MatrixXd mat);

    Eigen::MatrixXd readMatrixFromTxt(std::string file_name);
    void saveMatrixAsTxt(Eigen::MatrixXd mat, std::string file_name);
};

#endif