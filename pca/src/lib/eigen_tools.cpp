#include "eigen_tools.hpp"

EigenTools::EigenTools() {}

// Returns true if any of coefficients in mat is NaN
bool EigenTools::hasNaN(Eigen::MatrixXd mat) {
    bool ret = false;
    for(int c = 0; c < mat.cols(); c++) {
        for(int r = 0; r < mat.rows(); r++) {
            if(std::isnan(mat(r, c))) {
                ret = true;
                break;
            }
        }
    }
    return ret;
}

// Returns true if any of coefficients in mat is Inf or -Inf
bool EigenTools::hasInf(Eigen::MatrixXd mat) {
    bool ret = false;
    for(int c = 0; c < mat.cols(); c++) {
        for(int r = 0; r < mat.rows(); r++) {
            // Either positive infinity or negative infinity
            if(std::isinf(mat(r, c))) {
                ret = true;
                break;
            }
        }
    }
    return ret;
}

// Reads a matrix from a txt file
Eigen::MatrixXd EigenTools::readMatrixFromTxt(std::string file_name) {
    std::vector<std::list<std::string> > numbers;
    std::ifstream input_ifstream(file_name);
    std::string line;
    int row = 0;
    while(std::getline(input_ifstream, line)) {
        numbers.push_back(std::list<std::string>());
        std::stringstream ss(line);
        std::string number;
        while(std::getline(ss, number, ' ')) {
            numbers[row].push_back(number);
        }
        row += 1;
    }
    input_ifstream.close();
    
    int col = numbers[0].size();
    Eigen::MatrixXd mat(row, col);
    for(int r = 0; r < row; r++) {
        int c = 0;
        for(auto it = numbers[r].begin(); it != numbers[r].end(); ++it) {
            mat(r, c) = std::stod(*it);
            c += 1;
        }
    }

    return mat;
}

// Saves a matrix as a txt file
void EigenTools::saveMatrixAsTxt(Eigen::MatrixXd mat, std::string file_name) {
    Eigen::IOFormat txt_format(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n");

    std::ofstream txt_ofstream(file_name);
    if(txt_ofstream.is_open()) {
        txt_ofstream << mat.format(txt_format);
        txt_ofstream.close();
    } else {
        abort();
    }

    return;
}
