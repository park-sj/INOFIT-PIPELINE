#include <vector>
#include <set>
#include <list>
#include <cassert>
#include <string>
#include <algorithm>
#include <cmath>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <Eigen/Dense>
#include "mesh_pca.hpp"
#include "mesh_tools.hpp"
#include "region_tools.hpp"
#include "eigen_tools.hpp"

#ifndef REGION_PCA_MANAGER
#define REGION_PCA_MANAGER

class RegionPCAManager {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    RegionPCAManager();

    void construct(std::vector<Mesh*>& meshes, std::vector<std::set<int> >& region_vertex_mapping, int n_k, double smoothness_factor, std::vector<double>& stiffness_factor);
    
    Eigen::MatrixXd getMBlock(int v);
    Eigen::MatrixXd getSmoothnessRegularization();
    Eigen::MatrixXd getStiffnessRegularization();
    std::vector<std::set<int> >* getRegionVertexMapping();
    Eigen::VectorXd getMeanVector();
    Mesh getColoredMeanMesh();
    int getNkNr();

    void save(std::string f_M, std::string f_smoothness, std::string f_stiffness, std::string f_rvm, std::string f_cmm);
    void read(std::string f_M, std::string f_smoothness, std::string f_stiffness, std::string f_rvm, std::string f_cmm);

private:
    void saveMAsTxt(std::string file_name);
    void saveSmoothnessAsTxt(std::string file_name);
    void saveStiffnessAsTxt(std::string file_name);
    void saveRegionVertexMappingAsTxt(std::string file_name);
    void saveMeanAsOffWithColor(std::string file_name);

    void readMFromTxt(std::string file_name);
    void readSmoothnessFromTxt(std::string file_name);
    void readStiffnessFromTxt(std::string file_name);
    void readRegionVertexMappingFromTxt(std::string file_name);
    void readColoredMeanFromOff(std::string file_name);

    Eigen::MatrixXd M;
    Eigen::MatrixXd smoothness;
    Eigen::MatrixXd stiffness;
    Eigen::VectorXd mean_vector;

    std::vector<std::set<int> > region_vertex_mapping;

    Mesh colored_mean_mesh;
};

#endif