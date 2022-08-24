#include <vector>
#include <cassert>
#include <cstdlib>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <Eigen/Dense>
#include "mesh_tools.hpp"

#ifndef MESH_PCA_HPP
#define MESH_PCA_HPP

class MeshPCA {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    MeshPCA();
    void inputMeshes(std::vector<Mesh*>& meshes);
    void inputMeshMatrix(Eigen::MatrixXd mesh_matrix);
    void computeEigenvectors();
    Eigen::VectorXd approximate(Mesh& mesh, int n_k);
    Eigen::VectorXd approximate(Eigen::VectorXd mesh_vector, int n_k);
    void reconstructFrom(Eigen::VectorXd projection, Mesh& topology_mesh);
    Eigen::VectorXd reconstructFrom(Eigen::VectorXd projection);
    
    Eigen::VectorXd getMeanVector();
    void getMeanMesh(Mesh& topology_mesh);
    Eigen::MatrixXd getEigenvectorBlock(int i, int j, int r, int c);
    Eigen::VectorXd getEigenvalues();

private:
    Eigen::MatrixXd B;
    Eigen::VectorXd mean_vector;
    Eigen::MatrixXd P;
    Eigen::MatrixXd S;

};

#endif