#include <iostream>
#include <vector>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <Eigen/Dense>
#include "region_pca_manager.hpp"
#include "mesh_correspondence.hpp"
#include "mesh_tools.hpp"
#include "region_tools.hpp"
#include "eigen_tools.hpp"
#include "tools.hpp"

#ifndef DATA_DRIVEN_NONRIGID_ICP_HPP
#define DATA_DRIVEN_NONRIGID_ICP_HPP

class DataDrivenNonrigidICP {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    // To have a fixed-size vectorizable Eigen type as a member variable
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    DataDrivenNonrigidICP();
    void init(RegionPCAManager& rm, Mesh& target_mesh, std::vector<double>& region_weights);
    double iterateUntil(double threshold_error, int max_iteration, double threshold_distance, double per_iteration_smoothness_weight = 1.0, double per_iteration_stiffness_weight = 1.0);
    void printStatus(double threshold_distance, bool overall_only = false);
    Mesh getResultMesh();
    void getParameters(double& s, Eigen::Matrix3d& rot, Eigen::VectorXd& m, Eigen::Vector3d& t);

private:
    void computeVertexWeights(std::vector<double>& region_weights);
    void iterate(double threshold_distance, double per_iteration_smoothness_weight, double per_iteration_stiffness_weight);
    void computeCorrespondence();
    void update(double ds, Eigen::Vector3d dtheta, Eigen::VectorXd dm, Eigen::VectorXd dt);
    double computeMeanSquaredError();
    Eigen::Matrix3d getRotationForm(Eigen::Vector3d v);

    double s;
    Eigen::Matrix3d rot;    // Rotation is stored and accumulated as a matrix
    Eigen::VectorXd m;
    Eigen::Vector3d t;

    Eigen::VectorXd current_mesh_vector;
    Eigen::VectorXd correspondence_vector;
    RegionPCAManager *rm;

    MeshCorrespondence mc;

    std::vector<double> vertex_weights;

};

#endif