#include <cassert>
#include <cstdlib>
#include <vector>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <Eigen/Dense>

#ifndef MESH_TOOLS_HPP
#define MESH_TOOLS_HPP

class Direction {
public:
    Direction();
    Direction(std::string d);

    std::string getString();
    bool isOrthogonalTo(Direction direction);
    Direction oppositeDirection();
    Direction clockwiseDirection(Direction axis);
    Direction counterclockwiseDirection(Direction axis);
    Direction outerProductDirection(Direction direction_b);

private:
    void assertCorrectDirection();

    std::string d;
};

class MeshTools {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    MeshTools();

    Eigen::VectorXd meshToMeshVector(Mesh& mesh);
    void meshVectorToMesh(Eigen::VectorXd mesh_vector, Mesh& topology_mesh);
    Eigen::MatrixXd meshesToMeshMatrix(std::vector<Mesh*>& meshes);
    void meshMatrixToMeshes(Eigen::MatrixXd& mesh_matrix, std::vector<Mesh*>& topology_meshes);

    Eigen::VectorXd makeMeshMatrixCentered(Eigen::MatrixXd& mesh_matrix);
    std::vector<double> scaleMeshMatrixHeight(Eigen::MatrixXd& mesh_matrix, double h, Direction height_direction);
    void rotateMeshMatrixAxis(Eigen::MatrixXd& mesh_matrix, Direction face_from, Direction height_from, Direction face_to, Direction height_to);
    Eigen::MatrixXd translateMeshMatrixToOrigin(Eigen::MatrixXd& mesh_matrix);
    Eigen::MatrixXd translateMeshMatrixTipToOrigin(Eigen::MatrixXd& mesh_matrix, Direction d);
    void translateMeshMatrixEach(Eigen::MatrixXd& mesh_matrix, Eigen::MatrixXd translations);
    void translateMeshMatrixAll(Eigen::MatrixXd& mesh_matrix, Eigen::Vector3d translation);

    Eigen::VectorXd makeMeshesCentered(std::vector<Mesh*>& meshes);
    std::vector<double> scaleMeshesHeight(std::vector<Mesh*>& meshes, double h, Direction height_direction);
    void rotateMeshesAxis(std::vector<Mesh*>& meshes, Direction face_from, Direction height_from, Direction face_to, Direction height_to);
    Eigen::MatrixXd translateMeshesToOrigin(std::vector<Mesh*>& meshes);
    Eigen::MatrixXd translateMeshesTipToOrigin(std::vector<Mesh*>& meshes, Direction d);
    void translateMeshesEach(std::vector<Mesh*>& meshes, Eigen::MatrixXd translations);
    void translateMeshesAll(std::vector<Mesh*>& meshes, Eigen::Vector3d translation);

    double scaleMeshHeight(Mesh& mesh, double h, Direction height_direction);
    void rotateMeshAxis(Mesh& mesh, Direction face_from, Direction height_from, Direction face_to, Direction height_to);
    Eigen::Vector3d translateMeshToOrigin(Mesh& mesh);
    Eigen::Vector3d translateMeshTipToOrigin(Mesh& mesh, Direction d);
    void translateMesh(Mesh& mesh, Eigen::Vector3d translation);

private:
    void rotateMeshMatrixAxisOnce(Eigen::MatrixXd& mesh_matrix, Direction axis_from, Direction axis_to);
    void rotateMeshAxisOnce(Mesh& mesh, Direction axis_from, Direction axis_to);

    constexpr static const double INITIAL_MIN_VALUE = 987654321.0;
    constexpr static const double INITIAL_MAX_VALUE = -987654321.0;
};

#endif