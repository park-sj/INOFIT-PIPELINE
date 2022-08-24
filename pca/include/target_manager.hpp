#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "mesh_tools.hpp"
#include "region_tools.hpp"

#ifndef TARGET_MANAGER_HPP
#define TARGET_MANAGER_HPP

class TargetManager {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    TargetManager();
    ~TargetManager();

    void addMeshes(std::vector<std::string>& mesh_files, Direction face, Direction height);
    void addMeshes(std::vector<Mesh*>& meshes, Direction face, Direction height);
    void addMesh(std::string mesh_file, Direction face, Direction height);
    void addMesh(Mesh& mesh, Direction face, Direction height);

    void noseTipAlign();
    void alignBackToTargetMesh(Mesh& mesh, int target_mesh_index);
    void saveMesh(Mesh& mesh, std::string file_name, bool do_remove_color);

    std::vector<Mesh*> getMeshes();

    void clear();

private:
    std::vector<Mesh*> meshes;
    std::vector<Direction> initial_face_directions;
    std::vector<Direction> initial_height_directions;
    std::vector<Direction> current_face_directions;
    std::vector<Direction> current_height_directions;

    Eigen::MatrixXd translations;

};

#endif