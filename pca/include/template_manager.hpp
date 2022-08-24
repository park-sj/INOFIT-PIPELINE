#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "mesh_tools.hpp"
#include "region_tools.hpp"

#ifndef TEMPLATE_MANAGER_HPP
#define TEMPLATE_MANAGER_HPP

class TemplateManager {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    TemplateManager();
    ~TemplateManager();

    void addMeshes(std::vector<std::string>& mesh_files, Direction face, Direction height);
    void addMeshes(std::vector<Mesh*>& meshes, Direction face, Direction height);
    void addMesh(std::string mesh_file, Direction face, Direction height);
    void addMesh(Mesh& mesh, Direction face, Direction height);

    int arrangeRegions(std::set<int>& regions_to_remove, std::vector<std::set<int> >& regions_to_merge, int additional_region_overlap = 0);
    void noseTipAlign();

    std::vector<Mesh*> getMeshes();
    std::vector<std::set<int> >* getRegionVertexMapping(int i);

    void clear();

private:
    std::vector<Mesh*> meshes;
    std::vector<Direction> initial_face_directions;
    std::vector<Direction> initial_height_directions;
    std::vector<Direction> current_face_directions;
    std::vector<Direction> current_height_directions;

    std::vector<std::vector<std::set<int> > > region_vertex_mappings;

};

#endif