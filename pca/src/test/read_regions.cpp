#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "region_tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh mesh;

    if (argc != 2) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [input off file]" << std::endl;
        return EXIT_FAILURE;
    }

    // Reading mesh
    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream mesh_ifstream(argv[1]);
    mesh_ifstream >> mesh;
    mesh_ifstream.close();

    // Computing region-vertex mapping
    RegionTools rt;
    std::vector<std::set<int> > region_vertex_mapping;
    rt.readRegionsFromColor(mesh, region_vertex_mapping);

    std::cout << "The number of regions is " << region_vertex_mapping.size() << std::endl;
    for(int i = 0; i < region_vertex_mapping.size(); i++) {
        std::cout << "Region " << i << " has " << region_vertex_mapping[i].size() << " vertices" << std::endl;
    }

    return EXIT_SUCCESS;
}