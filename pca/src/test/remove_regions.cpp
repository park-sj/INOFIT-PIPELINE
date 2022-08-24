#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "region_tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh mesh;

    if (argc <= 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [input off file] [region 0] [region 1] [region 2] ..." << std::endl;
        return EXIT_FAILURE;
    }

    // Reading mesh
    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream mesh_ifstream(argv[1]);
    mesh_ifstream >> mesh;
    mesh_ifstream.close();

    std::set<int> regions_to_remove;
    for(int i = 2; i < argc; i++) {
        regions_to_remove.insert(atoi(argv[i]));
    }

    // Computing region-vertex mapping
    RegionTools rt;
    std::vector<std::set<int> > region_vertex_mapping;
    rt.readRegionsFromColor(mesh, region_vertex_mapping);

    // Removing regions
    for(auto it = regions_to_remove.begin(); it != regions_to_remove.end(); ++it) {
        rt.removeFacesInRegion(mesh, region_vertex_mapping, *it);
    }
    rt.removeActually(mesh, region_vertex_mapping);

    // Writing mesh
    int len = std::string(argv[1]).length();
    std::string output_file = std::string(argv[1]).substr(0, len-4) + "_removed.off";
    std::cout << "Writing " << output_file << std::endl;
    rt.writeCorrectOffFileWithFaceColor(mesh, output_file);

    return EXIT_SUCCESS;
}