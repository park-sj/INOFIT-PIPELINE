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
typedef CGAL::Color Color;

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

    std::set<int> regions_to_merge;
    for(int i = 2; i < argc; i++) {
        regions_to_merge.insert(atoi(argv[i]));
    }

    // Computing region-vertex mapping
    RegionTools rt;
    std::vector<std::set<int> > region_vertex_mapping;
    rt.readRegionsFromColor(mesh, region_vertex_mapping);

    // Merging specified regions
    std::cout << "Merging regions" << std::endl;
    rt.mergeRegions(mesh, region_vertex_mapping, regions_to_merge);

    // Coloring the merged region for visualization
    Color new_region_color = rt.getRegionColor(mesh, region_vertex_mapping, *(regions_to_merge.begin()));
    rt.setRegionColor(mesh, region_vertex_mapping, *(regions_to_merge.begin()), new_region_color);

    // Writing mesh
    int len = std::string(argv[1]).length();
    std::string output_file = std::string(argv[1]).substr(0, len-4) + "_merged.off";
    std::cout << "Writing " << output_file << std::endl;
    rt.writeCorrectOffFileWithFaceColor(mesh, output_file);

    return EXIT_SUCCESS;
}