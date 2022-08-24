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
    int region;
    int hops;

    if (argc != 4) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [input off file] [region] [hops]" << std::endl;
        return EXIT_FAILURE;
    }

    // Reading mesh
    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream mesh_ifstream(argv[1]);
    mesh_ifstream >> mesh;
    mesh_ifstream.close();

    // Reading arguments
    region = atoi(argv[2]);
    hops = atoi(argv[3]);

    // Computing region-vertex mapping
    RegionTools rt;
    std::vector<std::set<int> > region_vertex_mapping;
    rt.readRegionsFromColor(mesh, region_vertex_mapping);

    // Setting the colors of all regions to white
    Color color = rt.getRegionColor(mesh, region_vertex_mapping, region);
    rt.clearColorToWhite(mesh);

    // Widening the region
    rt.widenRegion(mesh, region_vertex_mapping, region, hops);

    // Setting the color of widened region to its original color
    rt.setRegionColor(mesh, region_vertex_mapping, region, color);

    // Writing mesh
    int len = std::string(argv[1]).length();
    std::string output_file = std::string(argv[1]).substr(0, len-4) + "_widened.off";
    std::cout << "Writing " << output_file << std::endl;
    rt.writeCorrectOffFileWithFaceColor(mesh, output_file);

    return EXIT_SUCCESS;
}