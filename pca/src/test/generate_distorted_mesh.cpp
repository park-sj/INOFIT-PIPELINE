#include <fstream>
#include <string>
#include <cstdlib>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "mesh_distortion.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh mesh;
    double distortion_rate;
    int N;

    if (argc != 4) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [input off file] [distortion rate] [number of output file]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream mesh_ifstream(argv[1]);
    mesh_ifstream >> mesh;
    mesh_ifstream.close();

    if (mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    // Tested distortion rate: 0.001
    distortion_rate = atof(argv[2]);
    N = atoi(argv[3]);

    // Generating a distorted mesh and writing to off file for N times
    MeshDistortion md;
    for(int i = 0; i < N; i++) {
        int len = std::string(argv[1]).length();
        std::string output_file = std::string(argv[1]).substr(0, len-4)
            + "_distorted_"
            + std::to_string(distortion_rate) + "_"
            + std::to_string(i) + ".off";

        std::cout << "Generating " << output_file << std::endl;
        Mesh output_mesh = mesh;
        md.distort(output_mesh, distortion_rate);

        std::cout << "Writing " << output_file << std::endl;
        std::ofstream mesh_ofstream(output_file);
        mesh_ofstream << output_mesh;
        mesh_ofstream.close();
    }
	
    return EXIT_SUCCESS;
}