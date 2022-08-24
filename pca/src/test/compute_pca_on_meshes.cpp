#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <Eigen/Dense>
#include "mesh_pca.hpp"
#include "mesh_distortion.hpp"
#include "tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh source_mesh, output_mesh;
    double distortion_rate;
    int n_m;
    int n_k;

    if (argc != 5) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source off file] [distortion rate] [number of distorted meshes] [number of principal components to be used in reconstruction]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream mesh_ifstream(argv[1]);
    mesh_ifstream >> source_mesh;
    output_mesh = source_mesh;
    mesh_ifstream.close();

    if (source_mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    // Tested distortion rate: 0.001
    distortion_rate = atof(argv[2]);
    n_m = atoi(argv[3]);
    n_k = atoi(argv[4]);

    if(n_k > n_m) {
        std::cout << "The number of principal components cannot exceed the number of meshes" << std::endl;
        std::cout << "Therefore, the number of principal components are set equal to the number of meshes" << std::endl;
        n_k = n_m;
    }

    // Generating distorted meshes based on the source mesh
    MeshDistortion md;
    std::vector<Mesh*> meshes;
    for(int i = 0; i < n_m; i++) {
        std::cout << "Generating distorted mesh " << i << std::endl;
        Mesh* mesh = new Mesh();
        *mesh = source_mesh;
        md.distort(*mesh, distortion_rate);
        meshes.push_back(mesh);
    }

    // Computing PCA
    std::cout << "Computing PCA from distorted meshes" << std::endl;
    MeshPCA pca;
    pca.inputMeshes(meshes);
    pca.computeEigenvectors();

    // Approximating the source mesh
    std::cout << "Approximating the source mesh based on computed PCA" << std::endl;
    auto projection = pca.approximate(output_mesh, n_k);

    // Reconstructing a mesh from the projection
    std::cout << "Reconstructing a mesh from the projection" << std::endl;
    pca.reconstructFrom(projection, output_mesh);

    // Generating the reconstructed mesh
    int len = std::string(argv[1]).length();
    std::string output_file = std::string(argv[1]).substr(0, len-4)
        + "_reconstructed_"
        + std::to_string(distortion_rate) + "_"
        + std::to_string(n_m) + "_"
        + std::to_string(n_k) + ".off";

    std::cout << "Writing to " << output_file << std::endl;
    std::ofstream mesh_ofstream(output_file);
    mesh_ofstream << output_mesh;
    mesh_ofstream.close();

    for(int i = 0; i < n_m; i++) {
        delete meshes[i];
    }

    return EXIT_SUCCESS;
}