#include <fstream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "mesh_alignment_3rd_party.hpp"
#include "tools.hpp"

// Should be float to work with OpenGR
typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh source_mesh, target_mesh, output_mesh;
    MeshAlignment3rdParty ma;

    if (argc != 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source file] [target off file]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading source mesh from " << std::string(argv[1]) << std::endl;
    std::ifstream source_ifstream(argv[1]);
    source_ifstream >> source_mesh;
    source_ifstream.close();

    if (source_mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading target mesh from " << std::string(argv[2]) << std::endl;
    std::ifstream target_ifstream(argv[2]);
    target_ifstream >> target_mesh;
    target_ifstream.close();

    if (target_mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[2]) << std::endl;
        return EXIT_FAILURE;
    }

    int len;
    std::string output_file;

    std::cout << "Aligning coarsely" << std::endl;
    output_mesh = source_mesh;
    ma.coarseAlignment(output_mesh, target_mesh);

    len = std::string(argv[2]).length();
    output_file = std::string(argv[2]).substr(0, len-4) + "_coarsely_aligned.off";
    std::cout << "Writing " << output_file << std::endl;
    std::ofstream coarse_ofstream(output_file);
    coarse_ofstream << output_mesh;
    coarse_ofstream.close();

    std::cout << "Aligning finely" << std::endl;
    output_mesh = source_mesh;
    ma.fineAlignment(output_mesh, target_mesh);

    len = std::string(argv[2]).length();
    output_file = std::string(argv[2]).substr(0, len-4) + "_finely_aligned.off";
    std::cout << "Writing " << output_file << std::endl;
    std::ofstream fine_ofstream(output_file);
    fine_ofstream << output_mesh;
    fine_ofstream.close();

    std::cout << "Aligning generally (coarsely first, then finely)" << std::endl;
    output_mesh = source_mesh;
    ma.generalAlignment(output_mesh, target_mesh);

    len = std::string(argv[2]).length();
    output_file = std::string(argv[2]).substr(0, len-4) + "_generally_aligned.off";
    std::cout << "Writing " << output_file << std::endl;
    std::ofstream general_ofstream(output_file);
    general_ofstream << output_mesh;
    general_ofstream.close();

    return EXIT_SUCCESS;
}