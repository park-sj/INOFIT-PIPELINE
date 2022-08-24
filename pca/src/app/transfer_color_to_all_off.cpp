#include <iostream>
#include <string>
#include <filesystem>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "region_tools.hpp"
#include "tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source directory] [target directory] [color reference off file]" << std::endl;
        return EXIT_FAILURE;
    }

    RegionTools rt;
    FileNameEditor fe;

    std::string source_dir = std::string(argv[1]);
    std::string target_dir = std::string(argv[2]);
    std::string color_reference_file = std::string(argv[3]);

    Mesh color_reference_mesh;
    std::ifstream color_mesh_ifstream(color_reference_file);
    color_mesh_ifstream >> color_reference_mesh;
    color_mesh_ifstream.close();

    if(color_reference_mesh.is_empty())
    {
        std::cerr << "Cannot open the color reference off file " + color_reference_file << std::endl;
        return EXIT_FAILURE;
    }

    int v_n = color_reference_mesh.number_of_vertices();
    int e_n = color_reference_mesh.number_of_edges();
    int f_n = color_reference_mesh.number_of_faces();

    for(const auto & entry : std::filesystem::recursive_directory_iterator(source_dir))
    {
        if(entry.is_directory()) continue;

        if(entry.path().extension().compare(".off") == 0)
        {
            Mesh mesh;
            std::string input_file = entry.path().string();
            std::ifstream mesh_ifstream(input_file);
            mesh_ifstream >> mesh;
            mesh_ifstream.close();

            if(v_n != mesh.number_of_vertices() || e_n != mesh.number_of_edges() || f_n != mesh.number_of_faces())
            {
                std::cout << "Topology is different: " << input_file << std::endl;
                continue;
            }

            int len = input_file.length();
            std::string output_file_name_only = fe.extractName(input_file);
            std::string output_file = fe.addPath(target_dir, output_file_name_only) + ".off";

            rt.applyColor(mesh, color_reference_mesh);
            rt.writeCorrectOffFileWithFaceColor(mesh, output_file);

            std::cout << "Color transfered to " + output_file << std::endl;
        }
    }

    return EXIT_SUCCESS;
}