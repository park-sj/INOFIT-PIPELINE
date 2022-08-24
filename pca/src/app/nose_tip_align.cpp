#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "mesh_tools.hpp"
#include "region_tools.hpp"
#include "tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source directory] [target directory] [face direction]" << std::endl;
        return EXIT_FAILURE;
    }

    MeshTools mt;
    RegionTools rt;
    FileNameEditor fe;

    std::string source_dir = std::string(argv[1]);
    std::string target_dir = std::string(argv[2]);
    std::string face_direction = std::string(argv[3]);
    
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

            if(mesh.is_empty())
            {
                std::cerr << "Cannot open the off file " + input_file << std::endl;
                continue;
            }

            int len = input_file.length();
            std::string output_file_name_only = fe.extractName(input_file);
            std::string output_file = fe.addPath(target_dir, output_file_name_only) + "_nosetip.off";

            mt.translateMeshTipToOrigin(mesh, Direction(face_direction));
            rt.writeCorrectOffFileWithFaceColor(mesh, output_file);

            std::cout << output_file << " generated" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}