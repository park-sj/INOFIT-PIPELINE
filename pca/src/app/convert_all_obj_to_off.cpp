#include <fstream>
#include <string>
#include <filesystem>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "tools.hpp"

typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source directory] [target directory]" << std::endl;
        return EXIT_FAILURE;
    }

    FileNameEditor fe;

    std::string source_dir = std::string(argv[1]);
    std::string target_dir = std::string(argv[2]);

    for(const auto & entry : std::filesystem::recursive_directory_iterator(source_dir))
    {
        if(entry.is_directory()) continue;

        if(entry.path().extension().compare(".obj") == 0)
        {
            Mesh mesh;
            std::string input_file = entry.path().string();
            if(!OpenMesh::IO::read_mesh(mesh, input_file)) {
                std::cerr << "Cannot read " << input_file << std::endl;
            } else {
                int len = input_file.length();
                std::string output_file_name_only = fe.extractName(input_file);
                std::string output_file = fe.addPath(target_dir, output_file_name_only) + ".off";

                OpenMesh::IO::write_mesh(mesh, output_file);

                std::cout << output_file << " generated" << std::endl;
            }
        }
    }

    return EXIT_SUCCESS;
}