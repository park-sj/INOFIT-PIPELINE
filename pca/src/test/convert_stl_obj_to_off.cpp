#include <fstream>
#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

int main(int argc, char* argv[]) {
    Mesh mesh;

    if (argc != 2) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source stl/obj file]" << std::endl;
        return EXIT_FAILURE;
    }

    if(!OpenMesh::IO::read_mesh(mesh, argv[1])) {
        std::cerr << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    int len = std::string(argv[1]).length();
    std::string output_file = std::string(argv[1]).substr(0, len-4) + ".off";

    OpenMesh::IO::write_mesh(mesh, output_file);

    return EXIT_SUCCESS;
}