#include <iostream>
#include <string>
#include "parameter_reader.hpp"
#include "template_manager.hpp"
#include "region_pca_manager.hpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    ParameterReader params;
    TemplateManager tm;
    RegionPCAManager rm;
    SimpleClock sc;

    if (argc != 2) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [parameter json file]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading parameters from " << std::string(argv[1]) << std::endl;
    if(!params.readForPreprocessing(std::string(argv[1]))) {
        std::cout << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    sc.start();
    std::cout << "Reading template meshes" << std::endl;
    tm.addMeshes(params.template_files, params.template_face_direction, params.template_height_direction);
    std::cout << "Reading template meshes took " << sc.stop() << "s" << std::endl;

    sc.start();
    std::cout << "Arranging mesh regions" << std::endl;
    tm.arrangeRegions(params.regions_to_remove, params.regions_to_merge, params.additional_region_overlap);
    std::cout << "Arranging mesh regions took " << sc.stop() << "s" << std::endl;

    if(params.is_template_aligned) {
        std::cout << "Template meshes are already aligned -> Alignment skipped" << std::endl;
    } else {
        sc.start();
        std::cout << "Aligning meshes" << std::endl;
        tm.noseTipAlign();
        std::cout << "Aligning meshes took " << sc.stop() << "s" << std::endl;
    }

    auto aligned_template_meshes = tm.getMeshes();
    auto region_vertex_mapping = *(tm.getRegionVertexMapping(0));

    sc.start();
    std::cout << "Computing per-region PCA" << std::endl;
    rm.construct(aligned_template_meshes, region_vertex_mapping, params.n_k, params.smoothness_factor, params.stiffness_factors);
    std::cout << "Computing per-region PCA took " << sc.stop() << "s" << std::endl;

    FileNameEditor fe;
    std::string name = fe.extractName(std::string(argv[1]));
    std::string matrix_dir = params.matrix_file_directory;
    std::string region_dir = params.region_file_directory;
    std::string mesh_dir = params.mesh_file_directory;

    std::string f_M = fe.addPath(matrix_dir, name) + "_M.txt";
    std::string f_smoothness = fe.addPath(matrix_dir, name) + "_smoothness.txt";
    std::string f_stiffness = fe.addPath(matrix_dir, name) + "_stiffness.txt";
    std::string f_rvm = fe.addPath(region_dir, name) + "_rvm.txt";
    std::string f_cmm = fe.addPath(mesh_dir, name) + "_cmm.off";

    rm.save(f_M, f_smoothness, f_stiffness, f_rvm, f_cmm);

    std::cout << "Matrix M is saved to " << f_M << std::endl;
    std::cout << "Smoothness matrix is saved to " << f_smoothness << std::endl;
    std::cout << "Stiffness matrix is saved to " << f_stiffness << std::endl;
    std::cout << "Region-vertex mapping is saved to " << f_rvm << std::endl;
    std::cout << "Colored mean mesh is saved to " << f_cmm << std::endl;

    std::cout << "Preprocessing is successfully completed" << std::endl;

    return EXIT_SUCCESS;
}