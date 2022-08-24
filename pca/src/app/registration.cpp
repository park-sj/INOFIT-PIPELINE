#include <iostream>
#include <string>
#include <utility>
#include <Eigen/Dense>
#include "parameter_reader.hpp"
#include "target_manager.hpp"
#include "region_pca_manager.hpp"
#include "data_driven_nonrigid_icp.hpp"
#include "tools.hpp"

int main(int argc, char* argv[]) {
    ParameterReader params;
    TargetManager tm;
    RegionPCAManager rm;
    SimpleClock sc;

    if (argc != 2) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [parameter json file]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading parameters from " << std::string(argv[1]) << std::endl;
    if(!params.readForRegistration(std::string(argv[1]))) {
        std::cout << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reconstructing matrices from saved files" << std::endl;
    rm.read(params.m_matrix_file, params.smoothness_matrix_file, params.stiffness_matrix_file, params.region_vertex_mapping_file, params.colored_mean_mesh_file);

    sc.start();
    std::cout << "Reading target meshes" << std::endl;
    tm.addMeshes(params.target_files, params.target_face_direction, params.target_height_direction);
    std::cout << "Reading target meshes took " << sc.stop() << "s" << std::endl;

    if(params.is_target_aligned) {
        std::cout << "Target meshes are already aligned -> Alignment skipped" << std::endl;
    } else {
        sc.start();
        std::cout << "Aligning meshes" << std::endl;
        tm.noseTipAlign();
        std::cout << "Aligning meshes took " << sc.stop() << "s" << std::endl;
    }

    auto aligned_target_meshes = tm.getMeshes();

    for(int i = 0; i < aligned_target_meshes.size(); i++) {
        DataDrivenNonrigidICP icp;
        
        std::cout << "Registration for target mesh " << params.target_files[i] << std::endl;

        FileNameEditor fe;
        std::string output_file_prefix = fe.extractName(params.target_files[i]) + "_" + fe.extractName(std::string(argv[1]));
        output_file_prefix = fe.addPath(params.output_file_directory, output_file_prefix);

        sc.start();
        std::cout << "Initializing ICP" << std::endl;
        icp.init(rm, *(aligned_target_meshes[i]), params.region_weights);
        std::cout << "Initializing ICP took " << sc.stop() << "s" << std::endl;

        for(int is = 0; is < params.n_is; is++) {
            sc.start();
            std::cout << "Iteration set (" << params.iterations_per_is << " iterations) " << is << ": " << std::endl;
            double current_error = icp.iterateUntil(params.threshold_error, params.iterations_per_is,
                params.threshold_distance_per_is[is], params.smoothness_weight_per_is[is], params.stiffness_weight_per_is[is]);
            
            double elapsed = sc.stop();

            std::cout << "Current status (with threshold distance " << params.threshold_distance_per_is[is] << ")" << std::endl;
            icp.printStatus(params.threshold_distance_per_is[is], true);
            
            std::cout << "Iteration set " << is << " took " << elapsed << "s" << std::endl;

            if(params.do_write_intermediate_mesh == true) {
                std::cout << "Generating intermediate mesh file" << std::endl;
                auto output_mesh = icp.getResultMesh();
                if(!params.is_template_aligned) tm.alignBackToTargetMesh(output_mesh, i);
                tm.saveMesh(output_mesh, output_file_prefix + "_is" + std::to_string(is) + ".off", params.do_remove_color);
            }

            if(current_error < params.threshold_error) break;
        }
        std::cout << "ALL iteration sets finished" << std::endl;
/*
        double s;
        Eigen::Matrix3d rot;
        Eigen::VectorXd m;
        Eigen::Vector3d t;
        icp.getParameters(s, rot, m, t);
        std::cout << "Parameters:" << std::endl;
        std::cout << "(1) s" << std::endl;
        std::cout << s << std::endl;
        std::cout << "(2) rot" << std::endl;
        std::cout << rot << std::endl;
        std::cout << "(3) m" << std::endl;
        std::cout << m << std::endl;
        std::cout << "(4) t" << std::endl;
        std::cout << t << std::endl;
*/
        std::cout << "Writing to " << output_file_prefix + "_final.off" << std::endl;
        auto output_mesh = icp.getResultMesh();
        if(!params.is_target_aligned) tm.alignBackToTargetMesh(output_mesh, i);
        tm.saveMesh(output_mesh, output_file_prefix + "_final.off", params.do_remove_color);
    }

    std::cout << "Registration finished" << std::endl;

    return EXIT_SUCCESS;
}
