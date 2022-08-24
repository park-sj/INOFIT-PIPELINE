#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cassert>
#include <filesystem>
#include <json/json.h>
#include "mesh_tools.hpp"

#ifndef PARAMETER_READER_HPP
#define PARAMETER_READER_HPP

class ParameterReader {
public:
    ParameterReader();

    bool readForPreprocessing(std::string json_file);
    bool readForRegistration(std::string json_file);

public:
    // Attributes related to template files
    int n_m;
    std::vector<std::string> template_files;
    bool is_template_aligned;
    Direction template_face_direction;
    Direction template_height_direction;

    // Attributes related to mesh regions
    int n_r_initial;
    std::set<int> regions_to_remove;
    int n_r_removed;
    std::vector<std::set<int> > regions_to_merge;
    int n_r;
    int additional_region_overlap;

    // Attributes related to the construction of matrix M
    double smoothness_factor;
    std::vector<double> stiffness_factors;
    int n_k;

    // Attributes related to the saving of matrix M
    std::string matrix_file_directory;
    std::string region_file_directory;
    std::string mesh_file_directory;

    // Attributes related to the loading of matrix M
    std::string m_matrix_file;
    std::string smoothness_matrix_file;
    std::string stiffness_matrix_file;
    std::string region_vertex_mapping_file;
    std::string colored_mean_mesh_file;

    // Attributes related to target files
    int n_o;
    std::vector<std::string> target_files;
    bool is_target_aligned;
    Direction target_face_direction;
    Direction target_height_direction;

    // Attribute related to region weights
    // int n_r is used again
    std::vector<double> region_weights;

    // Attributes related to iteration
    int n_is;
    int iterations_per_is;
    std::vector<double> smoothness_weight_per_is;
    std::vector<double> stiffness_weight_per_is;
    std::vector<double> threshold_distance_per_is;
    double threshold_error;

    // Attributes related to output files
    std::string output_file_directory;
    bool do_write_intermediate_mesh;
    bool do_remove_color;
};

#endif