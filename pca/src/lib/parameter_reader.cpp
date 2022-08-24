#include "parameter_reader.hpp"

ParameterReader::ParameterReader() {}

// Reads preprocessing parameters from json file of specific format
// Returns false if reading process fails
// Stored parameters can be accessed freely from the outside of the class object (public member variables)
bool ParameterReader::readForPreprocessing(std::string json_file) {
    std::ifstream json_ifstream(json_file);
    Json::CharReaderBuilder builder;
    Json::Value value;
    JSONCPP_STRING errs;

    bool ret = Json::parseFromStream(builder, json_ifstream, &value, &errs);
    json_ifstream.close();
    if(ret == false) return ret;

    // Attributes related to template files
    n_m = value["NumberOfSourceFiles"].asInt();
    std::string source_file_dir = value["SourceFilesDirectory"].asString();
    for(const auto & entry : std::filesystem::recursive_directory_iterator(source_file_dir))
    {
        if(entry.is_directory()) continue;
        if(entry.path().extension().compare(".off") == 0) template_files.push_back(entry.path().string());
    }
    assert(n_m > 0 && n_m == template_files.size());
    is_template_aligned = value["IsAligned"].asBool();
    template_face_direction = Direction(value["FaceDirection"].asString());
    template_height_direction = Direction(value["HeightDirection"].asString());

    // Attributes related to mesh regions
    n_r_initial = value["NumberOfInitialRegions"].asInt();
    for(int i = 0; i < value["RegionsToRemove"].size(); i++) {
        int r = value["RegionsToRemove"][i].asInt();
        assert(r >= 0 && r < n_r_initial);
        regions_to_remove.insert(r);
    }
    n_r_removed = value["NumberOfRegionsAfterRemoval"].asInt();
    assert(n_r_removed > 0 && n_r_removed <= n_r_initial);
    for(int i = 0; i < value["RegionsToMerge"].size(); i++) {
        regions_to_merge.push_back(std::set<int>());
        for(int j = 0; j < value["RegionsToMerge"][i].size(); j++) {
            int r = value["RegionsToMerge"][i][j].asInt();
            assert(r >= 0 && r < n_r_initial);
            assert(regions_to_remove.find(r) == regions_to_remove.end());
            regions_to_merge[i].insert(r);
        }
    }
    n_r = value["NumberOfRegionsAfterMerge"].asInt();
    assert(n_r > 0 && n_r <= n_r_removed);
    additional_region_overlap = value["AdditionalRegionOverlap"].asInt();
    assert(additional_region_overlap >= 0);

    // Attributes related to the construction of matrix M
    smoothness_factor = value["SmoothnessFactor"].asDouble();
    assert(n_r == value["StiffnessFactorPerRegion"].size());
    for(int i = 0; i < n_r; i++) {
        stiffness_factors.push_back(value["StiffnessFactorPerRegion"][i].asDouble());
    }
    n_k = value["NumberOfPrincipalComponents"].asInt();
    assert(n_k > 0 && n_k <= n_m);

    // Attributes related to the saving of matrix M and color reference mean mesh
    matrix_file_directory = value["MatrixFileDirectory"].asString();
    region_file_directory = value["RegionFileDirectory"].asString();
    mesh_file_directory = value["MeshFileDirectory"].asString();

    return ret;
}

// Reads registration parameters from json file of specific format
// Returns false if reading process fails
// Stored parameters can be accessed freely from the outside of the class object (public member variables)
bool ParameterReader::readForRegistration(std::string json_file) {
    std::ifstream json_ifstream(json_file);
    Json::CharReaderBuilder builder;
    Json::Value value;
    JSONCPP_STRING errs;

    bool ret = Json::parseFromStream(builder, json_ifstream, &value, &errs);
    json_ifstream.close();
    if(ret == false) return ret;

    // Attributes related to the loading of matrix M and color reference mean mesh
    m_matrix_file = value["MMatrixFile"].asString();
    smoothness_matrix_file = value["SmoothnessMatrixFile"].asString();
    stiffness_matrix_file = value["StiffnessMatrixFile"].asString();
    region_vertex_mapping_file = value["RegionVertexMappingFile"].asString();
    colored_mean_mesh_file = value["ColoredMeanMeshFile"].asString();

    // Attributes related to target files
    n_o = value["NumberOfTargetFiles"].asInt();
    std::string target_file_dir = value["TargetFilesDirectory"].asString();
    for(const auto & entry : std::filesystem::recursive_directory_iterator(target_file_dir))
    {
        if(entry.is_directory()) continue;
        if(entry.path().extension().compare(".off") == 0) target_files.push_back(entry.path().string());
    }
    assert(n_o> 0 && n_o == target_files.size());
    is_target_aligned = value["IsAligned"].asBool();
    target_face_direction = Direction(value["FaceDirection"].asString());
    target_height_direction = Direction(value["HeightDirection"].asString());

    // Attribute related to region weights
    n_r = value["NumberOfRegions"].asInt();
    assert(n_r > 0 && n_r == value["RegionWeights"].size());
    for(int i = 0; i < n_r; i++) {
        region_weights.push_back(value["RegionWeights"][i].asDouble());
    }

    // Attributes related to iteration
    n_is = value["NumberOfIterationSets"].asInt();
    iterations_per_is = value["NumberOfIterationsPerIterationSet"].asInt();
    assert(n_is > 0 && iterations_per_is > 0);
    assert(n_is == value["SmoothnessWeightPerIterationSet"].size());
    assert(n_is == value["StiffnessWeightPerIterationSet"].size());
    assert(n_is == value["ThresholdDistancePerIterationSet"].size());
    for(int i = 0; i < n_is; i++) {
        smoothness_weight_per_is.push_back(value["SmoothnessWeightPerIterationSet"][i].asDouble());
        stiffness_weight_per_is.push_back(value["StiffnessWeightPerIterationSet"][i].asDouble());
        threshold_distance_per_is.push_back(value["ThresholdDistancePerIterationSet"][i].asDouble());
    }
    threshold_error = value["ThresholdError"].asDouble();

    // Attributes related to output files
    output_file_directory = value["OutputFileDirectory"].asString();
    do_write_intermediate_mesh = value["IntermediateOutput"].asBool();
    do_remove_color = value["RemoveColor"].asBool();

    return ret;
}