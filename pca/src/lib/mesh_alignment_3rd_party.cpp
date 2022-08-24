#include "mesh_alignment_3rd_party.hpp"

MeshAlignment3rdParty::MeshAlignment3rdParty() {}

// Aligns source mesh to target mesh coarsely
// Implementation is based on OpenGR Super4PCS algorithm
// Parameters are being tested (not working well yet)
void MeshAlignment3rdParty::coarseAlignment(Mesh& source_mesh, Mesh& target_mesh) {
    // Pairs of points and normals
    std::vector<Pwn> pwns_source, pwns_target;

    // Computing normals
    // Computed normals are conveyed via xyz file
    std::string temp_file;
    temp_file = computeNormalsAsXYZ(source_mesh);
    std::ifstream source_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(source_xyz_ifstream, std::back_inserter(pwns_source), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    source_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    temp_file = computeNormalsAsXYZ(target_mesh);
    std::ifstream target_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(target_xyz_ifstream, std::back_inserter(pwns_target), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    target_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    // Coarse point set registration (alignment)
    CGAL::OpenGR::register_point_sets(pwns_target, pwns_source,
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap())
            .number_of_samples(500)
            .maximum_running_time(60)
            .accuracy(0.5)
            .maximum_normal_deviation(30),
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap())
    );
    
    // Writing the result point set to source mesh
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for(v_it = v_begin; v_it != v_end; ++v_it) {
        source_mesh.point(*v_it) = pwns_source[*v_it].first;
    }

    return;
}

// Aligns source mesh to target mesh finely
// Implementation is based on pointmatcher ICP algorithm
// Parameters are being tested (not working well yet)
bool MeshAlignment3rdParty::fineAlignment(Mesh& source_mesh, Mesh& target_mesh) {
    // Pairs of points and normals
    std::vector<Pwn> pwns_source, pwns_target;

    // Computing normals
    // Computed normals are conveyed via xyz file
    std::string temp_file;
    temp_file = computeNormalsAsXYZ(source_mesh);
    std::ifstream source_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(source_xyz_ifstream, std::back_inserter(pwns_source), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    source_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    temp_file = computeNormalsAsXYZ(target_mesh);
    std::ifstream target_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(target_xyz_ifstream, std::back_inserter(pwns_target), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    target_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    // Setting parameters
    using CGAL::pointmatcher::ICP_config;

    std::vector<ICP_config> point_set_target_filters;
    //point_set_target_filters.push_back(ICP_config{"MinDistDataPointsFilter", {{"minDist", "0.5"}}});
    //point_set_target_filters.push_back(ICP_config{"MaxDensityDataPointsFilter", {{"maxDensity", "1000000"}}});
    point_set_target_filters.push_back(ICP_config{"RandomSamplingDataPointsFilter", {{"prob", "0.05"}}});

    std::vector<ICP_config> point_set_source_filters;
    //point_set_source_filters.push_back(ICP_config{"MinDistDataPointsFilter", {{"minDist", "0.5"}}});
    //point_set_source_filters.push_back(ICP_config{"MaxDensityDataPointsFilter", {{"maxDensity", "1000000"}}});
    point_set_source_filters.push_back(ICP_config{"RandomSamplingDataPointsFilter", {{"prob", "0.05"}}});

    ICP_config matcher{"KDTreeMatcher", {{"knn", "1"}, {"epsilon", "3.0"}}};

    std::vector<ICP_config> outlier_filters;
    outlier_filters.push_back(ICP_config{"TrimmedDistOutlierFilter", {{"ratio", "0.75"}}});

    ICP_config error_minimizer{"PointToPointErrorMinimizer"};

    std::vector<ICP_config> transformation_checkers;
    transformation_checkers.push_back(ICP_config{"CounterTransformationChecker", {{"maxIterationCount", "150"}}});
    transformation_checkers.push_back(ICP_config{"DifferentialTransformationChecker", {{"minDiffRotErr", "0.001"}, {"minDiffTransErr", "0.01"}, {"smoothLength", "4"}}});

    ICP_config inspector{"NullInspector"};

    ICP_config logger{"NullLogger"};

    // Fine point set registration (alignment)
    bool converged = CGAL::pointmatcher::register_point_sets(pwns_target, pwns_source,
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap())
            .point_set_filters(point_set_target_filters)
            .matcher(matcher)
            .outlier_filters(outlier_filters)
            .error_minimizer(error_minimizer)
            .transformation_checkers(transformation_checkers)
            .inspector(inspector)
            .logger(logger),
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap())
            .point_set_filters(point_set_source_filters)
    );

    // Writing the result point set to source mesh
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for(v_it = v_begin; v_it != v_end; ++v_it) {
        source_mesh.point(*v_it) = pwns_source[*v_it].first;
    }

    return converged;
}

// Aligns source mesh to target mesh (coarsely first, then finely)
// Implementation is based on OpenGR Super4PCS algorithm + pointmatcher ICP algorithm
// Currently uses default parameters
bool MeshAlignment3rdParty::generalAlignment(Mesh& source_mesh, Mesh& target_mesh) {
    // Pairs of points and normals
    std::vector<Pwn> pwns_source, pwns_target;

    // Computing normals
    // Computed normals are conveyed via xyz file
    std::string temp_file;
    temp_file = computeNormalsAsXYZ(source_mesh);
    std::ifstream source_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(source_xyz_ifstream, std::back_inserter(pwns_source), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    source_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    temp_file = computeNormalsAsXYZ(target_mesh);
    std::ifstream target_xyz_ifstream(temp_file);
    CGAL::read_xyz_points(target_xyz_ifstream, std::back_inserter(pwns_target), CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()));
    target_xyz_ifstream.close();
    removeTemporaryFile(temp_file);

    // Coarse point set registration (alignment)
    K::Aff_transformation_3 T = std::get<0>(
        CGAL::OpenGR::compute_registration_transformation(pwns_target, pwns_source,
            CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()),
            CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()))
    );

    // Fine point set registration (alignment)
    bool converged = CGAL::pointmatcher::register_point_sets(pwns_target, pwns_source,
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()),
        CGAL::parameters::point_map(PointMap()).normal_map(NormalMap()).transformation(T)
    );

    // Writing the result point set to source mesh
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for(v_it = v_begin; v_it != v_end; ++v_it) {
        source_mesh.point(*v_it) = pwns_source[*v_it].first;
    }

    return converged;
}

// Computes vertex normals and stores the mesh as an xyz format file (6 coordinates "x y z nx ny nz" per line)
// The stored xyz file will be fed to alignment functions to deliver vertex normal data
std::string MeshAlignment3rdParty::computeNormalsAsXYZ(Mesh& mesh) {
    auto vnormals = mesh.add_property_map<VertexDescriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals);

    std::string output_file = "mesh/temp/temp.xyz";
    std::ofstream output_ofstream(output_file);
    
    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for(v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        output_ofstream << p[0] << " " << p[1] << " " << p[2] << " ";
        Vector n = vnormals[*v_it];
        output_ofstream << n[0] << " " << n[1] << " " << n[2] << std::endl;
    }
    output_ofstream.close();

    return output_file;
}

// Removes a file
// Used to remove temporary ply file which is used to deliver vertex normal data
void MeshAlignment3rdParty::removeTemporaryFile(std::string filename) {
    std::remove(filename.c_str());
    return;
}