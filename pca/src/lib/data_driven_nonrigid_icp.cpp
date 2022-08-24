#include "data_driven_nonrigid_icp.hpp"

DataDrivenNonrigidICP::DataDrivenNonrigidICP() {}

// Initializes the DataDrivenNonrigidICP object
// Target mesh needs not to have the same topology with template meshes
// Assumes that
// (1) all template meshes have the same topology
// (2) template meshes are properly scaled and aligned
// (3) the target mesh is properly scaled and aligned compared to template meshes
void DataDrivenNonrigidICP::init(RegionPCAManager& rm, Mesh& target_mesh, std::vector<double>& region_weights) {
    this->rm = &rm;

    int n_k_n_r = this->rm->getNkNr();

    s = 1.0;
    rot = Eigen::Matrix3d::Identity(3, 3);
    m = Eigen::VectorXd::Zero(n_k_n_r);
    t = Eigen::Vector3d(0.0, 0.0, 0.0);
    current_mesh_vector = this->rm->getMeanVector();

    RegionTools rt;
    Mesh reference_mesh = this->rm->getColoredMeanMesh();
    computeVertexWeights(region_weights);

    mc.build_AABB(target_mesh);

    return;
}

// Iterates until the error goes below the threshold_error or until the iteration reaches max_iteration
// Returns the resulting error
double DataDrivenNonrigidICP::iterateUntil(double threshold_error, int max_iteration, double threshold_distance, double per_iteration_smoothness_weight, double per_iteration_stiffness_weight) {
    double error;

    for(int i = 0; i < max_iteration; i++) {
        iterate(threshold_distance, per_iteration_smoothness_weight, per_iteration_stiffness_weight);
        std::cout << "* " << std::flush;
        error = computeMeanSquaredError();
        if(error < threshold_error) {
            break;
        }
    }
    std::cout << std::endl;

    return error;
}

// Prints per-region total error, partial error (only for point pairs within threshold distance), and the ratio of point pairs within threshold distance
void DataDrivenNonrigidICP::printStatus(double threshold_distance, bool overall_only) {
    int n_d = current_mesh_vector.size();
    int n_v = n_d/3;
    int n_r = rm->getRegionVertexMapping()->size();
    
    std::vector<double> squared_distance_errors;

    // Computing squared distance error for all vertices
    for(int i = 0; i < n_d; i += 3) {
        double squared_distance_error = (current_mesh_vector.segment<3>(i) + t - correspondence_vector.segment<3>(i)).squaredNorm();
        squared_distance_errors.push_back(squared_distance_error);
    }

    if(!overall_only) {
        std::vector<double> total_sums(n_r, 0.0);
        std::vector<double> partial_sums(n_r, 0.0);
        std::vector<int> close_points_counts(n_r, 0);

        // Computing total sum, partial sum, and close points count for each region
        for(int i = 0; i < n_r; i++) {
            for(auto it = (*(rm->getRegionVertexMapping()))[i].begin(); it != (*(rm->getRegionVertexMapping()))[i].end(); ++it) {
                double squared_distance_error = squared_distance_errors[*it];
                if(squared_distance_error < threshold_distance * threshold_distance) {
                    partial_sums[i] += squared_distance_error;
                    close_points_counts[i] += 1;
                }
                total_sums[i] += squared_distance_error;
            }
        }

        // Printing summary for regions
        for(int i = 0; i < n_r; i++) {
            std::cout << "Region " << i << ":" << std::endl;
            std::cout << "Total MSE: " << total_sums[i] / (*(rm->getRegionVertexMapping()))[i].size()
                << ", Partial MSE: " << partial_sums[i] / close_points_counts[i]
                << ", Ratio of close points: " << double(close_points_counts[i]) / (*(rm->getRegionVertexMapping()))[i].size()
                << " (" << close_points_counts[i] << "/" << (*(rm->getRegionVertexMapping()))[i].size() << ")" << std::endl;
        }
    }

    // Computing total sum, partial sum, and close points count of the whole mesh
    double overall_total_sum = 0.0;
    double overall_partial_sum = 0.0;
    double overall_close_points_count = 0;

    for(int i = 0; i < n_v; i++) {
        double squared_distance_error = squared_distance_errors[i];
        if(squared_distance_error < threshold_distance * threshold_distance) {
            overall_partial_sum += squared_distance_error;
            overall_close_points_count += 1;
        }
        overall_total_sum += squared_distance_error;
    }

    // Printing an overall summary
    std::cout << "Overall:" << std::endl;
    std::cout << "Total MSE: " << overall_total_sum/n_v
        << ", Partial MSE: " << overall_partial_sum/overall_close_points_count
        << ", Ratio of close points: " << double(overall_close_points_count)/n_v
        << " (" << overall_close_points_count << "/" << n_v << ")" << std::endl;
    
    return;
}

// Returns current mesh
DataDrivenNonrigidICP::Mesh DataDrivenNonrigidICP::getResultMesh() {
    MeshTools mt;
    Mesh output_mesh = rm->getColoredMeanMesh();
    mt.meshVectorToMesh(current_mesh_vector, output_mesh);
    return output_mesh;
}

// Gets parameters
void DataDrivenNonrigidICP::getParameters(double& s, Eigen::Matrix3d& rot, Eigen::VectorXd& m, Eigen::Vector3d& t) {
    s = this->s;
    rot = this->rot;
    m = this->m;
    t = this->t;
    return;
}

// Computes per-vertex weights from per-region weights
// Assumes that member variables currnet_mesh_vector and region_vertex_mapping are properly set
// Also assumes that any vertex should be included in at least one region
void DataDrivenNonrigidICP::computeVertexWeights(std::vector<double>& region_weights) {
    int n_d = current_mesh_vector.size();
    int n_r = rm->getRegionVertexMapping()->size();
    int n_v = n_d/3;

    // Computing vertex-region mapping using region-vertex mapping
    std::vector<std::list<int> > vertex_region_mapping;
    for(int i = 0; i < n_v; i++) {
        vertex_region_mapping.push_back(std::list<int>());
    }

    for(int i = 0; i < n_r; i++) {
        for(auto it = (*(rm->getRegionVertexMapping()))[i].begin(); it != (*(rm->getRegionVertexMapping()))[i].end(); ++it) {
            vertex_region_mapping[*it].push_back(i);
        }
    }

    // Computing per-vertex weights from per-region weights
    vertex_weights.clear();
    for(int i = 0; i < n_v; i++) {
        double sum = 0.0;
        for(auto it = vertex_region_mapping[i].begin(); it != vertex_region_mapping[i].end(); ++it) {
            sum += region_weights[*it];
        }
        vertex_weights.push_back(sum/vertex_region_mapping[i].size());
    }

    return;
}

// Iterates to accumulatively register the template surface on target point cloud
// Per-iteration smoothness weight or stiffness weight (global) can be applied to implement varying smoothness and stiffness
void DataDrivenNonrigidICP::iterate(double threshold_distance, double per_iteration_smoothness_weight, double per_iteration_stiffness_weight) {
    int n_k_n_r = rm->getNkNr();
    //SimpleClock sc;
    
    // Finding correspondence between p+t and q
    //sc.start();
    computeCorrespondence();
    //sc.stopAndPrint("Computing correspondence");

    // Getting regularization matrices in advance
    // Even if the smoothness matrix has zero row (empty), it works correctly
    Eigen::MatrixXd smoothness = rm->getSmoothnessRegularization();
    Eigen::MatrixXd stiffness = rm->getStiffnessRegularization();

    // Building the matrix K and vector b
    int n_d = current_mesh_vector.size();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_d + smoothness.rows() + stiffness.rows(), 7 + n_k_n_r);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n_d + smoothness.rows() + stiffness.rows());

    //K.col(0).segment(0, n_d) = current_mesh_vector;
    for(int i = 0; i < n_d; i += 3) {
        // Rejecting outliers
        Eigen::VectorXd dist = correspondence_vector.segment<3>(i) - current_mesh_vector.segment<3>(i) - t;
        if(dist.squaredNorm() > threshold_distance * threshold_distance) {
            continue;
        }

        double vw = vertex_weights[i/3];
        b.segment<3>(i) = dist * vw;
        //K.block<3, 3>(i, 1) = -1 * getRotationForm(current_mesh_vector.segment<3>(i)) * vw;
        //K.block<3, 3>(i, 4) = Eigen::Matrix3d::Identity(3, 3) * vw;
        K.block(i, 7, 3, n_k_n_r) = s * rot * rm->getMBlock(i/3) * vw;
    }

    // Appending regularization matrices to K
    // Per-iteration weights are multiplied here
    K.block(n_d, 7, smoothness.rows(), n_k_n_r) = per_iteration_smoothness_weight * smoothness;
    K.block(n_d + smoothness.rows(), 7, stiffness.rows(), n_k_n_r) = per_iteration_stiffness_weight * stiffness;

    // Ensuring that K and b has no Nan or Inf
    EigenTools et;
    assert(et.hasNaN(K) == false && et.hasInf(K) == false);
    assert(et.hasNaN(b) == false && et.hasInf(b) == false);

    //sc.start();
    
    // Solving Kx = b with K^T*Kx = K^T*b using fullPivLu()
    Eigen::MatrixXd KTK = K.transpose() * K;
    Eigen::VectorXd KTb = K.transpose() * b;
    Eigen::VectorXd x = KTK.fullPivLu().solve(KTb);

    //sc.stopAndPrint("Solving the linear system");

    // Updating parameters s, rot, m, t, and new p (current_mesh_vector)
    update(x(0), x.segment<3>(1), x.segment(7, n_k_n_r), x.segment<3>(4));

    return;
}

// Computes current correspondence between template vertices and target points
// Current actual template vertices can be computed by adding t (translation)
void DataDrivenNonrigidICP::computeCorrespondence() {
    int n_d = current_mesh_vector.size();
    Eigen::VectorXd actual_points = current_mesh_vector;

    for(int i = 0; i < n_d; i += 3) {
        actual_points.segment<3>(i) += t;
    }

    correspondence_vector = mc.computeOnlyVP_AABB(actual_points);
    return;
}

// Updates parameters (scale, rotation, shape, and translation) according to the result of iteration
// Also updates current_mesh_vector with the updated parameters d, rot, m, and t
void DataDrivenNonrigidICP::update(double ds, Eigen::Vector3d dtheta, Eigen::VectorXd dm, Eigen::VectorXd dt) {
    // Updating parameters
    //s = (1 + ds) * s;
    //rot = (Eigen::Matrix3d::Identity(3, 3) + getRotationForm(dtheta)) * rot;
    m = dm + m;
    //t = dt + t;

    // Updating current mesh vertices
    int n_d = current_mesh_vector.size();
    Eigen::VectorXd mean_vector = rm->getMeanVector();
    for(int i = 0; i < n_d; i += 3) {
        current_mesh_vector.segment<3>(i) = s * rot * (rm->getMBlock(i/3) * m + mean_vector.segment<3>(i));
    }

    return;
}

// Computes total error (average of the squared norm of (actual point - corresponding target point)
double DataDrivenNonrigidICP::computeMeanSquaredError() {
    int n_d = current_mesh_vector.size();
    double sum = 0.0;

    for(int i = 0; i < n_d; i += 3) {
        double squared_distance_error = (current_mesh_vector.segment<3>(i) + t - correspondence_vector.segment<3>(i)).squaredNorm();
        sum += squared_distance_error;
    }

    return sum/(n_d/3);
}

// Returns the rotation-like matrix form (which is -[p]_x according to the terms of the paper)
Eigen::Matrix3d DataDrivenNonrigidICP::getRotationForm(Eigen::Vector3d v) {
    Eigen::Matrix3d mat = Eigen::Matrix3d::Zero(3, 3);
    mat(1, 2) = -v(0);
    mat(2, 1) = v(0);
    mat(0, 2) = v(1);
    mat(2, 0) = -v(1);
    mat(0, 1) = -v(2);
    mat(1, 0) = v(2);
    return mat;
}