#include "region_pca_manager.hpp"

RegionPCAManager::RegionPCAManager() {}

// Given vector of meshes and region-vertex mapping,
// computes PCA for each region and constructs matrix M (paper notation) which is composed of the specific combination of eigenvectors from different regions
// Only first n_k eigenvectors are used as principal components
// Regularization matrices (smoothness and stiffness) are also constructed in this function
// Smoothness factor and stiffness factors (per region) are applied to their construction
// Assumes that
// (1) at least two vertices are included in each region
// (2) every vertex is assigned a region
// (3) two neighboring regions have some vertices in common
// (4) the vertex number starts from 0 and ends with [number of vertices - 1]
void RegionPCAManager::construct(std::vector<Mesh*>& meshes, std::vector<std::set<int> >& region_vertex_mapping, int n_k, double smoothness_factor, std::vector<double>& stiffness_factor) {
    MeshTools mt;
    
    // Constructing a mesh matrix with meshes as column vectors
    Eigen::MatrixXd mesh_matrix = mt.meshesToMeshMatrix(meshes);
    int n_d = mesh_matrix.rows();
    int n_m = mesh_matrix.cols();
    int n_r = region_vertex_mapping.size();
    int n_v = n_d/3;

    assert(n_r == stiffness_factor.size());

    // Making the column vectors of the mesh matrix centered and computing the mean vector
    mean_vector = mt.makeMeshMatrixCentered(mesh_matrix);

    // Computing vertex-region mapping using region-vertex mapping
    std::vector<std::list<int> > vertex_region_mapping;
    for(int i = 0; i < n_v; i++) {
        vertex_region_mapping.push_back(std::list<int>());
    }

    for(int i = 0; i < n_r; i++) {
        for(auto it = region_vertex_mapping[i].begin(); it != region_vertex_mapping[i].end(); ++it) {
            vertex_region_mapping[*it].push_back(i);
        }
    }

    // Counting the number of vertices that are included in exactly two regions
    // Vertices included in 3, 4, 5, ... regions are ignored (as paper mentioned)
    // (preperation for the smoothness matrix construction)
    std::vector<int> vertices_in_two_regions;
    for(int i = 0; i < n_v; i++) {
        if(vertex_region_mapping[i].size() == 2) {
            vertices_in_two_regions.push_back(i);
        }
    }

    // Constructing the matrix M, the smoothness matrix, and the stiffness matrix at the same time
    M = Eigen::MatrixXd::Zero(n_d, n_k * n_r);
    smoothness = Eigen::MatrixXd::Zero(vertices_in_two_regions.size() * 3, n_k * n_r);
    stiffness = Eigen::MatrixXd::Zero(n_k * n_r, n_k * n_r);
    for(int i = 0; i < n_r; i++) {
        int region_size = region_vertex_mapping[i].size();
        Eigen::MatrixXd region_matrix(region_size * 3, n_m);

        // Constructing a region matrix which has only the rows whose corresponding vertices are in the region
        int idx_rm = 0;
        for(auto it_rv = region_vertex_mapping[i].begin(); it_rv != region_vertex_mapping[i].end(); ++it_rv) {
            int p = *it_rv;
            region_matrix.block(idx_rm*3, 0, 3, n_m) = mesh_matrix.block(p*3, 0, 3, n_m);
            ++idx_rm;
        }

        // Doing PCA on the region matrix
        MeshPCA pca;
        pca.inputMeshMatrix(region_matrix);
        pca.computeEigenvectors();

        // Constructing the stiffness matrix here
        Eigen::VectorXd eigenvalues = pca.getEigenvalues();
        for(int j = 0; j < n_k; j++) {
            stiffness(i*n_k + j, i*n_k + j) = 1.0 / sqrt(eigenvalues(j)) * stiffness_factor[i];
        }

        // Constructing the matrix M and smoothness here
        int idx_ev = 0;
        for(auto it_rv = region_vertex_mapping[i].begin(); it_rv != region_vertex_mapping[i].end(); ++it_rv) {
            int p = *it_rv;

            // Copying the computed eigenvectors to appropriate block of matrix M
            // If a vertex are included in more than one region,
            // its eigenvectors are divided by the number of its regions to balance the effect among vertices
            M.block(p*3, i*n_k, 3, n_k) = pca.getEigenvectorBlock(idx_ev*3, 0, 3, n_k) / double(vertex_region_mapping[p].size());

            // Copying the computed eigenvectors to appropriate block of smoothness matrix
            // Only the vertices included in exactly two regions are considered
            if(vertex_region_mapping[p].size() == 2) {
                auto it_tr = std::find(vertices_in_two_regions.begin(), vertices_in_two_regions.end(), p);
                int idx_sm = std::distance(vertices_in_two_regions.begin(), it_tr);
                assert(idx_sm < vertices_in_two_regions.size());

                if(*(vertex_region_mapping[p].begin()) == i) {
                    smoothness.block(idx_sm*3, i*n_k, 3, n_k) = pca.getEigenvectorBlock(idx_ev*3, 0, 3, n_k) * smoothness_factor;
                } else {
                    smoothness.block(idx_sm*3, i*n_k, 3, n_k) = (-1.0) * pca.getEigenvectorBlock(idx_ev*3, 0, 3, n_k) * smoothness_factor;
                }
            }

            ++idx_ev;
        }
    }

    // Storing mean mesh with color
    colored_mean_mesh = *(meshes[0]);
    mt.meshVectorToMesh(mean_vector, colored_mean_mesh);

    // Storing region-vertex mapping
    this->region_vertex_mapping = region_vertex_mapping;

    return;
}

// Given a vertex id, returns a corresponding block of matrix M
// The block size is 3 x (n_k * n_r)
Eigen::MatrixXd RegionPCAManager::getMBlock(int v) {
    return M.block(v*3, 0, 3, M.cols());
}

// Returns smoothness regularization matrix
Eigen::MatrixXd RegionPCAManager::getSmoothnessRegularization() {
    return smoothness;
}

// Returns stiffness regularization matrix
Eigen::MatrixXd RegionPCAManager::getStiffnessRegularization() {
    return stiffness;
}

// Returns a pointer to the region vertex mapping vector
std::vector<std::set<int> >* RegionPCAManager::getRegionVertexMapping() {
    return &region_vertex_mapping;
}

// Returns the mean vector
Eigen::VectorXd RegionPCAManager::getMeanVector() {
    return mean_vector;
}

// Returns the reference mean mesh
RegionPCAManager::Mesh RegionPCAManager::getColoredMeanMesh() {
    return colored_mean_mesh;
}

// Returns the value of n_k * n_r
int RegionPCAManager::getNkNr() {
    return M.cols();
}

// Saves the matrix M, smoothness matrix, stiffness matrix, region-vertex mapping, and colored mean mesh to files
void RegionPCAManager::save(std::string f_M, std::string f_smoothness, std::string f_stiffness, std::string f_rvm, std::string f_cmm) {
    saveMAsTxt(f_M);
    saveSmoothnessAsTxt(f_smoothness);
    saveStiffnessAsTxt(f_stiffness);
    saveRegionVertexMappingAsTxt(f_rvm);
    saveMeanAsOffWithColor(f_cmm);
    return;
}

// Reads the matrix M, smoothness matrix, stiffness matrix, region-vertex mapping, and colored mean mesh from files
void RegionPCAManager::read(std::string f_M, std::string f_smoothness, std::string f_stiffness, std::string f_rvm, std::string f_cmm) {
    readMFromTxt(f_M);
    readSmoothnessFromTxt(f_smoothness);
    readStiffnessFromTxt(f_stiffness);
    readRegionVertexMappingFromTxt(f_rvm);
    readColoredMeanFromOff(f_cmm);
    return;
}

// Saves the matrix M as a txt file
void RegionPCAManager::saveMAsTxt(std::string file_name) {
    EigenTools et;
    et.saveMatrixAsTxt(M, file_name);
    return;
}

// Saves the smoothness matrix as a txt file
void RegionPCAManager::saveSmoothnessAsTxt(std::string file_name) {
    EigenTools et;
    et.saveMatrixAsTxt(smoothness, file_name);
    return;
}

// Saves the stiffness matrix as a txt file
void RegionPCAManager::saveStiffnessAsTxt(std::string file_name) {
    EigenTools et;
    et.saveMatrixAsTxt(stiffness, file_name);
    return;
}

// Saves the region-vertex mapping as a txt file
void RegionPCAManager::saveRegionVertexMappingAsTxt(std::string file_name) {
    std::ofstream output_ofstream(file_name);
    output_ofstream << region_vertex_mapping.size() << std::endl;
    for(int r = 0; r < region_vertex_mapping.size(); r++) {
        output_ofstream << region_vertex_mapping[r].size() << " ";
        for(auto it = region_vertex_mapping[r].begin(); it != region_vertex_mapping[r].end(); ++it) {
            output_ofstream << *it << " ";
        }
        output_ofstream << std::endl;
    }
    output_ofstream.close();
    return;
}

// Saves the mean vector as an off file with color
void RegionPCAManager::saveMeanAsOffWithColor(std::string file_name) {
    RegionTools rt;
    rt.writeCorrectOffFileWithFaceColor(colored_mean_mesh, file_name);
    return;
}

// Reads the matrix M from a txt file
void RegionPCAManager::readMFromTxt(std::string file_name) {
    EigenTools et;
    M = et.readMatrixFromTxt(file_name);
    return;
}

// Reads the smoothness matrix from a txt file
void RegionPCAManager::readSmoothnessFromTxt(std::string file_name) {
    EigenTools et;
    smoothness = et.readMatrixFromTxt(file_name);
    return;
}

// Reads the Stiffness matrix from a txt file
void RegionPCAManager::readStiffnessFromTxt(std::string file_name) {
    EigenTools et;
    stiffness = et.readMatrixFromTxt(file_name);
    return;
}

// Reads the region-vertex mapping as a txt file
void RegionPCAManager::readRegionVertexMappingFromTxt(std::string file_name) {
    std::ifstream input_ifstream(file_name);
    int n_r;
    input_ifstream >> n_r;
    for(int i = 0; i < n_r; i++) {
        region_vertex_mapping.push_back(std::set<int>());

        int n_s;
        input_ifstream >> n_s;
        for(int j = 0; j < n_s; j++) {
            int v;
            input_ifstream >> v;
            region_vertex_mapping[i].insert(v);
        }
    }
    input_ifstream.close();

    return;
}

// Reads the mean vector from an mean mesh off file
void RegionPCAManager::readColoredMeanFromOff(std::string file_name) {
    std::ifstream mesh_ifstream(file_name);
    mesh_ifstream >> colored_mean_mesh;
    mesh_ifstream.close();

    MeshTools mt;
    mean_vector = mt.meshToMeshVector(colored_mean_mesh);

    return;
}