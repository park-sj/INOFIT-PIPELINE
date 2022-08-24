#include "mesh_pca.hpp"

MeshPCA::MeshPCA() {}

// Converts each mesh into a single column vector and stores them as a single matrix
// Assumes that the input meshes are properly scaled and aligned
void MeshPCA::inputMeshes(std::vector<Mesh*>& meshes) {
    MeshTools mt;
    
    // Constructing matrix B with meshes as column vectors
    B = mt.meshesToMeshMatrix(meshes);    

    // Making the column vectors of B centered and computing the mean vector
    mean_vector = mt.makeMeshMatrixCentered(B);

    return;
}

// Gets a mesh matrix as input
// Assumes that the input meshes are properly scaled and aligned
void MeshPCA::inputMeshMatrix(Eigen::MatrixXd mesh_matrix) {
    B = mesh_matrix;

    // Making the column vectors of B centered and computing the mean vector
    MeshTools mt;
    mean_vector = mt.makeMeshMatrixCentered(B);
    return;
}

// Computes eigenvectors of given input meshes
void MeshPCA::computeEigenvectors() {
    // Here B^T*B is used instead of B*B^T (the PCA transpose trick)
    Eigen::MatrixXd C = B.transpose() * B;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(C);
    if(solver.info() != Eigen::Success) abort();

    // Computing unit eigenvectors (the columns of P are eigenvectors)
    P = B * solver.eigenvectors();  // Multiply B from left since the transpose trick is used
    P = P.rowwise().reverse().eval();  // The order of columns is reversed
                                       // to ensure that the columns are in eigenvalue-decreasing order
    P.colwise().normalize();

    // Storing computed eigenvalues (n_m eigenvalues as a column vector)
    S = solver.eigenvalues();
    S = S.reverse().eval();  // The order of eigenvalues is reversed to be decreasing order

    return;
}

// Approximates a mesh to an n_k-dimensional point
// The maximum possible value of n_k is the number of meshes used to compute PCA
// If n_k is greater than the possible maximum value, n_k is set to the maximum value
Eigen::VectorXd MeshPCA::approximate(Mesh& mesh, int n_k) {
    MeshTools mt;
    return approximate(mt.meshToMeshVector(mesh), n_k);
}

// Approximates a mesh vector to an n_k-dimensional point
// The maximum possible value of n_k is the number of meshes used to compute PCA
// If n_k is greater than the possible maximum value, n_k is set to the maximum value
Eigen::VectorXd MeshPCA::approximate(Eigen::VectorXd mesh_vector, int n_k) {
    if(n_k > B.cols()) n_k = B.cols();

    // Checking whether the size of mesh_vector is equal to 3 * the number of vertices of internal meshes
    assert(mesh_vector.size() == B.rows());

    // Centering the vector
    mesh_vector -= mean_vector;

    return P.transpose().block(0, 0, n_k, B.rows()) * mesh_vector;
}

// Reconstructs a mesh from the projection
// Assumes that the topology_mesh already has all edges and faces set appropriately
void MeshPCA::reconstructFrom(Eigen::VectorXd projection, Mesh& topology_mesh) {
    MeshTools mt;
    mt.meshVectorToMesh(reconstructFrom(projection), topology_mesh);
    return;
}

// Reconstructs a mesh vector from the projection
Eigen::VectorXd MeshPCA::reconstructFrom(Eigen::VectorXd projection) {
    // Checking whether the size of projection is not greater than the number of meshes used to compute PCA
    assert(projection.size() <= B.cols());

    int n_k = projection.size();
    return (P.block(0, 0, B.rows(), n_k) * projection) + mean_vector;
}

// Returns mean mesh vector
Eigen::VectorXd MeshPCA::getMeanVector() {
    return mean_vector;
}

// Sets the vertices of topology_mesh according to the mean mesh vector
void MeshPCA::getMeanMesh(Mesh& topology_mesh) {
    MeshTools mt;
    mt.meshVectorToMesh(mean_vector, topology_mesh);
    return;
}

// Returns a specific block of P matrix (a matrix whose columns are eigenvectors)
// Eigenvectors are stored in eigenvalue-decreasing order
Eigen::MatrixXd MeshPCA::getEigenvectorBlock(int i, int j, int r, int c) {
    int n_d = P.rows();
    int n_m = P.cols();

    assert(0 <= i && i < n_d);
    assert(0 <= j && j < n_m);
    assert(0 <= r && i+r <= n_d);
    assert(0 <= c && j+c <= n_m);

    return P.block(i, j, r, c);
}

// Returns matrix S (a column vector of eigenvalues)
// Eigenvalues are stored in decreasing order
Eigen::VectorXd MeshPCA::getEigenvalues() {
    return S;
}