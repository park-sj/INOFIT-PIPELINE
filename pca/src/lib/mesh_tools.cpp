#include "mesh_tools.hpp"

// Default constructor is defined but should be avoided
Direction::Direction() {
    d = "X";
}

// A direction object should be constructed with a initial direction string (default constructor should be avoided)
// The direction d should be one of "X", "Y", "Z", "-X", "-Y", "-Z"
Direction::Direction(std::string d) {
    this->d = d;
    assertCorrectDirection();
}

// Returns stored direction string
std::string Direction::getString() {
    return d;
}

// Checks if the stored direction is orthogonal to given direction
bool Direction::isOrthogonalTo(Direction direction) {
    return this->d != direction.getString() && this->d != direction.oppositeDirection().getString();
}

// Returns the opposite direction object
Direction Direction::oppositeDirection() {
    std::string opposite_direction;
    if(this->d[0] == '-') {
        opposite_direction = this->d.substr(1, 1);
    } else {
        opposite_direction = "-" + this->d;
    }

    return Direction(opposite_direction);
}

// Returns 90-degree clockwise rotated direction object about the axis
Direction Direction::clockwiseDirection(Direction axis) {
    if(!isOrthogonalTo(axis)) return Direction(this->d);

    std::string from, to;
    std::string a = axis.getString();
    if(a == "X") from = "Z", to = "Y";
    else if(a == "Y") from = "X", to = "Z";
    else if(a == "Z") from = "Y", to = "X";
    else if(a == "-X") from = "Y", to = "Z";
    else if(a == "-Y") from = "Z", to = "X";
    else if(a == "-Z") from = "X", to = "Y";
    else abort();

    std::string clockwise_direction;
    int idx = 0;
    bool minus = false;

    if(this->d[idx] == '-') {
        minus = true;
        idx += 1;
    }

    if(this->d[idx] == from[0]) clockwise_direction = minus ? "-" + to : to;
    else if(this->d[idx] == to[0]) clockwise_direction = minus ? from : "-" + from;
    else abort();

    return Direction(clockwise_direction);
}

// Returns 90-degree counterclockwise rotated direction object about the axis
Direction Direction::counterclockwiseDirection(Direction axis) {
    if(!isOrthogonalTo(axis)) return Direction(this->d);

    return clockwiseDirection(axis).oppositeDirection();
}

// Returns outer product direction object
// Assumes that the stored direction and direction_b are orthogonal
Direction Direction::outerProductDirection(Direction direction_b) {
    assert(isOrthogonalTo(direction_b));

    int number_of_minus = 0;
    std::string a, b, c;

    if(this->d[0] == '-') {
        number_of_minus += 1;
        a = this->d.substr(1, 1);
    } else {
        a = this->d;
    }

    if(direction_b.getString()[0] == '-') {
        number_of_minus += 1;
        b = direction_b.getString().substr(1, 1);
    } else {
        b = direction_b.getString();
    }

    if(a == "X" && b == "Y") {
        if(number_of_minus % 2 == 0) c = "Z";
        else c = "-Z";
    } else if(a == "Y" && b == "X") {
        if(number_of_minus % 2 == 0) c = "-Z";
        else c = "Z";
    } else if(a == "Y" && b == "Z") {
        if(number_of_minus % 2 == 0) c = "X";
        else c = "-X";
    } else if(a == "Z" && b == "Y") {
        if(number_of_minus % 2 == 0) c = "-X";
        else c = "X";
    } else if(a == "Z" && b == "X") {
        if(number_of_minus % 2 == 0) c = "Y";
        else c = "-Y";
    } else if(a == "X" && b == "Z") {
        if(number_of_minus % 2 == 0) c = "-Y";
        else c = "Y";
    } else abort();

    return Direction(c);
}

// Asserts that stored direction d is one of "X", "Y", "Z", "-X", "-Y", "-Z"
void Direction::assertCorrectDirection() {
    assert(d == "X" || d == "Y" || d == "Z" || d == "-X" || d == "-Y" || d == "-Z");
    return;
}

MeshTools::MeshTools() {}

// Converts a mesh into a corresponding mesh vector
Eigen::VectorXd MeshTools::meshToMeshVector(Mesh& mesh) {
    Eigen::VectorXd mesh_vector(3 * mesh.number_of_vertices());

    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    int idx = 0;
    for(v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        mesh_vector(idx++) = p[0];
        mesh_vector(idx++) = p[1];
        mesh_vector(idx++) = p[2];
    }

    return mesh_vector;
}

// Converts a mesh vector into a corresponding mesh
// Assumes that the topology_mesh already has all edges and faces set appropriately
void MeshTools::meshVectorToMesh(Eigen::VectorXd mesh_vector, Mesh& topology_mesh) {
    // Checking whether the size of mesh_vector is equal to 3 * the number of vertices of given mesh
    assert(3 * topology_mesh.number_of_vertices() == mesh_vector.size());

    Mesh::Vertex_range v_range = topology_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    int idx = 0;
    for(v_it = v_begin; v_it != v_end; ++v_it) {
        double x = mesh_vector(idx++);
        double y = mesh_vector(idx++);
        double z = mesh_vector(idx++);
        topology_mesh.point(*v_it) = Point(x, y, z);
    }
    
    return;
}

// Constructs mesh matrix with meshes as column vectors
Eigen::MatrixXd MeshTools::meshesToMeshMatrix(std::vector<Mesh*>& meshes) {
    int n_m = meshes.size();
    int n_v = meshes[0]->number_of_vertices();
    int n_d = 3*n_v;

    // Checking whether all given meshes have the same number of vertices
    for(int i = 1; i < meshes.size(); i++) {
        assert(meshes[i]->number_of_vertices() == n_v);
    }

    // Constructing mesh matrix with meshes as column vectors
    Eigen::MatrixXd mat(n_d, n_m);
    mat = Eigen::MatrixXd(n_d, n_m);
    for(int i = 0; i < n_m; i++) {
        mat.col(i) = meshToMeshVector(*(meshes[i]));
    }

    return mat;
}

// Converts a mesh matrix into a corresponding meshes
// Assumes that the topology_meshes already has meshes with all edges and faces set appropriately
void MeshTools::meshMatrixToMeshes(Eigen::MatrixXd& mesh_matrix, std::vector<Mesh*>& topology_meshes) {
    assert(mesh_matrix.cols() == topology_meshes.size());

    int n_m = mesh_matrix.cols();
    for(int i = 0; i < n_m; i++) {
        meshVectorToMesh(mesh_matrix.col(i), *(topology_meshes[i]));
    }

    return;
}

// Makes the mesh matrix centered
// Returns mean mesh vector
Eigen::VectorXd MeshTools::makeMeshMatrixCentered(Eigen::MatrixXd& mesh_matrix) {
    int n_d = mesh_matrix.rows();
    int n_m = mesh_matrix.cols();
    Eigen::VectorXd mean_vector = Eigen::VectorXd::Zero(n_d);

    assert(n_d % 3 == 0);

    // Computing the mean vector
    for(int i = 0; i < n_m; i++) {
        mean_vector += mesh_matrix.col(i);
    }
    mean_vector /= n_m;

    // Making the column vectors of mesh_matrix centered
    mesh_matrix.colwise() -= mean_vector;

    return mean_vector;
}

// Scales the height of each mesh in the mesh matrix to given h
// The direction of the height should be one of "X", "Y", "Z", "-X", "-Y", "-Z"
// Returned scale_back_factors can be used to recover the sizes of meshes
std::vector<double> MeshTools::scaleMeshMatrixHeight(Eigen::MatrixXd& mesh_matrix, double h, Direction height_direction) {
    std::vector<double> scale_back_factors;
    int n_d = mesh_matrix.rows();
    int n_m = mesh_matrix.cols();

    assert(n_d % 3 == 0);
    assert(h > 0.0);
    
    int height_idx;
    std::string d = height_direction.getString();
    if(d == "X" || d == "-X") height_idx = 0;
    else if(d == "Y" || d == "-Y") height_idx = 1;
    else if(d == "Z" || d == "-Z") height_idx = 2;
    else abort();

    for(int i = 0; i < n_m; i++) {
        // Computing the height of i-th column mesh vector
        double max_coordinate = INITIAL_MAX_VALUE;
        double min_coordinate = INITIAL_MIN_VALUE;
        for(int j = 0; j < n_d; j += 3) {
            double coordinate = mesh_matrix(j + height_idx, i);
            max_coordinate = coordinate > max_coordinate ? coordinate : max_coordinate;
            min_coordinate = coordinate < min_coordinate ? coordinate : min_coordinate;
        }

        double original_height = max_coordinate - min_coordinate;
        double scale = h / original_height;

        // Scaling all of x, y, z coordinates by the computed scale
        for(int j = 0; j < n_d; j += 3) {
            mesh_matrix(j, i) *= scale;
        }

        scale_back_factors.push_back(1.0 / scale);
    }

    return scale_back_factors;
}

// Rotates mesh matrix according to their face and height directions
// Assumes that meshes already have the same face and height directions
// The directions of face and height should be orthogonal
void MeshTools::rotateMeshMatrixAxis(Eigen::MatrixXd& mesh_matrix, Direction face_from, Direction height_from, Direction face_to, Direction height_to) {
    assert(face_from.isOrthogonalTo(height_from) && face_to.isOrthogonalTo(height_to));

    Direction current_face = face_from;
    Direction current_height = height_from;

    // Rotating face axis
    while(current_face.getString() != face_to.getString()) {
        Direction rotation_axis = current_height;
        if(current_face.oppositeDirection().getString() != face_to.getString()) {
            rotation_axis = current_face.outerProductDirection(face_to);
        }
        rotateMeshMatrixAxisOnce(mesh_matrix, current_face, current_face.counterclockwiseDirection(rotation_axis));
        current_face = current_face.counterclockwiseDirection(rotation_axis);
        current_height = current_height.counterclockwiseDirection(rotation_axis);
    }

    // Rotating height axis (face axes are already matched)
    while(current_height.getString() != height_to.getString()) {
        rotateMeshMatrixAxisOnce(mesh_matrix, current_height, current_height.counterclockwiseDirection(current_face));
        current_height = current_height.counterclockwiseDirection(current_face);
    }

    return;
}

// Translates each mesh vector in mesh matrix so that the average point is at the origin
// Returns the amount of translations of each mesh vector that can be used later to recover their initial positions
Eigen::MatrixXd MeshTools::translateMeshMatrixToOrigin(Eigen::MatrixXd& mesh_matrix) {
    int n_m = mesh_matrix.cols();
    int n_d = mesh_matrix.rows();
    Eigen::MatrixXd translations(3, n_m);

    assert(n_d % 3 == 0);

    for(int i = 0; i < n_m; i++) {
        Eigen::Vector3d mean_point(0.0, 0.0, 0.0);
        for(int j = 0; j < n_d; j += 3) {
            mean_point += mesh_matrix.col(i).segment<3>(j);
        }
        mean_point /= (double)(n_d / 3);
        for(int j = 0; j < n_d; j += 3) {
            mesh_matrix.col(i).segment<3>(j) -= mean_point;
        }
        translations.col(i) = mean_point;
    }

    return translations;
}

// Translates each mesh vector in mesh matrix so that the tip point is at the origin
// Returns the amount of translations of each mesh vector that can be used later to recover their initial positions
Eigen::MatrixXd MeshTools::translateMeshMatrixTipToOrigin(Eigen::MatrixXd& mesh_matrix, Direction d) {
    int n_m = mesh_matrix.cols();
    int n_d = mesh_matrix.rows();
    Eigen::MatrixXd translations(3, n_m);

    assert(n_d % 3 == 0);

    for(int i = 0; i < n_m; i++) {
        Eigen::Vector3d tip_point = mesh_matrix.col(i).segment<3>(0);
        if(d.getString() == "X") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(0) > tip_point(0)) tip_point = p;
            }
        } else if(d.getString() == "-X") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(0) < tip_point(0)) tip_point = p;
            }
        } else if(d.getString() == "Y") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(1) > tip_point(1)) tip_point = p;
            }
        } else if(d.getString() == "-Y") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(1) < tip_point(1)) tip_point = p;
            }
        } else if(d.getString() == "Z") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(2) > tip_point(2)) tip_point = p;
            }
        } else if(d.getString() == "-Z") {
            for(int j = 0; j < n_d; j += 3) {
                Eigen::Vector3d p = mesh_matrix.col(i).segment<3>(j);
                if(p(2) < tip_point(2)) tip_point = p;
            }
        } else {
            abort();
        }

        for(int j = 0; j < n_d; j += 3) {
            mesh_matrix.col(i).segment<3>(j) -= tip_point;
        }

        translations.col(i) = tip_point;
    }

    return translations;
}

// Translates each mesh in mesh matrix by each 3D column vector specified in translations matrix
void MeshTools::translateMeshMatrixEach(Eigen::MatrixXd& mesh_matrix, Eigen::MatrixXd translations) {
    int n_m = mesh_matrix.cols();
    int n_d = mesh_matrix.rows();

    assert(n_d % 3 == 0);
    assert(translations.rows() == 3 && translations.cols() == n_m);

    for(int i = 0; i < n_m; i++) {
        for(int j = 0; j < n_d; j += 3) {
            mesh_matrix.col(i).segment<3>(j) += translations.col(i);
        }
    }

    return;
}

// Translates all meshes in mesh matrix by given 3D column vector
void MeshTools::translateMeshMatrixAll(Eigen::MatrixXd& mesh_matrix, Eigen::Vector3d translation) {
    int n_m = mesh_matrix.cols();
    int n_d = mesh_matrix.rows();

    assert(n_d % 3 == 0);

    for(int i = 0; i < n_m; i++) {
        for(int j = 0; j < n_d; j += 3) {
            mesh_matrix.col(i).segment<3>(j) += translation;
        }
    }

    return;
}

// Makes the meshes centered
// Returns mean mesh vector
Eigen::VectorXd MeshTools::makeMeshesCentered(std::vector<Mesh*>& meshes) {
    Eigen::MatrixXd mesh_matrix = meshesToMeshMatrix(meshes);
    Eigen::VectorXd mean_vector = makeMeshMatrixCentered(mesh_matrix);

    for(int i = 0; i < meshes.size(); i++) {
        meshVectorToMesh(mesh_matrix.col(i), *(meshes[i]));
    }

    return mean_vector;
}

// Scales the height of each mesh to given h
// Returned scale_back_factors can be used to recover the sizes of meshes
std::vector<double> MeshTools::scaleMeshesHeight(std::vector<Mesh*>& meshes, double h, Direction height_direction) {
    std::vector<double> scale_back_factors;

    for(int i = 0; i < meshes.size(); i++) {
        double scale_back_factor = scaleMeshHeight(*(meshes[i]), h, height_direction);
        scale_back_factors.push_back(scale_back_factor);
    }

    return scale_back_factors;
}

// Rotates meshes according to their face and height directions
// Assumes that meshes already have the same face and height directions
// The directions of face and height should be orthogonal
void MeshTools::rotateMeshesAxis(std::vector<Mesh*>& meshes, Direction face_from, Direction height_from, Direction face_to, Direction height_to) {
    for(int i = 0; i < meshes.size(); i++) {
        rotateMeshAxis(*(meshes[i]), face_from, height_from, face_to, height_to);
    }
    return;
}

// Translates each mesh so that the average point is at the origin
// Returns the amount of translations of each mesh that can be used later to recover their initial positions
Eigen::MatrixXd MeshTools::translateMeshesToOrigin(std::vector<Mesh*>& meshes) {
    Eigen::MatrixXd translations(3, meshes.size());

    for(int i = 0; i < meshes.size(); i++) {
        translations.col(i) = translateMeshToOrigin(*(meshes[i]));
    }

    return translations;
}

// Translates each mesh so that the tip point is at the origin
// Returns the amount of translations of each mesh that can be used later to recover their initial positions
Eigen::MatrixXd MeshTools::translateMeshesTipToOrigin(std::vector<Mesh*>& meshes, Direction d) {
    Eigen::MatrixXd translations(3, meshes.size());

    for(int i = 0; i < meshes.size(); i++) {
        translations.col(i) = translateMeshTipToOrigin(*(meshes[i]), d);
    }

    return translations;
}

// Translates each mesh by each 3D column vector specified in translations matrix
void MeshTools::translateMeshesEach(std::vector<Mesh*>& meshes, Eigen::MatrixXd translations) {
    assert(translations.rows() == 3 && translations.cols() == meshes.size());

    for(int i = 0; i < meshes.size(); i++) {
        translateMesh(*(meshes[i]), translations.col(i));
    }

    return;
}

// Translates all meshes by given 3D column vector
void MeshTools::translateMeshesAll(std::vector<Mesh*>& meshes, Eigen::Vector3d translation) {
    for(int i = 0; i < meshes.size(); i++) {
        translateMesh(*(meshes[i]), translation);
    }

    return;
}

// Scales the height of mesh to given h
// Returned "1.0 / scale" can be used to recover the size of the mesh
double MeshTools::scaleMeshHeight(Mesh& mesh, double h, Direction height_direction) {
    assert(h > 0.0);
    
    int height_idx;
    std::string d = height_direction.getString();
    if(d == "X" || d == "-X") height_idx = 0;
    else if(d == "Y" || d == "-Y") height_idx = 1;
    else if(d == "Z" || d == "-Z") height_idx = 2;
    else abort();

    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    double max_coordinate = INITIAL_MAX_VALUE;
    double min_coordinate = INITIAL_MIN_VALUE;

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double coordinate = (mesh.point(*v_it))[height_idx];
        max_coordinate = coordinate > max_coordinate ? coordinate : max_coordinate;
        min_coordinate = coordinate < min_coordinate ? coordinate : min_coordinate;
    }

    double original_height = max_coordinate - min_coordinate;
    double scale = h / original_height;

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        double x = p[0] * scale;
        double y = p[1] * scale;
        double z = p[2] * scale;
        mesh.point(*v_it) = Point(x, y, z);
    }

    return 1.0 / scale;
}

// Rotates a mesh according to its face and height directions
// The directions of face and height should be orthogonal
void MeshTools::rotateMeshAxis(Mesh& mesh, Direction face_from, Direction height_from, Direction face_to, Direction height_to) {
    assert(face_from.isOrthogonalTo(height_from) && face_to.isOrthogonalTo(height_to));

    Direction current_face = face_from;
    Direction current_height = height_from;

    // Rotating face axis
    while(current_face.getString() != face_to.getString()) {
        Direction rotation_axis = current_height;
        if(current_face.oppositeDirection().getString() != face_to.getString()) {
            rotation_axis = current_face.outerProductDirection(face_to);
        }
        rotateMeshAxisOnce(mesh, current_face, current_face.counterclockwiseDirection(rotation_axis));
        current_face = current_face.counterclockwiseDirection(rotation_axis);
        current_height = current_height.counterclockwiseDirection(rotation_axis);
    }

    // Rotating height axis (face axes are already matched)
    while(current_height.getString() != height_to.getString()) {
        rotateMeshAxisOnce(mesh, current_height, current_height.counterclockwiseDirection(current_face));
        current_height = current_height.counterclockwiseDirection(current_face);
    }

    return;
}

// Translates given mesh so that the average point is at the origin
// Returns the amount of translation that can be used later to recover the initial position
Eigen::Vector3d MeshTools::translateMeshToOrigin(Mesh& mesh) {
    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    Eigen::Vector3d mean_point(0.0, 0.0, 0.0);
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        mean_point += Eigen::Vector3d(p[0], p[1], p[2]);
    }
    mean_point /= mesh.num_vertices();
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        double x = p[0] - mean_point(0);
        double y = p[1] - mean_point(1);
        double z = p[2] - mean_point(2);
        mesh.point(*v_it) = Point(x, y, z);
    }

    return mean_point;
}

// Translates given mesh so that the tip point is at the origin
// Returns the amount of translation that can be used later to recover the initial position
Eigen::Vector3d MeshTools::translateMeshTipToOrigin(Mesh& mesh, Direction d) {
    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    int tip_point_index = 0;
    if(d.getString() == "X") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[0] > tip[0]) tip_point_index = int(*v_it);
        }
    } else if(d.getString() == "-X") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[0] < tip[0]) tip_point_index = int(*v_it);
        }
    } else if(d.getString() == "Y") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[1] > tip[1]) tip_point_index = int(*v_it);
        }
    } else if(d.getString() == "-Y") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[1] < tip[1]) tip_point_index = int(*v_it);
        }
    } else if(d.getString() == "Z") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[2] > tip[2]) tip_point_index = int(*v_it);
        }
    } else if(d.getString() == "-Z") {
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            Point p = mesh.point(*v_it);
            Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
            if(p[2] < tip[2]) tip_point_index = int(*v_it);
        }
    } else {
        abort();
    }
    
    Point tip = mesh.point(Mesh::Vertex_index(tip_point_index));
    Eigen::Vector3d translation(-tip[0], -tip[1], -tip[2]);
    translateMesh(mesh, translation);

    return -1.0 * translation;
}

// Translates mesh by given 3D column vector
void MeshTools::translateMesh(Mesh& mesh, Eigen::Vector3d translation) {
    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        double x = p[0] + translation(0);
        double y = p[1] + translation(1);
        double z = p[2] + translation(2);
        mesh.point(*v_it) = Point(x, y, z);
    }

    return;
}

// Rotates mesh matrix from the axis_from to the axis_to
// Assumes that axis_from and axis_to are orthogonal
void MeshTools::rotateMeshMatrixAxisOnce(Eigen::MatrixXd& mesh_matrix, Direction axis_from, Direction axis_to) {
    int n_m = mesh_matrix.cols();
    int n_d = mesh_matrix.rows();

    assert(n_d % 3 == 0);
    assert(axis_from.isOrthogonalTo(axis_to));

    // Any single 90-degree rotation can be represented as one of X to Y, Y to X, Y to Z, Z to Y, Z to X, X to Z (without minus direction)
    Direction current_from = axis_from;
    Direction current_to = axis_to;
    Direction rotation_axis = current_from.outerProductDirection(current_to);
    while(current_from.getString()[0] == '-' || current_to.getString()[0] == '-') {
        current_from = current_from.clockwiseDirection(rotation_axis);
        current_to = current_to.clockwiseDirection(rotation_axis);
    }

    int axis_from_idx, axis_to_idx;
    if(current_from.getString() == "X") axis_from_idx = 0;
    else if(current_from.getString() == "Y") axis_from_idx = 1;
    else if(current_from.getString() == "Z") axis_from_idx = 2;
    else abort();
    if(current_to.getString() == "X") axis_to_idx = 0;
    else if(current_to.getString() == "Y") axis_to_idx = 1;
    else if(current_to.getString() == "Z") axis_to_idx = 2;
    else abort();

    for(int i = 0; i < n_d; i += 3) {
        Eigen::MatrixXd points = mesh_matrix.block(i, 0, 3, n_m);
        Eigen::MatrixXd new_points(3, n_m);
        for(int j = 0; j < 3; j++) {
            if(j == axis_from_idx) {
                new_points.row(j) = -1.0 * points.row(axis_to_idx);
            } else if(i == axis_to_idx) {
                new_points.row(j) = points.row(axis_from_idx);
            } else {
                new_points.row(j) = points.row(i);
            }
        }
        mesh_matrix.block(i, 0, 3, n_m) = new_points;
    }

    return;
}

// Rotates a mesh from the axis_from to the axis_to
// Assumes that axis_from and axis_to are orthogonal
void MeshTools::rotateMeshAxisOnce(Mesh& mesh, Direction axis_from, Direction axis_to) {
    assert(axis_from.isOrthogonalTo(axis_to));

    // Any single 90-degree rotation can be represented as one of X to Y, Y to X, Y to Z, Z to Y, Z to X, X to Z (without minus direction)
    Direction current_from = axis_from;
    Direction current_to = axis_to;
    Direction rotation_axis = current_from.outerProductDirection(current_to);
    while(current_from.getString()[0] == '-' || current_to.getString()[0] == '-') {
        current_from = current_from.clockwiseDirection(rotation_axis);
        current_to = current_to.clockwiseDirection(rotation_axis);
    }

    int axis_from_idx, axis_to_idx;
    if(current_from.getString() == "X") axis_from_idx = 0;
    else if(current_from.getString() == "Y") axis_from_idx = 1;
    else if(current_from.getString() == "Z") axis_from_idx = 2;
    else abort();
    if(current_to.getString() == "X") axis_to_idx = 0;
    else if(current_to.getString() == "Y") axis_to_idx = 1;
    else if(current_to.getString() == "Z") axis_to_idx = 2;
    else abort();

    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point p = mesh.point(*v_it);
        double new_p[3];
        for(int i = 0; i < 3; i++) {
            if(i == axis_from_idx) {
                new_p[i] = -p[axis_to_idx];
            } else if(i == axis_to_idx) {
                new_p[i] = p[axis_from_idx];
            } else {
                new_p[i] = p[i];
            }
        }
        mesh.point(*v_it) = Point(new_p[0], new_p[1], new_p[2]);
    }

    return;
}