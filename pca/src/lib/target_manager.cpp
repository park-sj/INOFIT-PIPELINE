#include "target_manager.hpp"

TargetManager::TargetManager() {}

// Deletes dynamically allocated Mesh objects on destruction
TargetManager::~TargetManager() {
    clear();
}

// Reads mesh files and stores meshes with their face and height directions
// Assumes that the face and height directions are orthogonal
void TargetManager::addMeshes(std::vector<std::string>& mesh_files, Direction face, Direction height) {
    for(int i = 0; i < mesh_files.size(); i++) {
        addMesh(mesh_files[i], face, height);
    }
    return;
}

// Stores meshes with their face and height directions
// Assumes that the face and height directions are orthogonal
void TargetManager::addMeshes(std::vector<Mesh*>& meshes, Direction face, Direction height) {
    for(int i = 0; i < meshes.size(); i++) {
        addMesh(*(meshes[i]), face, height);
    }
    return;
}

// Reads a mesh file and stores a mesh with its face and height direction
// Assumes that the face and height direction are orthogonal
void TargetManager::addMesh(std::string mesh_file, Direction face, Direction height) {
    assert(face.isOrthogonalTo(height));

    Mesh* m = new Mesh();
    std::ifstream mesh_ifstream(mesh_file);
    mesh_ifstream >> *m;
    mesh_ifstream.close();

    assert(m->is_empty() == false);

    meshes.push_back(m);
    initial_face_directions.push_back(face);
    initial_height_directions.push_back(height);
    current_face_directions.push_back(face);
    current_height_directions.push_back(height);

    return;
}

// Stores a mesh with its face and height direction
// Assumes that the face and height direction are orthogonal
void TargetManager::addMesh(Mesh& mesh, Direction face, Direction height) {
    assert(face.isOrthogonalTo(height));

    Mesh* m = new Mesh();
    *m = mesh;

    meshes.push_back(m);
    initial_face_directions.push_back(face);
    initial_height_directions.push_back(height);
    current_face_directions.push_back(face);
    current_height_directions.push_back(height);

    return;
}

// Aligns based on the position of nose tip
// At first, rotates template meshes to face Z and height Y directions,
// then translates template meshes so that the nose tips of template meshes are at the origin
// Assumes that meshes have ordinary noses
// Note that scaling is not included (maintaining the initial scale)
void TargetManager::noseTipAlign() {
    MeshTools mt;

    for(int i = 0; i < meshes.size(); i++) {
        mt.rotateMeshAxis(*(meshes[i]), current_face_directions[i], current_height_directions[i], Direction("Z"), Direction("Y"));
        current_face_directions[i] = Direction("Z");
        current_height_directions[i] = Direction("Y");
    }

    translations = mt.translateMeshesTipToOrigin(meshes, Direction("Z"));
    
    return;
}

// Translates and rotates given mesh to the initial position and direction of the original target mesh
void TargetManager::alignBackToTargetMesh(Mesh& mesh, int target_mesh_index) {
    MeshTools mt;

    mt.translateMesh(mesh, translations.col(target_mesh_index));

    mt.rotateMeshAxis(mesh, current_face_directions[target_mesh_index], current_height_directions[target_mesh_index],
        initial_face_directions[target_mesh_index], initial_height_directions[target_mesh_index]);

    return;
}

// Translates and rotates given mesh to the initial position and direction of the original target mesh,
// then saves the mesh as an off file
void TargetManager::saveMesh(Mesh& mesh, std::string file_name, bool do_remove_color) {
    if(do_remove_color) {
        RegionTools rt;
        rt.removeColorProperty(mesh);
    }

    std::ofstream output_ofstream(file_name);
    output_ofstream << mesh;
    output_ofstream.close();

    return;
}

// Returns a vector of mesh pointers by value
std::vector<TargetManager::Mesh*> TargetManager::getMeshes() {
    return meshes;
}

// Deletes dynamically allocated Mesh objects
void TargetManager::clear() {
    for(auto it = meshes.begin(); it != meshes.end(); ++it) {
        delete *it;
    }

    meshes.clear();
    initial_face_directions.clear();
    initial_height_directions.clear();
    current_face_directions.clear();
    current_height_directions.clear();
    return;
}