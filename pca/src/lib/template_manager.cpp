#include "template_manager.hpp"

TemplateManager::TemplateManager() {}

// Deletes dynamically allocated Mesh objects on destruction
TemplateManager::~TemplateManager() {
    clear();
}

// Reads mesh files and stores meshes with their face and height directions
// Assumes that the face and height directions are orthogonal
void TemplateManager::addMeshes(std::vector<std::string>& mesh_files, Direction face, Direction height) {
    for(int i = 0; i < mesh_files.size(); i++) {
        addMesh(mesh_files[i], face, height);
    }
    return;
}

// Stores meshes with their face and height directions
// Assumes that the face and height directions are orthogonal
void TemplateManager::addMeshes(std::vector<Mesh*>& meshes, Direction face, Direction height) {
    for(int i = 0; i < meshes.size(); i++) {
        addMesh(*(meshes[i]), face, height);
    }
    return;
}

// Reads a mesh file and stores a mesh with its face and height direction
// Assumes that the face and height direction are orthogonal
// While region-vertex mapping for given mesh is stored individually,
// it is assumed that all added meshes have the same topology and region-vertex mapping
void TemplateManager::addMesh(std::string mesh_file, Direction face, Direction height) {
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

    RegionTools rt;
    region_vertex_mappings.push_back(std::vector<std::set<int> >());
    rt.readRegionsFromColor(*m, region_vertex_mappings.back());

    return;
}

// Stores a mesh with its face and height direction
// Assumes that the face and height direction are orthogonal
// While region-vertex mapping for given mesh is stored individually,
// it is assumed that all added meshes have the same topology and region-vertex mapping
void TemplateManager::addMesh(Mesh& mesh, Direction face, Direction height) {
    assert(face.isOrthogonalTo(height));

    Mesh* m = new Mesh();
    *m = mesh;

    meshes.push_back(m);
    initial_face_directions.push_back(face);
    initial_height_directions.push_back(height);
    current_face_directions.push_back(face);
    current_height_directions.push_back(height);

    RegionTools rt;
    region_vertex_mappings.push_back(std::vector<std::set<int> >());
    rt.readRegionsFromColor(*m, region_vertex_mappings.back());

    return;
}

// Reads, removes, merges, and expands regions of the stored meshes
// Returns the number of resulting regions
// Assumes that all the topologies and regions of the stored meshes are the same
// Colors are updated to stored meshes to represent modified regions
// Note that region overlap is not reflected to the color (impossible to represent)
int TemplateManager::arrangeRegions(std::set<int>& regions_to_remove, std::vector<std::set<int> >& regions_to_merge, int additional_region_overlap) {
    RegionTools rt;

    // Region indexes are changed every time some regions are merged
    // So the regions indexes are tracked
    std::vector<int> region_index_tracking;
    for(int i = 0; i < region_vertex_mappings[0].size(); i++) {
        region_index_tracking.push_back(i);
    }

    // Merging specified regions
    // Newly merged regions are appropriately colored for each mesh
    for(int i = 0; i < regions_to_merge.size(); i++) {
        std::set<int> regions_to_merge_with_index_corrected;
        for(auto it = regions_to_merge[i].begin(); it != regions_to_merge[i].end(); ++it) {
            regions_to_merge_with_index_corrected.insert(region_index_tracking[*it]);
        }

        for(int j = 0; j < meshes.size(); j++) {
            rt.mergeRegionsAndSetColor(*(meshes[j]), region_vertex_mappings[j], regions_to_merge_with_index_corrected);
        }
        
        int merged_region_index = *(regions_to_merge_with_index_corrected.begin());
        for(int j = 0; j < region_index_tracking.size(); j++) {
            if(regions_to_merge_with_index_corrected.find(region_index_tracking[j]) != regions_to_merge_with_index_corrected.end()) {
                region_index_tracking[j] = merged_region_index;
            } else {
                int index_decrease = 0;
                for(auto it = ++(regions_to_merge_with_index_corrected.begin()); it != regions_to_merge_with_index_corrected.end(); ++it) {
                    if(region_index_tracking[j] > *it) index_decrease += 1;
                    else break;
                }
                region_index_tracking[j] -= index_decrease;
            }
        }
    }

    // Removing specified regions
    // Corresponding faces are removed for each mesh
    // Note that after region removal all the vertex indexes are changed
    for(int i = 0; i < meshes.size(); i++) {
        for(auto it = regions_to_remove.begin(); it != regions_to_remove.end(); ++it) {
            rt.removeFacesInRegion(*(meshes[i]), region_vertex_mappings[i], region_index_tracking[*it]);
        }
        rt.removeActually(*(meshes[i]), region_vertex_mappings[i]);
    }

    // Expanding regions to have larger overlaps between regions
    // Note that overlapping regions are not expressed in color
    for(int i = 0; i < meshes.size(); i++) {
        rt.widenAllRegions(*(meshes[i]), region_vertex_mappings[i], additional_region_overlap);
    }
/*
    for(int i = 0; i < region_vertex_mappings[0].size(); i++) {
        for(auto it = region_vertex_mappings[0][i].begin(); it != region_vertex_mappings[0][i].end(); ++it) {
            assert(*it < 5847);
        }
    }
*/
    return region_vertex_mappings[0].size();
}

// Aligns based on the position of nose tip
// At first, rotates template meshes to face Z and height Y directions,
// then translates template meshes so that the nose tips of template meshes are at the origin
// Assumes that meshes have ordinary noses
// Note that scaling is not included (maintaining the initial scale)
void TemplateManager::noseTipAlign() {
    MeshTools mt;

    for(int i = 0; i < meshes.size(); i++) {
        mt.rotateMeshAxis(*(meshes[i]), current_face_directions[i], current_height_directions[i], Direction("Z"), Direction("Y"));
        current_face_directions[i] = Direction("Z");
        current_height_directions[i] = Direction("Y");
    }

    mt.translateMeshesTipToOrigin(meshes, Direction("Z"));
    
    return;
}

// Returns a vector of mesh pointers by value
std::vector<TemplateManager::Mesh*> TemplateManager::getMeshes() {
    return meshes;
}

// Returns a pointer to the region vertex mapping vector of i-th mesh
std::vector<std::set<int> >* TemplateManager::getRegionVertexMapping(int i) {
    assert(i >= 0 && i < region_vertex_mappings.size());
    return &(region_vertex_mappings[i]);
}

// Deletes dynamically allocated Mesh objects
void TemplateManager::clear() {
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