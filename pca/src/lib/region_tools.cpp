#include "region_tools.hpp"

RegionTools::RegionTools() {}

// Constructs region-vertex mapping in a mesh by its face color
// Regions can be overlapped with each other (at the border)
void RegionTools::readRegionsFromColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping) {
    region_vertex_mapping.clear();

    // If the given mesh has no face color information, all vertices belong to a single region
    bool colorExists = mesh.property_map<Mesh::Face_index, Color>("f:color").second;
    if(!colorExists) {
        region_vertex_mapping.push_back(std::set<int>());
        for(int i = 0; i < mesh.number_of_vertices(); i++) {
            region_vertex_mapping[0].insert(i);
        }
        return;
    }

    // Iterating all faces in mesh
    Mesh::Face_range f_range = mesh.faces();
    Mesh::Face_range::iterator f_begin, f_end, f_it;
    f_begin = f_range.begin();
    f_end = f_range.end();

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;
    auto fc_begin = fcolors.begin();
    auto fc_it = fc_begin;

    std::map<int, int> color_region_mapping;
    int region_number = 0;

    for (f_it = f_begin; f_it != f_end; ++f_it) {
        int color_integer = colorToInt(*fc_it);
        // If a new color is encountered, new region is created
        if(color_region_mapping.find(color_integer) == color_region_mapping.end()) {
            color_region_mapping[color_integer] = region_number++;
            region_vertex_mapping.push_back(std::set<int>());
        }
        ++fc_it;

        Mesh::Halfedge_index hi = mesh.halfedge(*f_it);
        CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> > fv_range = mesh.vertices_around_face(hi);
        CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> >::iterator fv_begin, fv_end, fv_it;
        fv_begin = fv_range.begin();
        fv_end = fv_range.end();

        for (fv_it = fv_begin; fv_it != fv_end; ++fv_it) {
            region_vertex_mapping[color_region_mapping[color_integer]].insert((int)*fv_it);
        }
    }

    return;
}

// Merges regions specified in regions_to_merge set
// Note that the region number of each region is reassigned after merge
// Colors of the mesh are not changed in this function
void RegionTools::mergeRegions(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, std::set<int> regions_to_merge) {
    // Checking if all regions to merge are in the range of current regions
    for(auto it = regions_to_merge.begin(); it != regions_to_merge.end(); ++it) {
        assert(*it >= 0 && *it < region_vertex_mapping.size());
    }

    // Nothing to do
    if(regions_to_merge.size() < 2) return;

    // Moving vertices from higher regions to the lowest region
    int lowest_region = *(regions_to_merge.begin());
    for(auto rit = regions_to_merge.rbegin(); rit != regions_to_merge.rend(); ++rit) {
        if(*rit != lowest_region) {
            region_vertex_mapping[lowest_region].insert(region_vertex_mapping[*rit].begin(), region_vertex_mapping[*rit].end());
            region_vertex_mapping.erase(region_vertex_mapping.begin() + *rit);
        }
    }

    return;
}

// Merges regions specified in regions_to_merge set
// Note that the region number of each region is reassigned after merge
// Color of the merged region is newly set to express that the vertices are in the same region
void RegionTools::mergeRegionsAndSetColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, std::set<int> regions_to_merge) {
    // Nothing to do
    if(regions_to_merge.size() < 2) return;
    
    // Keeping the color of the first region in regions_to_merge
    int first_region = *(regions_to_merge.begin());
    Color color = getRegionColor(mesh, region_vertex_mapping, first_region);

    // Merging regions to update region_vertex_mapping
    mergeRegions(mesh, region_vertex_mapping, regions_to_merge);

    // Actually updating the color of the merged region in mesh
    setRegionColor(mesh, region_vertex_mapping, first_region, color);

    return;
}

// Widens a region by the amount of hops (to make the overlapped region wider)
// Mesh colors are not changed in this function
void RegionTools::widenRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, int hops) {
    assert(region >= 0 && region < region_vertex_mapping.size());

    for(int h = 0; h < hops; h++) {
        Mesh::Vertex_range v_range = mesh.vertices();
        Mesh::Vertex_range::iterator v_begin, v_end, v_it;
        v_begin = v_range.begin();
        v_end = v_range.end();

        // Finding all the vertices on the region border
        std::vector<Mesh::Vertex_index> v;
        for (v_it = v_begin; v_it != v_end; ++v_it) {
            if(isVertexBorderOfRegion(mesh, region_vertex_mapping, region, *v_it)) {
                v.push_back(*v_it);
            }
        }

        // Vertices around a border vertex are newly included to the region
        for(int i = 0; i < v.size(); i++) {
            Mesh::Halfedge_index hi = mesh.halfedge(v[i]);
            CGAL::Iterator_range<CGAL::Vertex_around_target_iterator<Mesh> > vv_range = mesh.vertices_around_target(hi);
            CGAL::Iterator_range<CGAL::Vertex_around_target_iterator<Mesh> >::iterator vv_begin, vv_end, vv_it;
            vv_begin = vv_range.begin();
            vv_end = vv_range.end();

            for (vv_it = vv_begin; vv_it != vv_end; ++vv_it) {
                region_vertex_mapping[region].insert((int)*vv_it);
            }
        }
    }

    return;
}

// Widens all regions by the amount of hops (to make the overlapped region wider)
// Mesh colors are not changed in this function
void RegionTools::widenAllRegions(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int hops) {
    for(int r = 0; r < region_vertex_mapping.size(); r++) {
        widenRegion(mesh, region_vertex_mapping, r, hops);
    }
    return;
}

// Removes all faces in a region
// Note that this function does not update region_vertex_mapping
// Instead, this function removes corresponding faces from the mesh
// Based on the CGAL Surface_mesh implementation, here removing means that corresponding faces are internally marked as "removed",
// not changing the current indexes of faces, edges, halfedges, and vertices
// To actually remove the "removed" faces and its vertices from the mesh and off file, call removeActually() after this
void RegionTools::removeFacesInRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region) {
    assert(region >= 0 && region < region_vertex_mapping.size());

    // Iterating all faces in mesh
    Mesh::Face_range f_range = mesh.faces();
    Mesh::Face_range::iterator f_begin, f_end, f_it;
    f_begin = f_range.begin();
    f_end = f_range.end();

    for (f_it = f_begin; f_it != f_end; ++f_it) {
        if(isFaceInRegion(mesh, region_vertex_mapping, region, *f_it)) {
            mesh.remove_face(*f_it);
        }
    }

    return;
}

// Acually removes faces marked as "removed" from the mesh and off file
// Assumes that some faces are marked as "removed" by previous call of removeFacesInRegion()
// Also removes corresponding vertices whose incident faces are all removed (as well as corresponding halfedges and edges)
// This function changes the indexes of faces, edges, halfedges, and vertices,
// and also properly updates region_vertex_mapping
void RegionTools::removeActually(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping) {
    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    // Iterating vertices to find vertices to be removed
    for(v_it = v_begin; v_it != v_end; ++v_it) {
        Mesh::Halfedge_index hi = mesh.halfedge(*v_it);
        CGAL::Iterator_range<CGAL::Face_around_target_iterator<Mesh> > vf_range = mesh.faces_around_target(hi);
        CGAL::Iterator_range<CGAL::Face_around_target_iterator<Mesh> >::iterator vf_begin, vf_end, vf_it;
        vf_begin = vf_range.begin();
        vf_end = vf_range.end();

        bool should_be_removed = true;
        for(vf_it = vf_begin; vf_it != vf_end; ++vf_it) {
            if(vf_it->is_valid() && !mesh.is_removed(*vf_it)) {
                should_be_removed = false;
                break;
            }
        }

        if(should_be_removed) {
            mesh.remove_vertex(*v_it);
        }
    }

    // Actually removing marked faces and vertices
    // Note that corresponding halfedges and edges are not removed here (by CGAL implementation)
    mesh.collect_garbage();

    // Writing to a temporary off file and read it again
    // to update corresponding edges and halfedges correctly
    std::string temp_file = "mesh/temp/temp.off";
    writeCorrectOffFileWithFaceColor(mesh, temp_file);
    Mesh new_mesh;
    std::ifstream mesh_ifstream(temp_file);
    mesh_ifstream >> new_mesh;
    mesh_ifstream.close();
    removeTemporaryFile(temp_file);

    // Finally updating region_vertex_mapping and mesh appropriately
    readRegionsFromColor(new_mesh, region_vertex_mapping);
    mesh = new_mesh;

    return;
}

// Sets all face colors in a region to a specific color
void RegionTools::setRegionColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Color color) {
    assert(region >= 0 && region < region_vertex_mapping.size());

    // Iterating all faces in mesh
    Mesh::Face_range f_range = mesh.faces();
    Mesh::Face_range::iterator f_begin, f_end, f_it;
    f_begin = f_range.begin();
    f_end = f_range.end();

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;

    for (f_it = f_begin; f_it != f_end; ++f_it) {
        if(isFaceInRegion(mesh, region_vertex_mapping, region, *f_it)) {
            fcolors[*f_it] = color;
        }
    }

    return;
}

// Gets the color of the region
CGAL::Color RegionTools::getRegionColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region) {
    assert(region >= 0 && region < region_vertex_mapping.size());
    assert(region_vertex_mapping[region].size() > 0);

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;

    Mesh::Face_range f_range = mesh.faces();
    Mesh::Face_range::iterator f_begin, f_end, f_it;
    f_begin = f_range.begin();
    f_end = f_range.end();

    Color color;
    for (f_it = f_begin; f_it != f_end; ++f_it) {
        if(isFaceInRegion(mesh, region_vertex_mapping, region, *f_it)) {
            color = fcolors[*f_it];
            break;
        }
    }

    return color;
}

// Copies all the colors from the reference mesh
// Assumes that reference mesh has the same topology with given mesh
void RegionTools::applyColor(Mesh& mesh, Mesh& color_reference) {
    bool colorExists = color_reference.property_map<Mesh::Face_index, Color>("f:color").second;
    if(!colorExists)
    {
        removeColorProperty(mesh);
        return;
    }

    colorExists = mesh.property_map<Mesh::Face_index, Color>("f:color").second;
    if(!colorExists)
    {
        mesh.add_property_map<Mesh::Face_index, Color>("f:color", Color(255, 255, 255));
    }

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;
    auto fc_begin = fcolors.begin();
    auto fc_end = fcolors.end();
    auto fc_it = fc_begin;

    Mesh::Property_map<Mesh::Face_index, Color> fcolors_ref = color_reference.property_map<Mesh::Face_index, Color>("f:color").first;
    auto fc_ref_it = fcolors_ref.begin();

    for(fc_it = fc_begin; fc_it != fc_end; ++fc_it) {
        *fc_it = *fc_ref_it;
        ++fc_ref_it;
    }

    return;
}

// Sets all face colors to RGB 255 255 255 (white)
void RegionTools::clearColorToWhite(Mesh& mesh) {
    // If the given mesh has no face color information, there is nothing to do
    bool colorExists = mesh.property_map<Mesh::Face_index, Color>("f:color").second;
    if(!colorExists) return;

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;
    auto fc_begin = fcolors.begin();
    auto fc_end = fcolors.end();
    auto fc_it = fc_begin;

    for (fc_it = fc_begin; fc_it != fc_end; ++fc_it) {
        *fc_it = Color(255, 255, 255);
    }

    return;
}

// Removes color property from the mesh
void RegionTools::removeColorProperty(Mesh& mesh) {
    // If the given mesh has no face color information, there is nothing to do
    bool colorExists = mesh.property_map<Mesh::Face_index, Color>("f:color").second;
    if(!colorExists) return;

    Mesh::Property_map<Mesh::Face_index, Color> fcolors = mesh.property_map<Mesh::Face_index, Color>("f:color").first;
    mesh.remove_property_map(fcolors);

    return;
}

// Writes face-colored "OFF" file
// Note that default CGAL implementation outputs OFF file with its first line "COFF" when a mesh has face colors,
// which is problematic to read the OFF file again from CGAL
void RegionTools::writeCorrectOffFileWithFaceColor(Mesh& mesh, std::string file_name) {
    std::ofstream mesh_ofstream(file_name);
    mesh_ofstream << mesh;
    mesh_ofstream.seekp(0);
    mesh_ofstream << "OFF ";
    mesh_ofstream.close();
    return;
}

// Converts an RGB color to an one-to-one mapping integer
int RegionTools::colorToInt(Color color) {
    return int(color[0]) + 256*int(color[1]) + 256*256*int(color[2]);
}

// Converts a color integer to an one-to-one mapping RGB color
CGAL::Color RegionTools::intToColor(int color) {
    assert(color >= 0 && color < 256*256*256);

    return Color(color%256, (color/256)%256, color/(256*256));
}

// Checks if a given vertex vi is on the border of the region
bool RegionTools::isVertexBorderOfRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Mesh::Vertex_index vi) {
    assert(region >= 0 && region < region_vertex_mapping.size());
    
    // If a vertex vi is in the region but at least one of its neighbor vertices is not in the same region,
    // then the vertex vi is on the border of the region
    bool ret = false;
    if(region_vertex_mapping[region].find((int)vi) != region_vertex_mapping[region].end()) {
        Mesh::Halfedge_index hi = mesh.halfedge(vi);
        CGAL::Iterator_range<CGAL::Vertex_around_target_iterator<Mesh> > vv_range = mesh.vertices_around_target(hi);
        CGAL::Iterator_range<CGAL::Vertex_around_target_iterator<Mesh> >::iterator vv_begin, vv_end, vv_it;
        vv_begin = vv_range.begin();
        vv_end = vv_range.end();

        for (vv_it = vv_begin; vv_it != vv_end; ++vv_it) {
            if(region_vertex_mapping[region].find((int)*vv_it) == region_vertex_mapping[region].end()) {
                ret = true;
                break;
            }
        }
    }

    return ret;
}

// Checks if a given face fi is in the region (that is, all vertices of a given face are in the region)
bool RegionTools::isFaceInRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Mesh::Face_index fi) {
    assert(region >= 0 && region < region_vertex_mapping.size());

    Mesh::Halfedge_index hi = mesh.halfedge(fi);
    CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> > fv_range = mesh.vertices_around_face(hi);
    CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> >::iterator fv_begin, fv_end, fv_it;
    fv_begin = fv_range.begin();
    fv_end = fv_range.end();

    bool ret = true;
    for (fv_it = fv_begin; fv_it != fv_end; ++fv_it) {
        if(region_vertex_mapping[region].find(*fv_it) == region_vertex_mapping[region].end()) {
            ret = false;
            break;
        }
    }

    return ret;
}

// Removes a file
// Used to remove temporary off file
void RegionTools::removeTemporaryFile(std::string filename) {
    std::remove(filename.c_str());
    return;
}