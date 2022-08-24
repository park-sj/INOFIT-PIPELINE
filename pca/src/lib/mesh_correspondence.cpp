#include "mesh_correspondence.hpp"

MeshCorrespondence::MeshCorrespondence() {}

// Finds vertex to vertex correspondence via k-d tree
// If any target vertices are not close enough to a source vertex, corresponding vertex index is set to -1
void MeshCorrespondence::computeVV_KD(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, int> >& v) {
    // Vertex-Point Property Map
    VertexPointPmap vppmap = get(CGAL::vertex_point, target_mesh);

    // Building a k-d tree
    KD_Tree tree(target_mesh.vertices().begin(), target_mesh.vertices().end(), Splitter(), SearchTraits(vppmap));

    // Finding the nearest target point for each source mesh vertices
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    std::vector<Mesh::Vertex_index> nearest_points;
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        std::back_insert_iterator<std::vector<Mesh::Vertex_index> > back_it(nearest_points);
        FuzzySphere fs(source_mesh.point(*v_it), KD_MAX_DISTANCE, 0.0, SearchTraits(vppmap));

        tree.search(back_it, fs);

        double min_distance = INITIAL_MIN_VALUE;
        int min_id = -1;
        for (auto it = nearest_points.begin(); it != nearest_points.end(); ++it) {
            double distance = CGAL::squared_distance(source_mesh.point(*v_it), target_mesh.point(*it));
            if(distance < min_distance) {
                min_distance = distance;
                min_id = int(*it);
            }
        }
        v.push_back(std::make_pair(int(*v_it), min_id));

        nearest_points.clear();
    }

    return;
}

// Finds vertex to vertex correspondence via AABB tree
void MeshCorrespondence::computeVV_AABB(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, int> >& v) {
    // Building a AABB tree
    AABB_Tree tree(target_mesh.faces().begin(), target_mesh.faces().end(), target_mesh);

    // Building a secondary search structure (internal k-d tree) of AABB tree
    tree.accelerate_distance_queries();

    // Finding the nearest target point for each source mesh vertices
    // NOTE: The nearest vertex is not always one of vertices of the nearest triangle
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point_and_primitive_id pp = tree.closest_point_and_primitive(source_mesh.point(*v_it));
        Point closest_point = pp.first;
        Mesh::Face_index fi = pp.second;

        Mesh::Halfedge_index hi = target_mesh.halfedge(fi);
        CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> > fv_range = target_mesh.vertices_around_face(hi);
        CGAL::Iterator_range<CGAL::Vertex_around_face_iterator<Mesh> >::iterator fv_begin, fv_end, fv_it;
        fv_begin = fv_range.begin();
        fv_end = fv_range.end();

        double min_distance = INITIAL_MIN_VALUE;
        int min_id = -1;
        for (fv_it = fv_begin; fv_it != fv_end; ++fv_it) {
            double distance = CGAL::squared_distance(closest_point, target_mesh.point(*fv_it));
            if(distance < min_distance) {
                min_distance = distance;
                min_id = int(*fv_it);
            }
        }

        v.push_back(std::make_pair(int(*v_it), min_id));
    }

    return;
}

// Finds vertex to point correspondence via AABB tree
void MeshCorrespondence::computeVP_AABB(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, Point> >& v) {
    // Building a AABB tree
    AABB_Tree tree(target_mesh.faces().begin(), target_mesh.faces().end(), target_mesh);

    // Building a secondary search structure (internal k-d tree) of AABB tree
    tree.accelerate_distance_queries();

    // Finding the nearest target point for each source mesh vertices
    // Target point can not only be a vertex but also any point on face
    Mesh::Vertex_range v_range = source_mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        Point closest_point = tree.closest_point(source_mesh.point(*v_it));

        v.push_back(std::make_pair(int(*v_it), closest_point));
    }

    return;
}

// Finds vertex to point correspondence via AABB tree
// Returns the closest points in the form of column vector
Eigen::VectorXd MeshCorrespondence::computeVP_AABB(Eigen::VectorXd source_mesh_vector, Mesh& target_mesh) {
    int n_d = source_mesh_vector.size();
    Eigen::VectorXd closest_points_vector(n_d);

    assert(n_d % 3 == 0);

    // Building a AABB tree
    AABB_Tree tree(target_mesh.faces().begin(), target_mesh.faces().end(), target_mesh);

    // Building a secondary search structure (internal k-d tree) of AABB tree
    tree.accelerate_distance_queries();

    // Finding the nearest target point for each source mesh points
    // Target point can not only be a vertex but also any point on face
    for(int i = 0; i < n_d; i += 3) {
        double x = source_mesh_vector(i);
        double y = source_mesh_vector(i+1);
        double z = source_mesh_vector(i+2);

        Point closest_point = tree.closest_point(Point(x, y, z));

        closest_points_vector(i) = closest_point[0];
        closest_points_vector(i+1) = closest_point[1];
        closest_points_vector(i+2) = closest_point[2];
    }

    return closest_points_vector;
}

// Builds an AABB tree with target_mesh to later accelerate queries without redundant build operations
void MeshCorrespondence::build_AABB(Mesh& target_mesh) {
    // Building a AABB tree
    aabb_tree.insert(target_mesh.faces().begin(), target_mesh.faces().end(), target_mesh);
    aabb_tree.build(target_mesh);

    // Building a secondary search structure (internal k-d tree) of AABB tree
    aabb_tree.accelerate_distance_queries();

    return;
}

// Finds vertex to point correspondence via AABB tree
// Returns the closest points in the form of column vector
// Assumes that internal aabb_tree is already built via "build_AABB" function
Eigen::VectorXd MeshCorrespondence::computeOnlyVP_AABB(Eigen::VectorXd source_mesh_vector) {
    int n_d = source_mesh_vector.size();
    Eigen::VectorXd closest_points_vector(n_d);

    assert(n_d % 3 == 0);

    // Finding the nearest target point for each source mesh points
    // Target point can not only be a vertex but also any point on face
    for(int i = 0; i < n_d; i += 3) {
        double x = source_mesh_vector(i);
        double y = source_mesh_vector(i+1);
        double z = source_mesh_vector(i+2);

        Point closest_point = aabb_tree.closest_point(Point(x, y, z));

        closest_points_vector(i) = closest_point[0];
        closest_points_vector(i+1) = closest_point[1];
        closest_points_vector(i+2) = closest_point[2];
    }

    return closest_points_vector;
}