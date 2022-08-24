#include <vector>
#include <utility>

// Defined to prevent the warning "The practice of declaring the Bind placeholders (_1, _2, ...)
// in the global namespace is deprecated." from the Boost library
#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Splitters.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/squared_distance_3.h>
#include <Eigen/Dense>

#ifndef MESH_CORRESPONDENCE_H
#define MESH_CORRESPONDENCE_H

class MeshCorrespondence {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

    typedef boost::property_map<Mesh, CGAL::vertex_point_t>::type VertexPointPmap;
    typedef CGAL::Search_traits_3<K> SearchTraitsBase;
    typedef CGAL::Search_traits_adapter<Mesh::Vertex_index, VertexPointPmap, SearchTraitsBase> SearchTraits;
    typedef CGAL::Midpoint_of_max_spread<SearchTraits> Splitter;
    typedef CGAL::Kd_tree<SearchTraits, Splitter> KD_Tree;
    typedef CGAL::Fuzzy_sphere<SearchTraits> FuzzySphere;

    typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> AABB_Primitive;
    typedef CGAL::AABB_traits<K, AABB_Primitive> AABB_Traits;
    typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;
    typedef AABB_Tree::Point_and_primitive_id Point_and_primitive_id;

public:
    MeshCorrespondence();
    void computeVV_KD(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, int> >& v);
    void computeVV_AABB(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, int> >& v);
    void computeVP_AABB(Mesh& source_mesh, Mesh& target_mesh, std::vector<std::pair<int, Point> >& v);
    Eigen::VectorXd computeVP_AABB(Eigen::VectorXd source_mesh_vector, Mesh& target_mesh);

    void build_AABB(Mesh& target_mesh);
    Eigen::VectorXd computeOnlyVP_AABB(Eigen::VectorXd source_mesh_vector);

private:
    AABB_Tree aabb_tree;

    constexpr static const double INITIAL_MIN_VALUE = 987654321.0;
    constexpr static const double KD_MAX_DISTANCE = 0.1;
};

#endif