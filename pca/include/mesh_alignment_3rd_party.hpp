#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <iterator>
#include <cstdio>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/property_map.h>

#include <CGAL/OpenGR/compute_registration_transformation.h>
#include <CGAL/OpenGR/register_point_sets.h>
#include <CGAL/pointmatcher/compute_registration_transformation.h>
#include <CGAL/pointmatcher/register_point_sets.h>

#ifndef MESH_ALIGNMENT_3RD_PARTY_HPP
#define MESH_ALIGNMENT_3RD_PARTY_HPP

class MeshAlignment3rdParty {
private:
    // Should be float to work with OpenGR
    typedef CGAL::Simple_cartesian<float> K;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef CGAL::Surface_mesh<Point> Mesh;

    typedef boost::graph_traits<Mesh>::vertex_descriptor VertexDescriptor;

    typedef std::pair<Point, Vector> Pwn;
    typedef CGAL::First_of_pair_property_map<Pwn> PointMap;
    typedef CGAL::Second_of_pair_property_map<Pwn> NormalMap;

public:
    MeshAlignment3rdParty();
    void coarseAlignment(Mesh& source_mesh, Mesh& target_mesh);
    bool fineAlignment(Mesh& source_mesh, Mesh& target_mesh);
    bool generalAlignment(Mesh& source_mesh, Mesh& target_mesh);

private:
    std::string computeNormalsAsXYZ(Mesh& mesh);
    void removeTemporaryFile(std::string filename);

};

#endif