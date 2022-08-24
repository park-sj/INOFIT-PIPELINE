#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "tools.hpp"

#ifndef MESH_DISTORTION_H
#define MESH_DISTORTION_H

class MeshDistortion {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;

public:
    MeshDistortion();
    void distort(Mesh& mesh, double distortion_rate);

private:
    constexpr static const double INITIAL_MIN_VALUE = 987654321.0;
    constexpr static const double INITIAL_MAX_VALUE = -987654321.0;
};

#endif