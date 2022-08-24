#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#ifndef REGION_TOOLS_HPP
#define REGION_TOOLS_HPP

class RegionTools {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef K::Point_3 Point;
    typedef CGAL::Surface_mesh<Point> Mesh;
    typedef CGAL::Color Color;

public:
    RegionTools();

    void readRegionsFromColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping);
    void mergeRegions(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, std::set<int> regions_to_merge);
    void mergeRegionsAndSetColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, std::set<int> regions_to_merge);
    void widenRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, int hops);
    void widenAllRegions(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int hops);
    void removeFacesInRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region);
    void removeActually(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping);
    void setRegionColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Color color);
    Color getRegionColor(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region);
    void applyColor(Mesh& mesh, Mesh& color_reference);
    void clearColorToWhite(Mesh& mesh);
    void removeColorProperty(Mesh& mesh);

    void writeCorrectOffFileWithFaceColor(Mesh& mesh, std::string file_name);

private:
    int colorToInt(Color color);
    Color intToColor(int color);

    bool isVertexBorderOfRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Mesh::Vertex_index vi);
    bool isFaceInRegion(Mesh& mesh, std::vector<std::set<int> >& region_vertex_mapping, int region, Mesh::Face_index fi);

    void removeTemporaryFile(std::string filename);
    
};

#endif