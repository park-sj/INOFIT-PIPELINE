#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include "mesh_correspondence.hpp"
#include "tools.hpp"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

int main(int argc, char* argv[]) {
    Mesh source_mesh, target_mesh;

    if (argc != 3) {
        std::cerr << "Usage: " + std::string(argv[0]) + " [source off file] [target off file]" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading " << std::string(argv[1]) << std::endl;
    std::ifstream source_ifstream(argv[1]);
    source_ifstream >> source_mesh;
    source_ifstream.close();

    if (source_mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[1]) << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Reading " << std::string(argv[2]) << std::endl;
    std::ifstream target_ifstream(argv[2]);
    target_ifstream >> target_mesh;
    target_ifstream.close();


    if (target_mesh.is_empty()) {
        std::cerr << "Cannot read " << std::string(argv[2]) << std::endl;
        return EXIT_FAILURE;
    }

    SimpleClock sc;
    double kd_time, aabb_time, aabb_time_vp;
    std::vector<std::pair<int, int> > correspondence_kd, correspondence_aabb;
    std::vector<std::pair<int, Point> > correspondence_aabb_vp;
    MeshCorrespondence mc;

    // Finding closest vertices by k-d tree
    std::cout << "Finding closest vertices via k-d tree" << std::endl;
    sc.start();
    mc.computeVV_KD(source_mesh, target_mesh, correspondence_kd);
    kd_time = sc.stop();

    // Finding closest vertices by AABB tree
    std::cout << "Finding closest vertices via AABB tree" << std::endl;
    sc.start();
    mc.computeVV_AABB(source_mesh, target_mesh, correspondence_aabb);
    aabb_time = sc.stop();

    // Finding closest points by AABB tree
    std::cout << "Finding closest points via AABB tree" << std::endl;
    sc.start();
    mc.computeVP_AABB(source_mesh, target_mesh, correspondence_aabb_vp);
    aabb_time_vp = sc.stop();

    // Comparing the three methods by squared distance
    for(int i = 0; i < source_mesh.num_vertices(); i++) {
        std::cout << correspondence_kd[i].first << " : ";
        if(correspondence_kd[i].second != -1) {
            std::cout << CGAL::squared_distance(source_mesh.point((Mesh::Vertex_index)correspondence_kd[i].first),
                    target_mesh.point((Mesh::Vertex_index)correspondence_kd[i].second)) << " (kd), ";
        } else {
            std::cout << "x (kd), ";
        }
        std::cout << CGAL::squared_distance(source_mesh.point((Mesh::Vertex_index)correspondence_aabb[i].first),
                target_mesh.point((Mesh::Vertex_index)correspondence_aabb[i].second)) << " (aabb), ";
        std::cout << CGAL::squared_distance(source_mesh.point((Mesh::Vertex_index)correspondence_aabb_vp[i].first),
                correspondence_aabb_vp[i].second) << " (aabb vertex-point)";
        std::cout << std::endl;
    }

    std::cout << "k-d tree (closest vertices) took " << kd_time << "s" << std::endl;
    std::cout << "AABB tree (closest vertices) took " << aabb_time << "s" << std::endl;
    std::cout << "AABB tree (closest points) took " << aabb_time_vp << "s" << std::endl;

    return EXIT_SUCCESS;
}