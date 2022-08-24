#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
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
    double kd_time, aabb_time;
    std::vector<std::pair<int, int> > correspondence_kd, correspondence_aabb;
    MeshCorrespondence mc;

    // Finding correspondence of mesh vertices by k-d tree
    std::cout << "Finding correspondence via k-d tree" << std::endl;
    sc.start();
    mc.computeVV_KD(source_mesh, target_mesh, correspondence_kd);
    kd_time = sc.stop();

    // Finding correspondence of mesh vertices by AABB tree
    std::cout << "Finding correspondence via AABB tree" << std::endl;
    sc.start();
    mc.computeVV_AABB(source_mesh, target_mesh, correspondence_aabb);
    aabb_time = sc.stop();

    // Comparing the two methods, k-d tree and AABB tree
    for(int i = 0; i < source_mesh.num_vertices(); i++) {
        std::cout << correspondence_kd[i].first << " : "
            << correspondence_kd[i].second << " (kd), "
            << correspondence_aabb[i].second << " (aabb)" << std::endl;
    }

    std::cout << "k-d tree took " << kd_time << "s" << std::endl;
    std::cout << "AABB tree took " << aabb_time << "s" << std::endl;

    return EXIT_SUCCESS;
}