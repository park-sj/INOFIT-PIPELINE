#include "mesh_distortion.hpp"

MeshDistortion::MeshDistortion() {}

void MeshDistortion::distort(Mesh& mesh, double distortion_rate) {
    // Computing min, max, sum, mean, std_dev for each coordinate of vertices
    double min[3] = { INITIAL_MIN_VALUE, INITIAL_MIN_VALUE, INITIAL_MIN_VALUE };
    double max[3] = { INITIAL_MAX_VALUE, INITIAL_MAX_VALUE, INITIAL_MAX_VALUE };
    double sum[3] = { 0.0, 0.0, 0.0 };
    double mean[3] = { 0.0, 0.0, 0.0 };
    double std_dev[3] = { 0.0, 0.0, 0.0 };

    Mesh::Vertex_range v_range = mesh.vertices();
    Mesh::Vertex_range::iterator v_begin, v_end, v_it;
    v_begin = v_range.begin();
    v_end = v_range.end();
    int v_n = mesh.number_of_vertices();

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double a;
        for (int i = 0; i < 3; i++) {
            a = mesh.point(*v_it)[i];
            min[i] = a < min[i] ? a : min[i];
            max[i] = a > max[i] ? a : max[i];
            sum[i] += a;
        }
    }

    for (int i = 0; i < 3; i++) {
        mean[i] = sum[i] / v_n;
    }

    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double a;
        for (int i = 0; i < 3; i++) {
            a = mesh.point(*v_it)[i];
            std_dev[i] += (a - mean[i]) * (a - mean[i]);
        }
    }

    for (int i = 0; i < 3; i++) {
        std_dev[i] = sqrt(std_dev[i]/(v_n-1.0));
    }

    // Converting the coordinates into z-score
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double a[3];
        for (int i = 0; i < 3; i++) {
            a[i] = (mesh.point(*v_it)[i] - mean[i]) / std_dev[i];
        }
        mesh.point(*v_it) = Point(a[0], a[1], a[2]);
    }

    // Distorting the points by small amount
    RandomNumber random;
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double a[3];
        for (int i = 0; i < 3; i++) {
            a[i] = mesh.point(*v_it)[i];
            a[i] = a[i] + (random.get(200) - 100.0) * 0.01 * distortion_rate;
        }
        mesh.point(*v_it) = Point(a[0], a[1], a[2]);
    }

    // Converting the z-score back to the original scale
    for (v_it = v_begin; v_it != v_end; ++v_it) {
        double a[3];
        for (int i = 0; i < 3; i++) {
            a[i] = mesh.point(*v_it)[i] * std_dev[i] + mean[i];
        }
        mesh.point(*v_it) = Point(a[0], a[1], a[2]);
    }

    return;
}