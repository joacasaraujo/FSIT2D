#ifndef SHAPE_H
#define SHAPE_H

#include "maths/Math.h"

class Shape {
    public:
    void generate_particle(std::string _file_name);

    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> coordinates;
    std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> conectivities;
    double radius;
    double inertia_moment;
    polygon cloud;
};

#endif