#ifndef LATTICE_H
#define LATTICE_H

#include "maths/Math.h"

class Lattice {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        std::vector<int>            neighbor_node;
        double                      sound_speed       = 0;
        double                      density_bc        = 0;
        Vector3r                    velocity_bc       = Vector3r::Zero();
        Matrix9r                    m                 = Matrix9r::Zero();                              
        Matrix9r                    m_inv             = Matrix9r::Zero(); 
        Matrix9r                    relaxation_matrix = Matrix9r::Zero();                           
        Vector9r                    f_tmp             = Vector9r::Zero();
        Vector9r                    f                 = Vector9r::Zero();
        const int                   number_of_nodes   = 9;
        const std::vector<double>   node_weight       = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
        const std::vector<int>      opposite_node     = {0, 3, 4, 1, 2, 7, 8, 5, 6};
        const std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> discrete_velocity = {{0,0,0}, {1, 0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1, 1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0}};
};

#endif