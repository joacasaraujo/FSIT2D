#ifndef STATE_H
#define STATE_H

#include "maths/math.h"

class State {
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Vector3r prevPos   = Vector3r::Zero();
    Vector3r pos       = Vector3r::Zero();
    Vector3r vel       = Vector3r::Zero();
    Vector3r force     = Vector3r::Zero();
    Vector3r permForce = Vector3r::Zero();
    Vector3r hydro_force = Vector3r::Zero();  // Hydrodynamic force from fluid

    double moment = 0.0;
    double rotVel = 0.0;
    double rot    = 0.0; 
    double hydro_torque = 0.0;  // Hydrodynamic torque from fluid
    
};

#endif