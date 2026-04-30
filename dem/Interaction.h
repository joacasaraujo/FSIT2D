#ifndef INTERACTION_H
#define INTERACTION_H

#define _USE_MATH_DEFINES
#include "maths/math.h"
#include "dem/body.h"
#include "scene/scene.h"

class Body;

class Interaction {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Interaction(std::weak_ptr<Body> _body1, std::weak_ptr<Body> _body2) : body1(_body1), body2(_body2) {};

	//Methods:
	bool CheckContact();
	void CalculateUnitVectorandContact();
	void CalculateForceAndShearIncrements();
	void ApplyFrictionLaw();

	//Smart pointers:
	std::weak_ptr<Body> body1;
	std::weak_ptr<Body> body2;

	//Variables:
	Vector3r  unit_normal_vector  = Vector3r::Zero();
	Vector3r  unit_shear_vector   = Vector3r::Zero();
	Vector3r  contact     = Vector3r::Zero();
	double normal_force = 0.0;
	double shear_force  = 0.0;

};

#endif