#include "interaction.h"

bool Interaction::CheckContact() {
	auto b1 = body1.lock();
	auto b2 = body2.lock();
	if ((b2->state->pos - b1->state->pos).norm() < b1->shape->radius + b2->shape->radius) {
		return true;
	}
	return false;
}

void Interaction::CalculateUnitVectorandContact() {
	auto b1 = body1.lock();
	auto b2 = body2.lock();

	unit_normal_vector = (b2->state->pos - b1->state->pos) / (b2->state->pos - b1->state->pos).norm();
	unit_shear_vector = { unit_normal_vector[1], -unit_normal_vector[0], 0.0 };

	// Contact point on surface of body1
	contact = b1->state->pos + unit_normal_vector * b1->shape->radius;
}

void Interaction::CalculateForceAndShearIncrements() {
	const Scene& S = S.get_Scene();
	auto b1 = body1.lock();
	auto b2 = body2.lock();

	// Calculate overlap
	double distance = (b2->state->pos - b1->state->pos).norm();
	double overlap = b1->shape->radius + b2->shape->radius - distance;

	// Relative velocity (linear + rotational contributions)
	Vector3r relVel = (b1->state->vel - b2->state->vel) 
	                - (b1->state->rotVel * b1->shape->radius - b2->state->rotVel * b2->shape->radius) * unit_shear_vector;

	double shear_increment = (relVel.dot(unit_shear_vector)) * S.dt_dem;

	// Normal force from current overlap
	normal_force = S.normal_stiffness * overlap;
	// Shear force accumulated incrementally
	shear_force  += S.shear_stiffness  * shear_increment;
	ASSERT(S.normal_stiffness > 0 && S.shear_stiffness > 0);
}

void Interaction::ApplyFrictionLaw() {
	const Scene& S = S.get_Scene();
	double maxShearForce = normal_force * tan(S.friction_angle * M_PI / 180.0);
	if (abs(shear_force) > maxShearForce) {
		if (shear_force > 0)  shear_force =  maxShearForce;
		if (shear_force < 0)  shear_force = -maxShearForce;
	}
}