#ifndef BODY_H
#define BODY_H

#include "maths/Math.h"
#include "property/Material.h"
#include "property/State.h"
#include "property/Shape.h"
#include "dem/Interaction.h"
#include "lbm/Cell.h"

class Interaction;

class Body {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Body(int _id, std::shared_ptr<Shape> _shape, std::shared_ptr<State> _state, std::shared_ptr<Material> _material) : id(_id), shape(_shape), state(_state), material(_material){}

	bool CheckInteraction(int _body_id);

	//Body variables:
	int    id;

	//Smart pointers:
	std::vector<std::weak_ptr<Interaction>> inter;
	std::vector<std::shared_ptr<Cell>> fluid_interactions;
	std::shared_ptr<State> state;
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;

	Vector3r blockedDOFs = { 1, 1, 1 };
	double blockedMDOFs  = 1;

	~Body() {};
};

#endif //BODY_H