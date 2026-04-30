#include "body.h"

bool Body::CheckInteraction(int _body_id) {
	for (auto& I : inter) {
		auto Il = I.lock();

		auto Ilb1 = Il->body1.lock();
		auto Ilb2 = Il->body2.lock();

		if ((Ilb1->id == id && Ilb2->id == _body_id) || (Ilb2->id == id && Ilb1->id == _body_id)) {
			return true;
		}
	}
	return false;
}