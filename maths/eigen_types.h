#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

// Tipos numéricos base
// Você pode mudar aqui o tipo global (float/double) se quiser no futuro
using Real = double;

// --------------------------------------------
// Vetores comuns
// --------------------------------------------
template<typename Scalar>
using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
using Vector3r = Vector3<Real>;

template<typename Scalar>
using Vector9 = Eigen::Matrix<Scalar, 9, 1>;
using Vector9r = Vector9<Real>;

// --------------------------------------------
// Matrizes comuns
// --------------------------------------------
template<typename Scalar>
using Matrix9 = Eigen::Matrix<Scalar, 9, 9>;
using Matrix9r = Matrix9<Real>;