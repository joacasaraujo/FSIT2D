#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> point;
typedef boost::geometry::model::segment<point> segment;
typedef boost::geometry::model::box<point> aabb;
typedef boost::geometry::model::polygon<point, true, true> polygon;