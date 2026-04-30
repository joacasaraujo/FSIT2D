#ifndef CELL_H
#define CELL_H

#include "lbm/lattice.h"
#include "property/state.h"
#include "property/material.h"

class Scene;

class Cell {
    public:
    Cell(int _id, std::shared_ptr<Lattice> _lattice, std::shared_ptr<Material> _material, std::shared_ptr<State> _state) : 
    id(_id), lattice(_lattice), material(_material), state(_state), is_solid(false), is_wall(false){}
    ~Cell(){}
    
    int  id;
    bool is_solid;
    bool is_wall;
    double solid_fraction = 0.0;
    polygon grid;

    std::shared_ptr<Lattice> lattice;
    std::shared_ptr<Material> material;
    std::shared_ptr<State> state;

    double CalculateEqFunction(double _density, Vector3r _velocity, int k);  //Calculate Equilibrium Function
    void CalculateDensityAndVelocity();
    void set_neighbor_node();
    void set_initial_condition();
    void set_mrt_parameters();
};
#endif