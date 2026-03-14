#ifndef MATERIAL_H
#define MATERIAL_H

class Material {
    public:
    
    int id = -1;
    double density = 0;
    double mass = 0;       // Computed from density and shape
    double inertia = 0;
};

class Fluid : public Material {
    public:
    
    double kinViscosity = 0;
};

#endif