#ifndef SCENE_H
#define SCENE_H

// #include "maths/Math.h"
#include "lbm/cell.h"
#include "dem/body.h"
#include "dem/interaction.h"
#include "engine/engine.h"
#include "maths/math.h"

class Body;
class Interaction;

class Scene {
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    static Scene& get_Scene();
    int  CalculateCellId(int i, int j) {return i + j * model_max_corner[0];}
    void AddRectangularCanal();
    void addParticle(std::string _fileName, double _density, Vector3r _velocity);
    void AddDisk(Vector3r _center, double _radius, double _density);
    void MoveToNextTimeStep();
    void SetRelaxationParamters(double _s1, double _s2, double _s3, double _s4, double _s5, double _s6, double _s7, double _s8, double _s9);
    // Generates a random disk packing; optional params control count/densities
    void addDiskPacking(double _Rmin, double _Rmax, double _wall_gap, double _inlet_buf, double _outlet_buf, double _non_overlap,
                        int _max_particles = 0, int _wall_particles = 0, double _rho_p = 1.0,  double _Rw_min = 0.0, double _Rw_max = 0.0);

    std::vector<std::shared_ptr<Cell>> cells;
    std::vector<std::shared_ptr<Body>> bodies;
    std::vector<std::shared_ptr<Interaction>> interactions;
    std::vector<std::shared_ptr<Engine>> engines;

    //Simulation Parameters:
    int      iter           = 0;
    double   time           = 0.0;
    Vector3r model_min_corner = Vector3r::Zero();
    Vector3r model_max_corner = Vector3r::Zero();
    Vector3r gravity        = {0.0, -9.81, 0.0};

    //LBM Parameters:
    int apply_velocity_bc = -1;
    int apply_density_bc = -1;
    double lattice_spacing      = 1.0;
    double dt_lbm   = 1.0;
    double relaxation_time     = 1.0;
    double initial_density = 1.0;
    double latticeSpeed = 1.0;
    std::string collision_operator = "MRT";
    std::string dem_coupling = "None";
    Vector3r initial_velocity = Vector3r::Zero();
    Vector3r GUO_fluid_forcing = Vector3r::Zero();
    Vector9r relaxation_parameters = Vector9r::Zero();

    //DEM Parameters:
    double dt_dem            = 1.0;
    double friction_angle    = 30;
    double local_damping     = 0.0;
    double factor_of_safety  = 0.3;
    double normal_stiffness  = 1e6;
    double shear_stiffness   = 0.5e6;
    double border_stiffness  = 1e6;
    int nIter = 0;

};

#endif