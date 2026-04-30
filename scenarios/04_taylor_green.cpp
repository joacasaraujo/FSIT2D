#include "scene/scene.h"
#include "scene/output.h"
#include <cmath>

static inline double sq(double x){ return x*x; }

int main(){
    Timer T;
    Scene& S = S.get_Scene();
    Output O;

    // ---------- Geometry ----------
    S.model_max_corner   = { 128, 128, 1.0 };
    const int H = static_cast<int>(S.model_max_corner[1]);
    
    // ---------- Parameters ---------- 
    const double u0      = 0.05;
    double Re            = 5;
    double kinViscosity  = u0 * H / Re;
    S.dt_lbm             = 1.0;
    S.latticeSpeed       = 1.0;
    S.relaxation_time    = 3.0 * kinViscosity + 0.5;
    S.initial_density    = 1.0;
    S.collision_operator = "BGK";

    // -------- Grid ---------
    S.AddRectangularCanal();
    for (auto& Cptr : S.cells) {
        Cptr->is_wall  = false;  
        Cptr->is_solid = false;
    }  

    // Engines (Do not change order!) 
    S.engines.clear();
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<FluidStreaming>());

    // -------- Taylor-Green Initialization --------
    const double Lx = S.model_max_corner[0], Ly = S.model_max_corner[1];
    const double kx = 2.0*M_PI / Lx;
    const double ky = 2.0*M_PI / Ly;

    for (auto& Cptr : S.cells){
        auto& C = *Cptr;
        const double x = C.state->pos[0];
        const double y = C.state->pos[1];

        const double ux0 = +u0 * std::sin(kx*x) * std::cos(ky*y);
        const double uy0 = -u0 * std::cos(kx*x) * std::sin(ky*y);
        S.initial_velocity[0] = ux0;
        S.initial_velocity[1] = uy0;

        for (int k = 0; k < C.lattice->number_of_nodes; ++k){
            C.lattice->f[k] = C.CalculateEqFunction(S.initial_density, S.initial_velocity, k);
        }
        C.CalculateDensityAndVelocity();
    }

    // -------- Solver ---------
    double nSteps = 10000;
    std::string file_name = "Taylo-Green: ";
    while (S.time < nSteps){
        if (S.iter % 100 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
    return 0;
}