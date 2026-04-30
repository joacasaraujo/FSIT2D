#include "scene/scene.h"
#include "scene/output.h"

int main(){

    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------- Geometry ----------
    double D                 = 20;
    S.model_max_corner       = { 20*D, 5*D + 2 , 1.0 };
    const int H              = static_cast<int>(S.model_max_corner[1]) - 2;
    
    // ---------- Parameters ----------
	double uMax              = 0.1;
	double Re                = 5;
    double kinViscosity      = uMax * H / Re;
    S.dt_lbm                 = 1.0;
    S.lattice_spacing        = 1.0; 
    S.relaxation_time        = 3.0 * kinViscosity + 0.5;
    S.initial_density        = 1.0;
    S.collision_operator     = "BGK";
    
    // ------ Guo Forcing ------
    const double F_x = 8.0 * kinViscosity * uMax / (H * H);
    S.GUO_fluid_forcing[0] = F_x;
    
    // ---------- Grid ----------
    S.AddRectangularCanal();
    
    // Engines (Do not change order!) 
    S.engines.clear();
    S.engines.push_back(std::make_shared<FluidStreaming>());
    S.engines.push_back(std::make_shared<ApplyFluidForcing>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    
    // -------- Solver ---------
    double nSteps = 50000;
    std::string file_name = "Poiseuille_GUO: ";
    while (S.time < nSteps){
        if (S.iter % 100 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
    return 0;
}