#include "scene/Scene.h"
#include "scene/Output.h"

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
    
    // ---------- Grid ----------
    S.AddRectangularCanal();
    
    // ------ Zou And He BC  -------
    S.apply_velocity_bc = 4;
    for (int j = 1; j <= H; ++j) {
        double y  = j - 1.0;                   
        double ux = 4.0*uMax*(y/H - (y*y)/(H*H));
        int id_in  = S.CalculateCellId(0, j);
        S.cells[id_in]->lattice->velocity_bc = Vector3r(ux, 0.0, 0.0);
    }
    
    // Engines (Do not change order!) 
    S.engines.clear();
    S.engines.push_back(std::make_shared<ZouAndHeBC>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<FluidStreaming>());
    
    // -------- Solver ---------
    double nSteps = 50000;
    std::string file_name = "Poiseuille_ZOU: ";
    while (S.time < nSteps){
        if (S.iter % 100 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
    return 0;
}