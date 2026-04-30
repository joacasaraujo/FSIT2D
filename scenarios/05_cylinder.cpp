#include "scene/scene.h"
#include "scene/output.h"

int main(){

    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // --------- Grid ----------
    std::string file_name = "Cylinder: ";
    double D = 20;                                                    
    S.model_max_corner = { 20*D, 5*D + 2 , 1.0 };
    const int H = static_cast<int>(S.model_max_corner[1]) - 2;
    
    // --------- LBM ----------- 
	double uMax              = 0.1;                                             
	double Re                = 500;
    double kinViscosity      = uMax * D / Re;
    S.dt_lbm                 = 1.0;
    S.lattice_spacing        = 1.0; 
    S.relaxation_time        = 3.0 * kinViscosity + 0.5;
    S.initial_density        = 1.0;
    S.initial_velocity       = Vector3r::Zero();
    S.collision_operator     = "MRT";
    S.dem_coupling           = "IMB";

    double s8 = 2/(1 + 6*kinViscosity);
    S.SetRelaxationParamters(0, 1.4, 1.4, 0.75, 1.2, 1, 1.2, s8, s8);
    
    // -------- Bodies ---------
    S.AddDisk(Vector3r(S.model_max_corner[0]/4 ,S.model_max_corner[1]/2, 0), D/2.0, 1.0);

    // -------- Walls ---------
    S.AddRectangularCanal();

    // ------ Zou and He BC + Poiseuille ------
    S.apply_velocity_bc = 4;
    for (int j = 1; j <= H; ++j) {
        double y  = j - 1.0;                   
        double ux = 4.0*uMax*(y/H - (y*y)/(H*H));
        int id_in  = S.CalculateCellId(0, j);
        S.cells[id_in ]->lattice->velocity_bc = Vector3r(ux, 0.0, 0.0);
    }

    // Engines (Do not change order!) 
    S.engines.clear();
    S.engines.push_back(std::make_shared<FluidCollision>());   
    S.engines.push_back(std::make_shared<LatticeSearch>());    
    S.engines.push_back(std::make_shared<ImbBoundary>());      
    S.engines.push_back(std::make_shared<FluidStreaming>());   
    S.engines.push_back(std::make_shared<ZouAndHeBC>());      

    // Solver
    while (S.time < 10000.0) {
        if (S.iter % 100 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
    return 0;
}