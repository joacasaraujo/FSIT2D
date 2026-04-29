#include "scene/Scene.h"
#include "scene/Output.h"

// Write in english, use english variables and comments

int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------- Geometry ----------
    double D                 = 20;
    S.model_max_corner       = { 20*D, 5*D + 2 , 1.0 };
    const int H              = static_cast<int>(S.model_max_corner[1]) - 2;

    // ---------- Parameters ----------
    double Re                = 5;
    double uMax              = 0.1;
    double kinViscosity      = uMax * H / Re;
    S.dt_lbm                 = 1.0;
    S.lattice_spacing        = 1.0; 
    S.relaxation_time        = 3.0 * kinViscosity + 0.5;
    S.initial_density        = 1.0;
    S.collision_operator     = "BGK";
    
    // -------- Bodies ---------
    // add vertical line of bodies close to the inlet to apply horizontal force
    for (int i = 0; i < H; ++i) {
        double x = 0.1 * D;
        double y = i + 0.5;
        S.AddDisk(Vector3r(x, y, 0), 0.5 * D, 1.0);
    }

    for (auto& B : S.bodies) {
        B->state->force = Vector3r(1, 0.0, 0.0);
    }

    // ---------- Grid ----------
    S.AddRectangularCanal();

    // -------- Engines ---------
    S.engines.push_back(std::make_shared<LatticeSearch>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<ImbBoundary>());
    S.engines.push_back(std::make_shared<FluidStreaming>());
    S.engines.push_back(std::make_shared<ZouAndHeBC>());

    std::vector<std::shared_ptr<Engine>> dem_engines;
    dem_engines.push_back(std::make_shared<ContactResolution>());
    dem_engines.push_back(std::make_shared<InteractionLoop>());
    dem_engines.push_back(std::make_shared<BodyLoop>());
    dem_engines.push_back(std::make_shared<Integrator>());
    dem_engines.push_back(std::make_shared<UpdateContact>());

    // -------- Solver ---------
    double nSteps = 50000;
    std::string file_name = "Soil_Column: ";
    while (S.time < nSteps){
        if (S.iter % 100 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
    return 0;
}