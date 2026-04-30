#include "scene/scene.h"
#include "scene/output.h"

int main(){
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------- Geometry ----------
    const int N = 256;                       
    S.model_max_corner = { double(N), double(N), 1.0 };

    // ---------- Parameters ----------
    const double u_lid        = 0.05;                
    const double Re           = 100.0;
    const double kinViscosity = u_lid * N / Re;    
    S.relaxation_time         = 3.0 * kinViscosity + 0.5;
    S.initial_density         = 1.0;
    S.collision_operator      = "MRT"; 
    
    if (S.collision_operator == "MRT"){
        //Razzaghian (2012) para MRT
        std::cout <<"Here" << std::endl;
        double s8 = 2/(1 + 6*kinViscosity);
        S.SetRelaxationParamters(0, 1.4, 1.4, 0.75, 1.2, 1, 1.2, s8, s8);
    } 

    // ---------- Grid ----------
    S.AddRectangularCanal();                 

    // Adjusting walls of cavity
    for (auto& C : S.cells){
        C->is_wall  = false;
        C->is_solid = false;
    }
    for (int i=0;i<N;++i){
        S.cells[S.CalculateCellId(i,0)]->is_wall   = true; // bottom   
    }
    for (int j=0;j<N;++j){
        S.cells[S.CalculateCellId(0,j)]->is_wall   = true; // left
        S.cells[S.CalculateCellId(N-1,j)]->is_wall = true; // right
    }

    // ---------- Engines ----------
    S.engines.clear();
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<FluidStreaming>());
    S.engines.push_back(std::make_shared<ZouAndHeBC>());
    
    // ------ Zou And He BC  -------
    S.apply_velocity_bc = 1;   
    for (int i = 1; i < N-1; ++i){
        int id = S.CalculateCellId(i, N-1);
        S.cells[id]->lattice->velocity_bc = Vector3r(u_lid, 0.0, 0.0);
    }

    // -------- Solver ---------
    double nSteps = 50000;
    std::string file_name = S.collision_operator + "_Lid-Driven: ";
    while (S.time < nSteps){
        if (S.iter % 500 == 0) {
            O.DisplaySimulationInfo();
            O.ExportFluidVtk(file_name);
        }
        S.MoveToNextTimeStep();
    }
        // print kinViscosity and relaxation_time
    std::cout << "kinViscosity: " << kinViscosity << std::endl;
    std::cout << "relaxation_time: " << S.relaxation_time << std::endl;
    return 0;
}
