#include "scene/scene.h"    
#include "scene/output.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>

static inline double sqr(double x){ return x*x; }

int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------------- Domain ----------------
    double d = 20.0;
    double Rref = 0.5 * d;
    const int Nx = 20*d;
    const int Ny = 10*d + 2;
    S.model_min_corner = {0.0, 0.0, 0.0};
    S.model_max_corner = {double(Nx), double(Ny), 1.0};
    
    // ---------------- Fluid ----------------
    const double nu    = 0.10;
    const double rho_f = 1.0;
    S.dt_lbm           = 1.0;
    S.lattice_spacing  = 1.0;
    S.initial_density  = 1.0;
    S.initial_velocity = Vector3r::Zero();
    S.relaxation_time  = 3.0 * nu + 0.5;   // tau
    S.collision_operator = "BGK";
    S.dem_coupling       = "IMB";
    double s8 = 2/(1 + 6*nu);
    S.SetRelaxationParamters(0, 1.4, 1.4, 0.75, 1.2, 1, 1.2, s8, s8);

    // Guo forcing (interpreted as acceleration g_x)
    const double g_x = 2e-6;
    S.GUO_fluid_forcing = Vector3r(g_x, 0.0, 0.0);
    S.gravity           = Vector3r::Zero();

    // ---------------- DEM (disks fixed) ----------------
    S.friction_angle   = 30.0;
    S.local_damping    = 0.2;
    S.normal_stiffness = 1.0e3;
    S.shear_stiffness  = 5.0e2;
    S.border_stiffness = 1.0e4;

    // ---------------- Grid and walls ----------------
    S.AddRectangularCanal();

    // Ensure all cells start as fluid (except walls you set)
    for (auto& C : S.cells) {
        C->is_wall  = false;
        C->is_solid = false;
    }

    // Top and bottom no-slip walls (LBM wall cells)
    for (int i = 0; i < Nx; ++i) {
        S.cells[S.CalculateCellId(i, 0    )]->is_wall = true;
        S.cells[S.CalculateCellId(i, Ny-1)]->is_wall = true;
    }

    //---------------- Particle Packing ----------------
    std::string fileName;
    int numberOfParticles = 39;
    for(int i = 0; i < numberOfParticles; ++i){
        fileName = "particles/particles_displaced_d_20/" + std::to_string(i) + ".txt"; //Path para adicionar particulas
        S.addParticle(fileName, 1, Vector3r(0,0,0));
    }

    // Block motion (fixed geometry)
    for (auto& B : S.bodies) {
        B->blockedDOFs  = Vector3r::Zero();
        B->blockedMDOFs = 0.0;
    }

    // "Geometric" porosity from disk areas over total domain (rough, not IMB-grid exact)
    double area_solid = 0.0;
    for (auto& B : S.bodies) area_solid += M_PI * sqr(B->shape->radius);
    const double eps_geom = 1.0 - area_solid / ((Nx - 2*Rref) * (Ny - 2*Rref));

    // Carman–Kozeny-like estimate (approx)
    const double k_th_geom = (eps_geom*eps_geom*eps_geom * d*d) / (180.0 * sqr(1.0 - eps_geom));

    // ---------------- Engines ----------------
    S.engines.clear();
    S.engines.push_back(std::make_shared<LatticeSearch>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<ImbBoundary>());
    S.engines.push_back(std::make_shared<ApplyFluidForcing>());
    S.engines.push_back(std::make_shared<FluidStreaming>());    

    std::vector<std::shared_ptr<Engine>> dem_engines;
    dem_engines.push_back(std::make_shared<ContactResolution>());
    dem_engines.push_back(std::make_shared<InteractionLoop>());
    dem_engines.push_back(std::make_shared<BodyLoop>());
    dem_engines.push_back(std::make_shared<Integrator>());
    dem_engines.push_back(std::make_shared<UpdateContact>());

    // DEM dt (bodies fixed, but keep stable)
    const double mass_ref = 2.5 * M_PI * sqr(Rref);
    const double dt_crit  = 0.1 * std::sqrt(mass_ref / S.normal_stiffness);
    const int dem_subcycles = std::max(1, int(S.dt_lbm / dt_crit) + 1);
    S.dt_dem = S.dt_lbm / dem_subcycles;

    // ---------------- Measurement window ----------------
    // Use same region where particles are placed (buffered zone)
    const int i0 = std::max(0,   int(std::floor(106)));
    const int i1 = std::min(Nx-1,int(std::floor(293)));
    const int j0 = std::max(0,   -int(std::floor(Rref)));
    const int j1 = std::min(Ny-1,int(std::floor(200)));

    auto compute_window_stats = [&]() {
        // Compute:
        // - U_pore: mean u_x over FLUID cells (not wall, not solid)
        // - eps_region: (#fluid cells)/(#total non-wall cells) in window
        double u_sum = 0.0;
        int n_fluid = 0;
        int n_total = 0;

        for (int j = j0; j <= j1; ++j){
            for (int i = i0; i <= i1; ++i){
                auto& C = *S.cells[S.CalculateCellId(i,j)];
                if (C.is_wall) continue;          // exclude walls from total area
                ++n_total;
                if (C.is_solid) continue;         // solid from IMB
                u_sum += C.state->vel[0];
                ++n_fluid;
            }
        }

        const double U_pore = (n_fluid > 0) ? (u_sum / n_fluid) : 0.0;
        const double eps_region = (n_total > 0) ? (double(n_fluid)/double(n_total)) : 0.0;

        return std::make_tuple(U_pore, eps_region, n_fluid, n_total);
    };

    // ---------------- Output ----------------
    std::string file_name = "Real_Particle_Packing.csv";
    std::ofstream csv(file_name);
    csv << "time,k,k_ergun,k_kozeny,C_eff,Porosity\n";

    const double max_time = 10000.0;
    const int output_interval = 100;
    const int vtk_interval    = 500;

    std::cout << "Permeability scenario (Darcy) with fixed disks\n";
    std::cout << "Window: i=[" << i0 << "," << i1 << "], j=[" << j0 << "," << j1 << "]\n";
    std::cout << "eps_geom=" << eps_geom << "  k_th_geom=" << k_th_geom << "\n";
    std::cout << "Guo forcing g_x=" << g_x << "  nu=" << nu << "  tau=" << S.relaxation_time << "\n";

    // ---------------- Time loop ----------------
    while (S.time < max_time) {

        // LBM step
        for (auto& E : S.engines) E->action();

        // DEM subcycles (bodies fixed; still update loops)
        for (int sub = 0; sub < dem_subcycles; ++sub)
            for (auto& E : dem_engines) E->action();

        // Output
        if (S.iter % output_interval == 0) {
            const double mu = rho_f * nu;

            auto [U_pore, eps_region, n_fluid, n_total] = compute_window_stats();

            //Kozeny-Carman equation
            double porosity = eps_region;
            const double k_kozeny = (porosity*porosity*porosity * d*d) / (180.0 * sqr(1.0 - porosity));
            const double k_ergun = (porosity*porosity*porosity * d*d) / (150 * sqr(1.0 - porosity));

            // Superficial (Darcy) velocity:
            const double U_super = eps_region * U_pore;

            
            // Darcy permeability (treat GUO forcing as acceleration g_x):
            // dp/dx = -rho*g_x  =>  U_super = (k/mu)*(rho*g_x)
            const double k = (g_x != 0.0) ? (mu * U_super / (rho_f * g_x)) : 0.0;
            
            // Kozeny constante efetivo inferido a partir do k medido
            const double C_eff = (porosity*porosity*porosity * d*d) / (k * sqr(1.0 - porosity));

            // KC “calibrado” que bate com o LBM (usando o C_eff do instante)
            const double k_kozeny_cal = (porosity*porosity*porosity * d*d) / (C_eff * sqr(1.0 - porosity));

            // Reynolds (both definitions, for transparency)
            const double Re_pore  = (nu != 0.0) ? (U_pore  * d / nu) : 0.0;
            const double Re_super = (nu != 0.0) ? (U_super * d / nu) : 0.0;

            csv << S.time << ","
                << k << ","
                << k_ergun << ","
                << k_kozeny << ","
                << C_eff << ","
                << eps_region << "\n";

            std::cout << "t=" << std::setw(7) << std::fixed << std::setprecision(0) << S.time
                      << "  U_pore="  << std::scientific << std::setprecision(3) << U_pore
                      << "  U_super=" << U_super
                      << "  k="       << k
                      << "  k_kozeny=" << k_kozeny
                      << "  k_kozeny_cal=" << k_kozeny_cal
                      << "  C_eff=" << C_eff
                      << "  eps_reg=" << std::fixed << std::setprecision(4) << eps_region
                      << "  porosity=" << porosity
                      << "  Re_pore=" << std::scientific << std::setprecision(3) << Re_pore
                      << std::endl;
        }

        if (S.iter % vtk_interval == 0) {
            O.ExportFluidVtk("(RPP)Permeability_Fluid_");
            O.ExportParticleVtk("(RPP)Permeability_Particles_");
        }

        S.time += S.dt_lbm;
        ++S.iter;
    }

    csv.close();
    std::cout << "Finished. Results: permeability_results.csv\n";
    return 0;
}
