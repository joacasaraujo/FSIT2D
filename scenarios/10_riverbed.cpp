#include "scene/scene.h"    
#include "scene/output.h"
#include <random>
#include <fstream>
#include <iomanip>

// Riverbed-like test:
// - Poiseuille flow via Zou & He BC (velocity inlet, density outlet)
// - Polydisperse mobile particles forming a bed
// - IMB LBM-DEM coupling

static double urand(std::mt19937& gen, double a, double b){
    std::uniform_real_distribution<double> dist(a,b);
    return dist(gen);
}

int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ------------------------------------------------------------
    // Domain
    // ------------------------------------------------------------
    const int Nx = 300;
    const int Ny = 120;
    S.model_min_corner = {0.0, 0.0, 0.0};
    S.model_max_corner = {double(Nx), double(Ny), 1.0};

    // ------------------------------------------------------------
    // Fluid parameters
    // ------------------------------------------------------------
    const double rho_f = 1.0;
    const double nu    = 0.06;              // kinematic viscosity
    S.dt_lbm           = 1.0;
    S.lattice_spacing  = 1.0;
    S.initial_density  = rho_f;
    S.initial_velocity = Vector3r::Zero();
    S.relaxation_time  = 3.0 * nu + 0.5;    // tau ~ 0.68
    S.collision_operator = "BGK";
    S.dem_coupling       = "IMB";
    double s_nu = 1.0 / S.relaxation_time;
    S.SetRelaxationParamters(0, 1.4, 1.4, 0, 1.2, 0, 1.2, s_nu, s_nu);

    // No body force; flow driven by inlet BC
    S.GUO_fluid_forcing = Vector3r::Zero();
    S.gravity           = Vector3r(0.0, -5.0e-6, 0.0); // only DEM uses it

    // ------------------------------------------------------------
    // DEM parameters
    // ------------------------------------------------------------
    S.friction_angle   = 30.0;
    S.local_damping    = 0.2;
    S.normal_stiffness = 1.0e3;
    S.shear_stiffness  = 5.0e2;
    S.border_stiffness = 1.0e4;

    // ------------------------------------------------------------
    // Grid and walls
    // ------------------------------------------------------------
    S.AddRectangularCanal();
    for (auto& C : S.cells){
        C->is_wall  = false;
        C->is_solid = false;
    }
    // Top/bottom walls
    for (int i = 0; i < Nx; ++i){
        S.cells[S.CalculateCellId(i, 0    )]->is_wall = true;
        S.cells[S.CalculateCellId(i, Ny-1)]->is_wall = true;
    }

    // ------------------------------------------------------------
    // Polydisperse bed generation (bottom 30% of height)
    // ------------------------------------------------------------
    std::mt19937 gen(42);
    double Rmin = 2.0;
    double Rmax = 6.0;
    double bed_top = 0.3 * Ny;
    int max_particles = 120;

    for (int n = 0; n < max_particles; ++n){
        double R = urand(gen, Rmin, Rmax);
        // rejection sampling for non-overlap
        bool placed = false;
        for (int trial = 0; trial < 200 && !placed; ++trial){
            double x = urand(gen, R+1, Nx - R - 1);
            double y = urand(gen, R+1, bed_top - R - 1);
            Vector3r pos(x, y, 0.0);
            // check overlap
            bool ok = true;
            for (auto& B : S.bodies){
                double dist = (pos - B->state->pos).head<2>().norm();
                if (dist < (R + B->shape->radius)*1.05){
                    ok = false; break;
                }
            }
            if (ok){
                double rho_p = 2.5; // sand-like
                S.AddDisk(pos, R, rho_p);
                placed = true;
            }
        }
    }

    // Initialize velocities of bodies
    for (auto& B : S.bodies){
        B->state->vel = Vector3r::Zero();
        B->state->rotVel = 0.0;
    }

    // ------------------------------------------------------------
    // Zou & He Poiseuille: inlet velocity profile, outlet density
    // ------------------------------------------------------------
    S.apply_velocity_bc = 4; // left inlet
    S.apply_density_bc  = 2; // right outlet
    double umax = 0.02;      // peak inlet velocity (reduced for stability)
    int H = Ny - 2;
    for (int j = 1; j <= Ny-2; ++j){
        double y = j - 1.0;
        double ux = 4.0 * umax * (y / H) * (1.0 - y / H);
        int id_in = S.CalculateCellId(0, j);
        S.cells[id_in]->lattice->velocity_bc = Vector3r(ux, 0.0, 0.0);
    }
    for (int j = 1; j <= Ny-2; ++j){
        int id_out = S.CalculateCellId(Nx-1, j);
        S.cells[id_out]->lattice->density_bc = S.initial_density;
    }

    // ------------------------------------------------------------
    // Engines
    // ------------------------------------------------------------
    S.engines.clear();
    // Recommended order: collision -> IMB -> streaming -> BC
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

    // DEM time step conservative
    double mass_ref = 2.5 * M_PI * Rmax * Rmax;
    double dt_crit = 0.1 * std::sqrt(mass_ref / S.normal_stiffness);
    int dem_subcycles = std::max(1, (int)(S.dt_lbm / dt_crit) + 1);
    S.dt_dem = S.dt_lbm / dem_subcycles;

    // ------------------------------------------------------------
    // Output
    // ------------------------------------------------------------
    std::ofstream csv("riverbed_results.csv");
    csv << "time,mean_y,Favg,Fmax,ux_center\n";

    double max_time = 80000.0;
    int output_interval = 500;
    int vtk_interval = 1000;

    std::cout << "Riverbed polydisperse scenario (Zou-He Poiseuille)\n";
    std::cout << "Particles: " << S.bodies.size() << "  R in [" << Rmin << "," << Rmax << "]\n";
    std::cout << "umax=" << umax << "  nu=" << nu << "  tau=" << S.relaxation_time << "\n";

    while (S.time < max_time){

        // Update IMB geometry
        for (auto& Bp : S.bodies){
            auto& B = *Bp;
            double R = B.shape->radius;
            B.shape->cloud.clear();
            for (double theta = 0.0; theta < 2.0 * M_PI; theta += 0.1) {
                double x = B.state->pos[0] + R * std::cos(theta);
                double y = B.state->pos[1] + R * std::sin(theta);
                boost::geometry::append(B.shape->cloud.outer(), point(x, y));
            }
            boost::geometry::correct(B.shape->cloud);
        }

        for (auto& E : S.engines) E->action();

        for (int sub = 0; sub < dem_subcycles; ++sub){
            for (auto& E : dem_engines) E->action();
        }

        if (S.iter % output_interval == 0){
            // Bed height (mean y of all particles)
            double sumy = 0.0, Fsum = 0.0, Fmax = 0.0;
            for (auto& B : S.bodies){
                sumy += B->state->pos[1];
                double Fm = B->state->hydro_force.norm();
                Fsum += Fm;
                Fmax = std::max(Fmax, Fm);
            }
            double mean_y = sumy / S.bodies.size();
            double Favg   = Fsum / S.bodies.size();

            int cx = Nx/2, cy = Ny/2;
            double ux_c = S.cells[S.CalculateCellId(cx, cy)]->state->vel[0];

            csv << S.time << "," << mean_y << "," << Favg << "," << Fmax << "," << ux_c << "\n";

            std::cout << "t=" << std::setw(7) << std::fixed << std::setprecision(0) << S.time
                      << "  mean_y=" << std::setprecision(2) << mean_y
                      << "  Favg="   << std::scientific << std::setprecision(3) << Favg
                      << "  Fmax="   << Fmax
                      << "  ux_c="   << ux_c
                      << std::endl;
        }

        if (S.iter % vtk_interval == 0){
            O.ExportFluidVtk("Riverbed_Fluid_");
            O.ExportParticleVtk("Riverbed_Particles_");
        }

        S.time += S.dt_lbm;
        ++S.iter;
    }

    csv.close();
    std::cout << "Finished. Results: riverbed_results.csv, VTKs: Riverbed_*" << std::endl;
    return 0;
}

