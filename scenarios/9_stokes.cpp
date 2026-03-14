#include "scene/Scene.h"
#include "scene/Output.h"
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <string>

// ---------------- Helpers ----------------
static void BuildDiskPolygon(std::shared_ptr<Body>& B, double R, double dtheta = 0.1) {
    // Safer than cloud.clear()
    B->shape->cloud.outer().clear();
    B->shape->cloud.inners().clear();

    for (double theta = 0.0; theta < 2.0 * M_PI; theta += dtheta) {
        double x = B->state->pos[0] + R * std::cos(theta);
        double y = B->state->pos[1] + R * std::sin(theta);
        boost::geometry::append(B->shape->cloud.outer(), point(x, y));
    }
    boost::geometry::correct(B->shape->cloud);
}

static double Mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
}

static double StdDev(const std::vector<double>& v, double m) {
    if (v.size() < 2) return 0.0;
    double acc = 0.0;
    for (double x : v) {
        double d = x - m;
        acc += d * d;
    }
    return std::sqrt(acc / (double)(v.size() - 1));
}

// ---------------- Main ----------------
int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------- Geometry (set Ly manually for each run) ----------
    const double R  = 6.0;
    const double D  = 2.0 * R;

    const double Lx = 100.0 * R;
    const double Ly = 40.0 * R;   // <-- mude aqui (20R, 40R, 80R, 160R, ...)

    S.model_min_corner = {0.0, 0.0, 0.0};
    S.model_max_corner = {Lx, Ly, 1.0};

    // ---------- Fluid & Particle Parameters ----------
    const double rho_f = 1.0;
    const double nu    = 0.10;
    const double mu    = rho_f * nu;
    const double rho_p = 1.20;
    const double delta_rho = rho_p - rho_f;

    // Stokes Parameters
    const double g = 5.555555555555557e-05;
    const double F_net = delta_rho * M_PI * R * R * g;  
    const double settling_velocity_3D = 2 * R * R * g * delta_rho / (9 * mu);
    const double settling_velocity_2D = R * g * delta_rho / (6 * mu);

    // ---------- LBM Parameters ----------
    S.dt_lbm           = 1.0;
    S.lattice_spacing  = 1.0;
    S.initial_density  = rho_f;
    S.initial_velocity = Vector3r::Zero();
    S.relaxation_time  = 3.0 * nu + 0.5;
    S.collision_operator = "MRT";
    S.dem_coupling       = "IMB";
    S.GUO_fluid_forcing  = Vector3r::Zero();
    double s8 = 2/(1 + 6*nu);
    S.SetRelaxationParamters(0, 1.4, 1.4, 0.75, 1.2, 1, 1.2, s8, s8);

    // ---------- DEM Parameters ----------
    const double mass = rho_p * M_PI * R * R;
    S.friction_angle   = 30.0;
    S.local_damping    = 0.0;
    S.normal_stiffness = 5.0e2;
    S.shear_stiffness  = 2.5e2;
    S.border_stiffness = 5.0e3;
    S.gravity = Vector3r(0.0, -g, 0.0);

    // DEM subcycling
    const double dt_crit = 0.2 * std::sqrt(mass / S.normal_stiffness);
    const int dem_subcycles = std::max(1, (int)(S.dt_lbm / dt_crit) + 1);
    S.dt_dem = S.dt_lbm / dem_subcycles;

    // ---------- Grid ----------
    S.AddRectangularCanal();
    for (auto& C : S.cells) {
        C->set_initial_condition();
    }

    // ---------- Particle ----------
    Vector3r pos0(Lx * 0.5, Ly * 0.75, 0.0);
    S.AddDisk(pos0, R, rho_p);
    auto& B = S.bodies[0];

    B->state->vel    = Vector3r::Zero();
    B->state->rotVel = 0.0;

    BuildDiskPolygon(B, R);

    // ---------- Engines ----------
    S.engines.clear();
    S.engines.push_back(std::make_shared<LatticeSearch>());
    S.engines.push_back(std::make_shared<ImbBoundary>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<FluidStreaming>());

    std::vector<std::shared_ptr<Engine>> dem_engines;
    dem_engines.push_back(std::make_shared<ContactResolution>());
    dem_engines.push_back(std::make_shared<InteractionLoop>());
    dem_engines.push_back(std::make_shared<BodyLoop>());
    dem_engines.push_back(std::make_shared<Integrator>());
    dem_engines.push_back(std::make_shared<UpdateContact>());

    // ---------- Output ----------
    std::ostringstream fname;
    fname << "stokes_2D_Ly_" << (int)Ly << "_R_" << (int)R << ".csv";
    std::ofstream csv(fname.str());
    csv << "time,y,vy,Re,Fy_hydro,Fy_net,F_net,res_force\n";

    // ---------- Robust terminal extraction ----------
    const int output_interval = 200;
    const int vtk_interval    = 200;  
    const double max_time     = 10000.0;

    const int N_tail = 30;           
    std::vector<double> vy_tail;
    vy_tail.reserve(N_tail);

    const double stop_y = 0.25 * Ly; 
    const double start_collect_time = 1000.0; 

    while (S.time < max_time) {

        BuildDiskPolygon(B, R);

        // LBM
        for (auto& E : S.engines) E->action();

        // DEM
        for (int sub = 0; sub < dem_subcycles; ++sub) {
            for (auto& E : dem_engines) E->action();
        }

        if (S.iter % output_interval == 0) {
            const double vy = -B->state->vel[1]; // Settling velocity
            const double Re = vy * D / nu;

            const double Fy_hydro = B->state->hydro_force[1];
            const double Fy_net   = B->state->force[1];

            double res_force = abs(Fy_hydro - F_net)/F_net;

            csv << S.time << ","
                << B->state->pos[1] << ","
                << vy << ","
                << Re << ","
                << Fy_hydro << ","
                << Fy_net << ","
                << F_net << ","
                << res_force << "\n";

            std::cout << "t=" << std::setw(7) << std::fixed << std::setprecision(0) << S.time
                      << "  y=" << std::setw(8) << std::setprecision(2) << B->state->pos[1]
                      << "  vy=" << std::scientific << std::setprecision(3) << vy
                      << "  vy_theory_2D=" << std::scientific << std::setprecision(3) << settling_velocity_2D
                      << "  vy_theory_3D=" << std::scientific << std::setprecision(3) << settling_velocity_3D
                      << "  Re=" << std::fixed << std::setprecision(4) << Re
                      << "  2d ratio=" << std::scientific << std::setprecision(3) << (vy / settling_velocity_2D)
                      << "  3d ratio=" << std::scientific << std::setprecision(3) << (vy / settling_velocity_3D)
                      << "  F_hydro_lbm=" << std::scientific << Fy_hydro
                      << "  F_net=" << std::scientific << F_net
                      << "  res_force=" << std::scientific << res_force
                      << std::endl;

            if (S.time >= start_collect_time) {
                if ((int)vy_tail.size() < N_tail) {
                    vy_tail.push_back(vy);
                } else {
                    vy_tail.erase(vy_tail.begin());
                    vy_tail.push_back(vy);
                }

                if ((int)vy_tail.size() == N_tail) {
                    const double m = Mean(vy_tail);
                    const double s = StdDev(vy_tail, m);
                    const double cv = (std::abs(m) > 1e-30) ? (s / std::abs(m)) : 0.0;

                    if (cv < 0.01) { // 1%
                        std::cout << "\n*** Terminal regime (robust) detected ***\n"
                                  << "Ut_mean = " << std::scientific << m
                                  << "  Ut_std = " << s
                                  << "  (std/mean=" << std::fixed << std::setprecision(2) << (100.0*cv) << "%)\n\n";
                        break;
                    }
                }
            }
        }

        if (S.iter % vtk_interval == 0) {
            O.ExportFluidVtk("Stokes2D_Fluid_"+std::to_string(int(Ly))+"_");
            O.ExportParticleVtk("Stokes2D_Particle_"+std::to_string(int(Ly))+"_");
        }

        if (B->state->pos[1] < stop_y) {
            std::cout << "\nParticle reached stop_y (avoid bottom effects). Stopping.\n\n";
            break;
        }

        S.time += S.dt_lbm;
        ++S.iter;
    }

    csv.close();

    // ---------- Final robust metrics ----------
    if (!vy_tail.empty()) {
        const double m = Mean(vy_tail);
        const double s = StdDev(vy_tail, m);
        const double cv = (std::abs(m) > 1e-30) ? (s / std::abs(m)) : 0.0;

        std::cout << "\n===== Robust 2D summary =====\n";
        std::cout << "Ly = " << Ly << "   (Ly/D=" << (Ly/D) << ", ln(Ly/D)=" << std::log(Ly/D) << ")\n";
        std::cout << "Ut_mean = " << std::scientific << m << "\n";
        std::cout << "Ut_std  = " << s << "  (std/mean=" << std::fixed << std::setprecision(2) << (100.0*cv) << "%)\n";
        std::cout << "F_net   = " << std::scientific << F_net << "\n";
        std::cout << "=============================\n\n";
    } else {
        std::cout << "\nNo tail samples collected (increase max_time or lower start_collect_time).\n";
    }

    return 0;
}
