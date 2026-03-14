/**
 * Drained triaxial test scenario (2D DEM)
 *
 * Simulates a strain-controlled drained triaxial test on a granular sample:
 * - Confining stress (sigma_3) from fixed lateral boundaries
 * - Axial strain imposed by moving the top boundary downward at constant rate
 * - No fluid coupling (dry / drained: no pore pressure)
 *
 * Output: axial strain, deviator stress q = sigma_1 - sigma_3, volumetric strain,
 *         and stress–strain curve (CSV + console).
 */

#include "scene/Scene.h"
#include "scene/Output.h"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

static inline double sqr(double x) { return x * x; }

int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ---------------- Domain (triaxial sample: 2D slice) ----------------
    const double Lx = 80.0;   // sample width  [lu]
    const double Ly = 160.0;  // sample height [lu] (2:1 aspect ratio)
    const int Nx = static_cast<int>(Lx);
    const int Ny = static_cast<int>(Ly);

    S.model_min_corner = {0.0, 0.0, 0.0};
    S.model_max_corner = {Lx, Ly, 1.0};

    // ---------------- LBM (inactive for dry test; grid still required) ----------------
    S.dt_lbm          = 1.0;
    S.lattice_spacing = 1.0;
    S.initial_density = 1.0;
    S.initial_velocity = Vector3r::Zero();
    S.relaxation_time = 1.0;
    S.collision_operator = "BGK";
    S.dem_coupling    = "None";  // dry / drained: no fluid coupling
    S.GUO_fluid_forcing = Vector3r::Zero();
    S.gravity         = Vector3r::Zero();

    // ---------------- DEM (granular contact) ----------------
    S.friction_angle   = 30.0;
    S.local_damping    = 0.15;
    S.normal_stiffness = 1.0e4;
    S.shear_stiffness  = 5.0e3;
    S.border_stiffness = 1.0e5;  // stiff boundaries for clear stress transmission

    // ---------------- Grid and walls ----------------
    S.AddRectangularCanal();
    for (auto& C : S.cells) {
        C->is_wall  = false;
        C->is_solid = false;
    }
    // All four sides as no-slip (for consistency; DEM border does the work)
    for (int i = 0; i < Nx; ++i) {
        S.cells[S.CalculateCellId(i, 0)]->is_wall     = true;
        S.cells[S.CalculateCellId(i, Ny - 1)]->is_wall = true;
    }
    for (int j = 0; j < Ny; ++j) {
        S.cells[S.CalculateCellId(0, j)]->is_wall      = true;
        S.cells[S.CalculateCellId(Nx - 1, j)]->is_wall = true;
    }
    for (auto& C : S.cells)
        C->set_initial_condition();

    // ---------------- Particle packing (soil grains inside sample) ----------------
    const double Rmin = 1.2;
    const double Rmax = 2.2;
    const double wall_gap = 0.5;
    const double inlet_buf = 2.0;
    const double outlet_buf = 2.0;
    const int max_particles = 400;
    const int wall_particles = 0;  // no extra wall particles for triaxial
    const double rho_p = 2.6;

    S.addDiskPacking(Rmin, Rmax, wall_gap, inlet_buf, outlet_buf, 1.001,
                     max_particles, wall_particles, rho_p, 0.0, 0.0);

    // All DOFs free for soil grains
    for (auto& B : S.bodies) {
        B->blockedDOFs  = Vector3r(1, 1, 1);
        B->blockedMDOFs = 1.0;
    }

    // ---------------- DEM time step (stability) ----------------
    const double Rref = 0.5 * (Rmin + Rmax);
    const double mass_ref = rho_p * M_PI * sqr(Rref);
    const double dt_crit = 0.1 * std::sqrt(mass_ref / S.normal_stiffness);
    const int dem_subcycles = std::max(1, static_cast<int>(S.dt_lbm / dt_crit) + 1);
    S.dt_dem = S.dt_lbm / dem_subcycles;

    // ---------------- Engines ----------------
    S.engines.clear();
    S.engines.push_back(std::make_shared<LatticeSearch>());
    S.engines.push_back(std::make_shared<FluidCollision>());
    S.engines.push_back(std::make_shared<ImbBoundary>());
    S.engines.push_back(std::make_shared<FluidStreaming>());

    std::vector<std::shared_ptr<Engine>> dem_engines;
    dem_engines.push_back(std::make_shared<ContactResolution>());
    dem_engines.push_back(std::make_shared<InteractionLoop>());
    dem_engines.push_back(std::make_shared<BodyLoop>());
    dem_engines.push_back(std::make_shared<Integrator>());
    dem_engines.push_back(std::make_shared<UpdateContact>());

    // ---------------- Drained triaxial: strain-controlled loading ----------------
    const double Ly0 = Ly;  // initial height
    const double strain_rate = 1e-5;  // axial strain per LBM step (slow = drained)

    // Reaction forces (recomputed from border penetration)
    auto compute_reactions = [&]() {
        double F_top = 0.0, F_bottom = 0.0, F_left = 0.0, F_right = 0.0;
        const double x_min = S.model_min_corner[0];
        const double x_max = S.model_max_corner[0];
        const double y_min = S.model_min_corner[1];
        const double y_max = S.model_max_corner[1];

        for (auto& B : S.bodies) {
            const double r = B->shape->radius;
            const double px = B->state->pos[0];
            const double py = B->state->pos[1];

            if ((py + r) > y_max) F_top   += S.border_stiffness * ((py + r) - y_max);
            if ((py - r) < y_min) F_bottom += S.border_stiffness * (y_min - (py - r));
            if ((px - r) < x_min) F_left   += S.border_stiffness * (x_min - (px - r));
            if ((px + r) > x_max) F_right  += S.border_stiffness * ((px + r) - x_max);
        }
        return std::make_tuple(F_top, F_bottom, F_left, F_right);
    };

    // Current sample dimensions (effective: from boundary positions)
    const double width = S.model_max_corner[0] - S.model_min_corner[0];
    double current_height = S.model_max_corner[1] - S.model_min_corner[1];

    // Stress: force per unit length (2D)
    auto stresses_from_reactions = [&](double F_top, double F_bottom, double F_left, double F_right) {
        current_height = S.model_max_corner[1] - S.model_min_corner[1];
        const double sigma_axial = (F_top + F_bottom) * 0.5 / width;   // vertical stress
        const double sigma_lat   = (F_left + F_right) * 0.5 / current_height;  // lateral stress
        // Convention: sigma_1 = major (axial in compression), sigma_3 = minor (confining)
        const double sigma_1 = std::max(sigma_axial, sigma_lat);
        const double sigma_3 = std::min(sigma_axial, sigma_lat);
        const double q = sigma_1 - sigma_3;  // deviator stress
        return std::make_tuple(sigma_axial, sigma_lat, q, sigma_1, sigma_3);
    };

    // ---------------- Output ----------------
    std::ofstream csv("drained_triaxial.csv");
    csv << "time,axial_strain,volumetric_strain,sigma_axial,sigma_lateral,q,sigma_1,sigma_3\n";

    const double max_axial_strain = 0.15;  // stop at 15% axial strain
    const int output_interval = 50;
    const int vtk_interval    = 500;

    std::cout << "Drained triaxial test (2D DEM, strain-controlled)\n";
    std::cout << "Sample: " << Lx << " x " << Ly << ", particles: " << S.bodies.size() << "\n";
    std::cout << "Strain rate: " << strain_rate << " /step, max axial strain: " << max_axial_strain << "\n\n";

    double axial_strain = 0.0;

    // ---------------- Time loop ----------------
    while (axial_strain < max_axial_strain) {

        // Impose axial strain: move top boundary down (drained = slow, no pore pressure)
        const double delta_h = strain_rate * Ly0 * S.dt_lbm;
        S.model_max_corner[1] -= delta_h;
        current_height = S.model_max_corner[1] - S.model_min_corner[1];
        axial_strain = (Ly0 - current_height) / Ly0;

        // DEM subcycles
        for (int sub = 0; sub < dem_subcycles; ++sub)
            for (auto& E : dem_engines)
                E->action();

        // Optional: update IMB geometry for particles (if you later enable coupling)
        for (auto& B : S.bodies) {
            auto& sh = *B->shape;
            sh.cloud.clear();
            const double r = sh.radius;
            for (double theta = 0.0; theta < 2.0 * M_PI; theta += 0.15) {
                double x = B->state->pos[0] + r * std::cos(theta);
                double y = B->state->pos[1] + r * std::sin(theta);
                boost::geometry::append(sh.cloud.outer(), point(x, y));
            }
            boost::geometry::correct(sh.cloud);
        }

        // LBM step (pass-through for dry case)
        for (auto& E : S.engines)
            E->action();

        // Output
        if (S.iter % output_interval == 0) {
            auto [F_top, F_bottom, F_left, F_right] = compute_reactions();
            auto [sigma_axial, sigma_lat, q, sigma_1, sigma_3] = stresses_from_reactions(F_top, F_bottom, F_left, F_right);

            const double V0 = Lx * Ly0;
            const double V  = width * current_height;
            const double volumetric_strain = (V0 - V) / V0;

            csv << S.time << ","
                << axial_strain << ","
                << volumetric_strain << ","
                << sigma_axial << ","
                << sigma_lat << ","
                << q << ","
                << sigma_1 << ","
                << sigma_3 << "\n";

            std::cout << "t=" << std::setw(8) << std::fixed << std::setprecision(0) << S.time
                      << "  e_ax=" << std::setw(6) << std::setprecision(4) << axial_strain
                      << "  e_vol=" << std::setw(6) << std::setprecision(4) << volumetric_strain
                      << "  q=" << std::scientific << std::setprecision(3) << q
                      << "  sigma_1=" << sigma_1
                      << "  sigma_3=" << sigma_3
                      << std::endl;
        }

        if (S.iter % vtk_interval == 0) {
            O.ExportFluidVtk("(Triax)Drained_Fluid_");
            O.ExportParticleVtk("(Triax)Drained_Particles_");
        }

        S.time += S.dt_lbm;
        ++S.iter;
    }

    csv.close();
    std::cout << "\nFinished. Results: drained_triaxial.csv\n";
    return 0;
}
