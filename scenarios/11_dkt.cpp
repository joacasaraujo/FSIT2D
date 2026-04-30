#include "scene/scene.h"    
#include "scene/output.h"
#include <iomanip>
#include <fstream>
#include <cmath>

// Drafting–Kissing–Tumbling de duas partículas em queda
// LBM–DEM–IMB 2D (discos)

int main() {
    Timer T;
    Scene& S = Scene::get_Scene();
    Output O;

    // ============================================================
    // 1) PARÂMETROS FÍSICOS
    // ============================================================

    double R  = 6.0;           // raio das partículas [lu]
    double D  = 2.0 * R;       // diâmetro

    // Domínio tipo coluna
    double Lx = 10.0 * D;      // ~10D de largura  (120)
    double Ly = 25.0 * D;      // ~25D de altura   (300)

    S.model_min_corner = {0.0, 0.0, 0.0};
    S.model_max_corner = {Lx, Ly, 1.0};

    // Fluido
    double rho_f = 1.0;
    // Viscosidade menor para aumentar Re e favorecer DKT (wake mais pronunciada)
    double nu    = 0.008;      // viscosidade cinemática (tau = 3*nu+0.5 ≈ 0.524)
    double mu    = rho_f * nu;

    // Partículas
    double rho_p     = 1.5;
    double delta_rho = rho_p - rho_f;

    // Queremos Re moderado (~45–60) para DKT mais evidente
    // U_t alvo ~ 0.08
    double U_target = 0.08;
    double g = 4.0 * mu * U_target / (delta_rho * R * R);

    double F_net_teo = delta_rho * M_PI * R * R * g; // peso efetivo 2D

    // ============================================================
    // 2) LBM
    // ============================================================

    S.dt_lbm           = 1.0;
    S.lattice_spacing  = 1.0;
    S.initial_density  = rho_f;
    S.initial_velocity = Vector3r::Zero();

    S.relaxation_time    = 3.0 * nu + 0.5;  // τ = 3ν + 0.5
    S.collision_operator = "BGK";
    S.dem_coupling       = "IMB";

    S.GUO_fluid_forcing = Vector3r::Zero();       // fluido quiescente
    S.gravity           = Vector3r(0.0, -g, 0.0); // gravidade só no DEM

    // ============================================================
    // 3) DEM (mais estável)
    // ============================================================

    S.friction_angle   = 30.0;
    S.local_damping    = 0.10;     // menos amortecimento para permitir DKT
    S.normal_stiffness = 8.0e2;
    S.shear_stiffness  = 4.0e2;
    S.border_stiffness = 3.0e3;

    double mass = rho_p * M_PI * R * R;

    // Subciclagem DEM fixa e conservadora
    // int dem_subcycles = 50;
    double dt_crit = 0.1 * std::sqrt(mass / S.normal_stiffness);  // 0.2 instead of 0.1
    int dem_subcycles = std::max(1, (int)(S.dt_lbm / dt_crit) + 1);
    S.dt_dem = S.dt_lbm / dem_subcycles;

    // ============================================================
    // 4) POSIÇÃO INICIAL – DKT (PARTÍCULAS PERTO!)
    // ============================================================

    double x_center   = 0.5 * Lx;
    // Posições ajustadas para evidenciar DKT:
    // - offset lateral maior para iniciar drafting/kissing mais cedo
    // - gap vertical menor para acelerar a aproximação (Kiss)
    double dx_offset  = 0.25 * R;   // offset lateral maior para quebrar simetria
    double y_upper    = 0.80 * Ly;  // partícula de cima (trailing)
    double gap        = 1.0 * D;    // gap menor para acelerar o kiss
    double y_lower    = y_upper - gap; // partícula de baixo (leading)

    // Disco 0: partícula de cima (trailing / upstream)
    Vector3r pos_upper(x_center + dx_offset, y_upper, 0.0);
    // Disco 1: partícula de baixo (leading / downstream, na esteira)
    Vector3r pos_lower(x_center - dx_offset, y_lower, 0.0);

    S.AddDisk(pos_upper, R, rho_p);   // body 0
    S.AddDisk(pos_lower, R, rho_p);   // body 1

    for (auto& B : S.bodies) {
        B->state->vel    = Vector3r::Zero();
        B->state->rotVel = 0.0;
    }

    // ============================================================
    // 5) MALHA FLUIDA / PAREDES
    // ============================================================

    S.AddRectangularCanal();

    for (auto& C : S.cells) {
        C->is_wall  = false;
        C->is_solid = false;
    }

    // paredes inferior e superior
    for (int i = 0; i < (int)Lx; ++i) {
        S.cells[S.CalculateCellId(i, 0)]->is_wall = true;
        S.cells[S.CalculateCellId(i, (int)Ly - 1)]->is_wall = true;
    }
    // paredes laterais
    for (int j = 0; j < (int)Ly; ++j) {
        S.cells[S.CalculateCellId(0, j)]->is_wall           = true;
        S.cells[S.CalculateCellId((int)Lx - 1, j)]->is_wall = true;
    }

    for (auto& C : S.cells) {
        C->set_initial_condition();
    }

    // ============================================================
    // 6) ENGINES
    // ============================================================

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

    // ============================================================
    // 7) SAÍDA
    // ============================================================

    std::ofstream csv("dkt_two_disks_results.csv");
    csv << "time,"
           "y_upper,vy_upper,Fy_upper,"
           "y_lower,vy_lower,Fy_lower,"
           "dy\n";

    std::cout << "============================================================\n";
    std::cout << "   Drafting–Kissing–Tumbling - duas partículas (LBM–DEM–IMB)\n";
    std::cout << "============================================================\n\n";
    std::cout << "R = " << R << ",  Lx = " << Lx << ", Ly = " << Ly << "\n";
    std::cout << "nu = " << nu << ", tau = " << S.relaxation_time << "\n";
    std::cout << "rho_p = " << rho_p << ", rho_f = " << rho_f << "\n";
    std::cout << "g = "   << g << ", F_net_teo = " << F_net_teo << "\n\n";
    std::cout << "DEM subcycles: " << dem_subcycles << "\n\n";

    // ============================================================
    // 8) LOOP DE TEMPO
    // ============================================================

    double max_time        = 100000.0;
    int    output_interval = 200;
    int    vtk_interval    = 2000;

    while (S.time < max_time) {

        // fusível simples para divergência
        for (auto& Bptr : S.bodies) {
            if (!std::isfinite(Bptr->state->pos[0]) ||
                !std::isfinite(Bptr->state->pos[1]) ||
                std::fabs(Bptr->state->pos[1]) > 1e6) {
                std::cout << "Divergência detectada em t = " << S.time << std::endl;
                goto end_simulation;
            }
        }

        // Atualizar geometria IMB de TODAS as partículas
        for (auto& Bptr : S.bodies) {
            auto& B = *Bptr;
            B.shape->cloud.clear();
            for (double theta = 0.0; theta < 2.0 * M_PI; theta += 0.1) {
                double x = B.state->pos[0] + R * std::cos(theta);
                double y = B.state->pos[1] + R * std::sin(theta);
                boost::geometry::append(B.shape->cloud.outer(), point(x, y));
            }
            boost::geometry::correct(B.shape->cloud);
        }

        // LBM
        for (auto& E : S.engines) {
            E->action();
        }

        // DEM (subciclagem)
        for (int sub = 0; sub < dem_subcycles; ++sub) {
            for (auto& E : dem_engines) {
                E->action();
            }
        }

        // Saída
        if (S.iter % output_interval == 0) {
            auto& B_up = *S.bodies[0];
            auto& B_lo = *S.bodies[1];

            double vy_up = -B_up.state->vel[1];  // positivo para baixo
            double vy_lo = -B_lo.state->vel[1];
            double dy    = B_up.state->pos[1] - B_lo.state->pos[1]; // distância vertical

            csv << S.time << ","
                << B_up.state->pos[1] << "," << vy_up << "," << B_up.state->hydro_force[1] << ","
                << B_lo.state->pos[1] << "," << vy_lo << "," << B_lo.state->hydro_force[1] << ","
                << dy << "\n";

            std::cout << "t=" << std::setw(7) << std::fixed << std::setprecision(0) << S.time
                      << "  y_up=" << std::setw(7) << std::setprecision(1) << B_up.state->pos[1]
                      << "  y_lo=" << std::setw(7) << B_lo.state->pos[1]
                      << "  vy_up=" << std::scientific << std::setprecision(3) << vy_up
                      << "  vy_lo=" << vy_lo
                      << "  dy="    << std::fixed << std::setprecision(2) << dy
                      << std::endl;
        }

        // VTK
        if (S.iter % vtk_interval == 0) {
            O.ExportFluidVtk("DKT_Fluid_");
            O.ExportParticleVtk("DKT_Particles_");
        }

        // parar se ambas perto do fundo
        bool all_low = true;
        for (auto& Bptr : S.bodies) {
            if (Bptr->state->pos[1] > 3.0 * R) {
                all_low = false;
                break;
            }
        }
        if (all_low) {
            std::cout << "\nPartículas próximas ao fundo. Encerrando.\n\n";
            break;
        }

        S.time += S.dt_lbm;
        ++S.iter;
    }

end_simulation:

    csv.close();
    std::cout << "Simulação DKT finalizada. Resultados em dkt_two_disks_results.csv\n";
    std::cout << "VTKs: DKT_Fluid_*.vtk, DKT_Particles_*.vtk\n";

    return 0;
}
