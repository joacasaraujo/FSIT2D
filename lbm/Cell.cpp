#include "cell.h"
#include "scene/scene.h"

double Cell::CalculateEqFunction(double rho, Vector3r u, int k) {
    const Scene&  S  = Scene::get_Scene();
    const double  c  = S.latticeSpeed;                 // usually 1.0
    const double  cs2 = (c*c)/3.0;                     // D2Q9
    const double  inv_cs2 = 1.0 / cs2;
    const double  inv_cs4 = inv_cs2 * inv_cs2;

    const double cu = lattice->discrete_velocity[k].dot(u);
    const double u2 = u.squaredNorm();

    return lattice->node_weight[k] * rho *
           (1.0 + inv_cs2*cu + 0.5*inv_cs4*cu*cu - 0.5*inv_cs2*u2);
}

void Cell::CalculateDensityAndVelocity() {
    const Scene&  S  = Scene::get_Scene();
    if (is_solid || is_wall || solid_fraction >= 0.99) {
        material->density = 0.0;       
        state->vel.setZero();
        return;
    }

    double   rho = 0.0;
    Vector3r mom = Vector3r::Zero();

    for (int k = 0; k < lattice->number_of_nodes; ++k) {
        const double fk = lattice->f[k];

        if (!std::isfinite(fk)) {
            std::cout << "NaN/Inf in f at cell " << id << " pos = " << state->pos.transpose() << " k = " << k << " fk = " << fk << "\n";
            std::abort();
        }
        
        rho += fk;
        mom += lattice->discrete_velocity[k] * fk;
    }

    material->density = rho;

    if (rho > 1e-14) {
        Vector3r mom_eff = mom + 0.5 * S.GUO_fluid_forcing * S.dt_lbm;
        state->vel = mom_eff / rho;
        state->vel[2] = 0.0; 
    } else {
        state->vel.setZero();
    }

    ASSERT_MSG(std::isfinite(material->density), "NaN density");
    ASSERT_MSG(std::isfinite(state->vel[0]) && std::isfinite(state->vel[1]), "NaN velocity");
}


void Cell::set_neighbor_node(){
    const Scene& S = Scene::get_Scene();

    lattice->neighbor_node.clear();
    lattice->neighbor_node.reserve(lattice->number_of_nodes);

    const int Nx = static_cast<int>(S.model_max_corner[0]);
    const int Ny = static_cast<int>(S.model_max_corner[1]);

    for (int k = 0; k < lattice->number_of_nodes; ++k){
        Vector3r aux = state->pos + lattice->discrete_velocity[k];

        int ax = static_cast<int>(aux[0]);
        int ay = static_cast<int>(aux[1]);

        if (ax < 0) ax = Nx - 1;
        if (ay < 0) ay = Ny - 1;
        if (ax >= Nx) ax = 0;
        if (ay >= Ny) ay = 0;

        lattice->neighbor_node.push_back(ax + ay * Nx);
    }
}

void Cell::set_initial_condition(){
    const Scene& S = Scene::get_Scene();
    for (int k = 0; k < lattice->number_of_nodes; ++k)
        lattice->f[k] = CalculateEqFunction(S.initial_density, S.initial_velocity, k);
    CalculateDensityAndVelocity();
}

void Cell::set_mrt_parameters(){
    Scene& S = Scene::get_Scene();

    lattice->m << 1,1,1,1,1,1,1,1,1,
                 -4,-1,-1,-1,-1,2,2,2,2,
                  4,-2,-2,-2,-2,1,1,1,1,
                  0,1,0,-1,0,1,-1,-1,1,
                  0,-2,0,2,0,1,-1,-1,1,
                  0,0,1,0,-1,1,1,-1,-1,
                  0,0,-2,0,2,1,1,-1,-1,
                  0,1,-1,1,-1,0,0,0,0,
                  0,0,0,0,0,1,-1,1,-1;
    lattice->m_inv = lattice->m.inverse();
    lattice->relaxation_matrix = S.relaxation_parameters.asDiagonal();
}