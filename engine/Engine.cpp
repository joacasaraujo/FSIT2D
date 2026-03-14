#include "Engine.h"
#include "scene/Scene.h"

void ZouAndHeBC::action() {
    Scene& S = Scene::get_Scene();

    const double div = 1.0 / 6.0;
    const double aux = 2.0 / 3.0;

    const int Nx = static_cast<int>(S.model_max_corner[0]);
    const int Ny = static_cast<int>(S.model_max_corner[1]);

    // ---------------- Velocity BCs ----------------
    switch (S.apply_velocity_bc) {
        // case 1: top (j = Ny-1)
        case 1: {
            for (int i = 0; i < S.model_max_corner[0]; ++i) {
                int id = S.CalculateCellId(i, S.model_max_corner[1] - 1);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                const double ux = C->lattice->velocity_bc[0];
                const double uy = 0.0;

                // Reconstruct rho using known populations (Zou & He)
                double rho = (C->lattice->f[0] + C->lattice->f[1] + C->lattice->f[3]
                            + 2.0 * (C->lattice->f[2] + C->lattice->f[5] + C->lattice->f[6]))
                            / (1.0 - uy);

                // Unknown distributions at north boundary: f4,f7,f8
                C->lattice->f[4] = C->lattice->f[2];
                C->lattice->f[8] = C->lattice->f[6]
                                 + 0.5 * (C->lattice->f[3] - C->lattice->f[1])
                                 + (1.0 / 6.0) * rho * ux;
                C->lattice->f[7] = C->lattice->f[5]
                                 + 0.5 * (C->lattice->f[1] - C->lattice->f[3])
                                 - (1.0 / 6.0) * rho * ux;

                C->material->density = rho;
                C->state->vel = Vector3r(ux, uy, 0.0);
            }
        } break;

        // case 2: right (i = Nx-1)
        case 2: {
            int i = Nx - 1;
            for (int j = 1; j <= Ny - 2; ++j) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double rho = (C->lattice->f[0] + C->lattice->f[2] + C->lattice->f[4]
                             + 2.0*(C->lattice->f[1] + C->lattice->f[5] + C->lattice->f[8]))
                             / (1.0 + C->lattice->velocity_bc[0]);

                C->lattice->f[3] = C->lattice->f[1] - aux * rho * C->lattice->velocity_bc[0];
                C->lattice->f[7] = C->lattice->f[5] - div * rho * C->lattice->velocity_bc[0]
                                   + 0.5 * (C->lattice->f[2] - C->lattice->f[4]);
                C->lattice->f[6] = C->lattice->f[8] - div * rho * C->lattice->velocity_bc[0]
                                   + 0.5 * (C->lattice->f[4] - C->lattice->f[2]);
                C->CalculateDensityAndVelocity();
            }
        } break;

        // case 3: bottom (j = 0)
        case 3: {
            int j = 0;
            for (int i = 1; i <= Nx - 2; ++i) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double rho = (C->lattice->f[0] + C->lattice->f[1] + C->lattice->f[3]
                             + 2.0*(C->lattice->f[4] + C->lattice->f[7] + C->lattice->f[8]))
                             / (1.0 - C->lattice->velocity_bc[1]);

                C->lattice->f[2] = C->lattice->f[4] + aux * rho * C->lattice->velocity_bc[1];
                C->lattice->f[5] = C->lattice->f[7] + div * rho * C->lattice->velocity_bc[1]
                                   - 0.5 * (C->lattice->f[1] - C->lattice->f[3]);
                C->lattice->f[6] = C->lattice->f[8] + div * rho * C->lattice->velocity_bc[1]
                                   - 0.5 * (C->lattice->f[3] - C->lattice->f[1]);
                C->CalculateDensityAndVelocity();
            }
        } break;

        // case 4: left (i = 0)  <-- Poiseuille inlet (ux)
        case 4: {
            int i = 0;
            for (int j = 1; j <= Ny - 2; ++j) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double rho = (C->lattice->f[0] + C->lattice->f[2] + C->lattice->f[4]
                             + 2.0*(C->lattice->f[3] + C->lattice->f[6] + C->lattice->f[7]))
                             / (1.0 - C->lattice->velocity_bc[0]);

                C->lattice->f[1] = C->lattice->f[3] + aux * rho * C->lattice->velocity_bc[0];
                C->lattice->f[5] = C->lattice->f[7] + div * rho * C->lattice->velocity_bc[0]
                                   - 0.5 * (C->lattice->f[2] - C->lattice->f[4]);
                C->lattice->f[8] = C->lattice->f[6] + div * rho * C->lattice->velocity_bc[0]
                                   - 0.5 * (C->lattice->f[4] - C->lattice->f[2]);
                C->CalculateDensityAndVelocity();
            }
        } break;
    }

    // ---------------- Density BCs ----------------
    switch (S.apply_density_bc) {
        // case 1: top density
        case 1: {
            int j = Ny - 1;
            for (int i = 1; i <= Nx - 2; ++i) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double vy = -1.0 + (C->lattice->f[0] + C->lattice->f[1] + C->lattice->f[3]
                                   + 2.0*(C->lattice->f[2] + C->lattice->f[5] + C->lattice->f[6]))
                                   / C->lattice->density_bc;

                C->lattice->f[4] = C->lattice->f[2] - aux * C->lattice->density_bc * vy;
                C->lattice->f[7] = C->lattice->f[5] - div * C->lattice->density_bc * vy
                                   + 0.5 * (C->lattice->f[1] - C->lattice->f[3]);
                C->lattice->f[8] = C->lattice->f[6] - div * C->lattice->density_bc * vy
                                   + 0.5 * (C->lattice->f[3] - C->lattice->f[1]);
                C->CalculateDensityAndVelocity();
            }
        } break;

        // case 2: right density  <-- Poiseuille outlet (rho)
        case 2: {
            int i = Nx - 1;
            for (int j = 1; j <= Ny - 2; ++j) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double vx = -1.0 + (C->lattice->f[0] + C->lattice->f[2] + C->lattice->f[4]
                                   + 2.0*(C->lattice->f[1] + C->lattice->f[5] + C->lattice->f[8]))
                                   / C->lattice->density_bc;

                C->lattice->f[3] = C->lattice->f[1] - aux * C->lattice->density_bc * vx;
                C->lattice->f[7] = C->lattice->f[5] - div * C->lattice->density_bc * vx
                                   + 0.5 * (C->lattice->f[2] - C->lattice->f[4]);
                C->lattice->f[6] = C->lattice->f[8] - div * C->lattice->density_bc * vx
                                   + 0.5 * (C->lattice->f[4] - C->lattice->f[2]);
                C->CalculateDensityAndVelocity();
            }
        } break;

        // case 3: bottom density
        case 3: {
            int j = 0;
            for (int i = 1; i <= Nx - 2; ++i) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double vy = -1.0 + (C->lattice->f[0] + C->lattice->f[1] + C->lattice->f[3]
                                   + 2.0*(C->lattice->f[4] + C->lattice->f[7] + C->lattice->f[8]))
                                   / C->lattice->density_bc;

                C->lattice->f[2] = C->lattice->f[4] - aux * C->lattice->density_bc * vy;
                C->lattice->f[5] = C->lattice->f[7] - div * C->lattice->density_bc * vy
                                   + 0.5 * (C->lattice->f[3] - C->lattice->f[1]);
                C->lattice->f[6] = C->lattice->f[8] - div * C->lattice->density_bc * vy
                                   + 0.5 * (C->lattice->f[1] - C->lattice->f[3]);
                C->CalculateDensityAndVelocity();
            }
        } break;

        // case 4: left density
        case 4: {
            int i = 0;
            for (int j = 1; j <= Ny - 2; ++j) {
                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];
                if (C->is_solid || C->is_wall) continue;

                double vx = -1.0 + (C->lattice->f[0] + C->lattice->f[2] + C->lattice->f[4]
                                   + 2.0*(C->lattice->f[3] + C->lattice->f[6] + C->lattice->f[7]))
                                   / C->lattice->density_bc;

                C->lattice->f[1] = C->lattice->f[3] + aux * C->lattice->density_bc * vx;
                C->lattice->f[5] = C->lattice->f[7] + div * C->lattice->density_bc * vx
                                   + 0.5 * (C->lattice->f[4] - C->lattice->f[2]);
                C->lattice->f[8] = C->lattice->f[6] + div * C->lattice->density_bc * vx
                                   + 0.5 * (C->lattice->f[2] - C->lattice->f[4]);
                C->CalculateDensityAndVelocity();
            }
        } break;
    }
}

void FluidCollision::action() {
    Scene& S = Scene::get_Scene();

    for (auto& C : S.cells) {
        if (C->is_solid || C->is_wall || C->solid_fraction > 0.0) {
            // No collision for solid/wall cells
            continue;
        }

        // Compute equilibrium distribution fEq
        const double rho = C->material->density;
        const Vector3r u = C->state->vel;
        Vector9r fEq = Vector9r::Zero();
        for (int k = 0; k < C->lattice->number_of_nodes; ++k) {
            fEq[k] = C->CalculateEqFunction(rho, u, k);
        }

        // MRT relaxation
        if (S.collision_operator == "MRT"){
            C->lattice->f -= C->lattice->m_inv *
                            C->lattice->relaxation_matrix *
                            C->lattice->m *
                            (C->lattice->f - fEq);
        }
        else {
            // BGK relaxation
            C->lattice->f -= (1.0 / S.relaxation_time) * (C->lattice->f - fEq);
        }
    }
}

void FluidStreaming::action() {
    Scene& S = Scene::get_Scene();

    // 1) Clear f_tmp
    for (auto& Cptr : S.cells) {
        Cell& C = *Cptr;
        C.lattice->f_tmp.setZero();
    }

    // 2) Push-stream with link-wise bounce-back
    for (size_t idx = 0; idx < S.cells.size(); ++idx) {
        Cell& C = *S.cells[idx];
        if (C.is_solid || C.is_wall) continue;

        for (int k = 0; k < C.lattice->number_of_nodes; ++k) {
            const int nb = C.lattice->neighbor_node[k];
            const int ko = C.lattice->opposite_node[k];

            Cell& N = *S.cells[nb];

            if (N.is_solid || N.is_wall) {
                // reflect back into current cell, opposite direction
                C.lattice->f_tmp[ko] += C.lattice->f[k];
            } else {
                // normal push to neighbor
                N.lattice->f_tmp[k] += C.lattice->f[k];
            }
        }
    }

    // 3) Swap and update macros (fluids only)
    for (auto& Cptr : S.cells) {
        Cell& C = *Cptr;

        std::swap(C.lattice->f, C.lattice->f_tmp);

        if (C.is_solid || C.is_wall) {
            C.material->density = 0.0;             
            C.state->vel = Vector3r::Zero();
            continue;
        }

        C.CalculateDensityAndVelocity();
    }
}

void ApplyFluidForcing::action(){
    Scene& S = Scene::get_Scene();
    const double cs2 = 1.0 / 3.0;
    for (auto& Cptr : S.cells) {
            auto& C = *Cptr;
            if (C.is_solid || C.is_wall) continue;
            Vector3r F = S.GUO_fluid_forcing;
            for (int k = 0; k < C.lattice->number_of_nodes; ++k) {
                const auto& ci = C.lattice->discrete_velocity[k];
                const double ciF = ci.dot(F);
                const double ciu = ci.dot(C.state->vel);
                const double uF  = C.state->vel.dot(F);
                const double term = (ciF/cs2) + (ciu*ciF)/(cs2*cs2) - (uF/cs2);
                C.lattice->f[k] += (1.0 - 0.5 / S.relaxation_time)
                                 * C.lattice->node_weight[k] * term;
            }
        }
}

void LatticeSearch::action() {
    Scene& S = Scene::get_Scene();
    if (S.bodies.empty()) return;

    // ---------- 0) Cache previous cell state (no class changes needed) ----------
    const std::size_t nCells = S.cells.size();
    std::vector<double> prev_chi(nCells, 0.0);
    std::vector<bool>   prev_is_solid(nCells, false);

    for (std::size_t id = 0; id < nCells; ++id) {
        auto& C = S.cells[id];
        prev_chi[id]      = C->solid_fraction;
        prev_is_solid[id] = C->is_solid;
    }

    // ---------- 1) Reset fields for recalculation (do NOT touch f here) ----------
    for (auto& C : S.cells) {
        C->solid_fraction = 0.0;

        // Keep static walls unchanged; reset only non-wall solids
        if (!C->is_wall) {
            C->is_solid = false;
        }
    }

    // ---------- 2) Recompute solid fractions by geometry intersection ----------
    aabb box;

    for (auto& B : S.bodies) {
        B->fluid_interactions.clear();

        // Bounding box of body polygon
        boost::geometry::envelope(B->shape->cloud, box);

        int aabb_x_min = static_cast<int>(box.min_corner().get<0>());
        int aabb_x_max = static_cast<int>(box.max_corner().get<0>() + 1);
        int aabb_y_min = static_cast<int>(box.min_corner().get<1>());
        int aabb_y_max = static_cast<int>(box.max_corner().get<1>() + 1);

        // Clamp to domain bounds
        aabb_x_min = std::max(aabb_x_min, static_cast<int>(S.model_min_corner[0]));
        aabb_y_min = std::max(aabb_y_min, static_cast<int>(S.model_min_corner[1]));
        aabb_x_max = std::min(aabb_x_max, static_cast<int>(S.model_max_corner[0]));
        aabb_y_max = std::min(aabb_y_max, static_cast<int>(S.model_max_corner[1]));

        for (int i = aabb_x_min; i < aabb_x_max; ++i) {
            for (int j = aabb_y_min; j < aabb_y_max; ++j) {

                int id = S.CalculateCellId(i, j);
                auto& C = S.cells[id];

                // Skip static walls entirely (no IMB work, no interactions list)
                if (C->is_wall) {
                    continue;
                }

                // Fast reject: if grid doesn't intersect particle polygon
                if (!boost::geometry::intersects(C->grid, B->shape->cloud)) {
                    continue;
                }

                // This cell participates in IMB for this body
                B->fluid_interactions.push_back(C);

                // Compute intersection area between cell and particle polygon
                std::vector<polygon> output;
                boost::geometry::intersection(C->grid, B->shape->cloud, output);

                double overlap_area = 0.0;
                for (const auto& poly : output) {
                    overlap_area += boost::geometry::area(poly);
                }

                const double cell_area = boost::geometry::area(C->grid);
                double phi = 0.0;

                if (overlap_area > 0.0) {
                    phi = overlap_area / cell_area;
                } else {
                    // Fallback: if cell center is inside particle -> full
                    point cell_center(C->state->pos[0], C->state->pos[1]);
                    if (boost::geometry::within(cell_center, B->shape->cloud)) {
                        phi = 1.0;
                    }
                }

                // Accumulate for multiple bodies (clamp to [0,1])
                phi = std::min(1.0, std::max(0.0, phi));
                C->solid_fraction = std::min(1.0, C->solid_fraction + phi);

                if (C->solid_fraction >= 0.999) {
                    C->solid_fraction = 1.0;
                    C->is_solid = true;
                }
            }
        }
    }

    // ---------- 3) Now apply reinitialization ONLY for cells that truly became fluid ----------
    // Condition: previously had solid content (or was_solid) and now is pure fluid (chi==0)
    for (std::size_t id = 0; id < nCells; ++id) {
        auto& C = S.cells[id];

        if (C->is_wall) continue;

        const double chi_old = prev_chi[id];
        const bool   solid_old = prev_is_solid[id];
        const double chi_new = C->solid_fraction;

        // "Became fluid": old influenced by solid, new is fully fluid
        if ((chi_old > 0.0 || solid_old) && chi_new == 0.0) {

            // Option A (hard reset): set to equilibrium at initial macros
            // for (int k = 0; k < C->lattice->number_of_nodes; ++k) {
            //     C->lattice->f[k] = C->CalculateEqFunction(S.initial_density, S.initial_velocity, k);
            // }
            // C->material->density = S.initial_density;
            // C->state->vel = S.initial_velocity;

            // Option B (softer reset): uncomment to blend instead of hard reset
            const double alpha = 0.3; // 0<alpha<=1
            Vector9r fEq = Vector9r::Zero();
            for (int k = 0; k < C->lattice->number_of_nodes; ++k) {
                fEq[k] = C->CalculateEqFunction(S.initial_density, S.initial_velocity, k);
            }
            C->lattice->f = (1.0 - alpha) * C->lattice->f + alpha * fEq;
            C->material->density = S.initial_density;
            C->state->vel = S.initial_velocity;
        }
    }
}

// void ImbBoundary::action(){
//     Scene& S = Scene::get_Scene();
    
//     if (S.bodies.empty()) return;

//     // Calculate solid fraction in each cell for every body
//     for (auto& B : S.bodies) {
//         B->state->hydro_force = Vector3r::Zero();
//         B->state->hydro_torque = 0.0;
        
//         for (auto& C : B->fluid_interactions) {
//             double tau = S.relaxation_time;
//             double chi = C->solid_fraction;

//             // --- basic safety on chi ---
//             if (chi <= 0.0) continue;  // Pure fluid cell, skip IMB
//             if (chi > 1.0) chi = 1.0;
//             double Bn = chi * (tau - 0.5) / ((1.0 - chi) + (tau - 0.5));  // weighting factor

//             // Skip solid/wall cells
//             if (C->is_solid || C->is_wall) continue;

//             // --- For fully solid cells (chi >= 0.99), set f to equilibrium with body velocity ---
//             // This ensures valid f values when the cell becomes fluid again
//             if (chi >= 0.99) {
//                 double rho_eq = S.initial_density;  // Use reference density
//                 Vector3r ub = B->state->vel;
//                 for (int i = 0; i < C->lattice->number_of_nodes; ++i) {
//                     C->lattice->f[i] = C->CalculateEqFunction(rho_eq, ub, i);
//                 }
//                 C->material->density = rho_eq;
//                 C->state->vel = ub;
//                 continue;
//             }

//             // --- macroscopic variables at the fluid cell ---
//             double   rho = C->material->density;
//             Vector3r uf  = C->state->vel;            // fluid velocity
//             Vector3r ub  = B->state->vel;            // body velocity

//             // Safety check: if macros are invalid, reset to equilibrium
//             if (!std::isfinite(rho) || rho < 0.1 || !std::isfinite(uf.norm())) {
//                 rho = S.initial_density;
//                 uf = ub;  // Use body velocity as reference
//                 for (int i = 0; i < C->lattice->number_of_nodes; ++i) {
//                     C->lattice->f[i] = C->CalculateEqFunction(rho, uf, i);
//                 }
//                 C->material->density = rho;
//                 C->state->vel = uf;
//             }

//             // --- compute equilibria and omega_s (Noble & Torczynski) ---
//             Vector9r feq_f      = Vector9r::Zero();  // f^eq(ρ, u_fluid, i)
//             Vector9r feq_b      = Vector9r::Zero();  // f^eq(ρ, u_body, i)
//             Vector9r feq_b_opp  = Vector9r::Zero();  // f^eq(ρ, u_body, opp)
//             Vector9r omega_s    = Vector9r::Zero();  // solid collision operator

//             for (int i = 0; i < C->lattice->number_of_nodes; ++i) {
//                 int opp = C->lattice->opposite_node[i];

//                 feq_f[i]     = C->CalculateEqFunction(rho, uf, i);      // eq at direction i with fluid vel
//                 feq_b[i]     = C->CalculateEqFunction(rho, ub, i);      // eq at direction i with body vel
//                 feq_b_opp[i] = C->CalculateEqFunction(rho, ub, opp);    // eq at direction opp with body vel

//                 // Noble & Torczynski:
//                 // Ω_i^s = [f_opp - f^eq(ρ,u_b,opp)] - [f_i - f^eq(ρ,u_f,i)]
//                 omega_s[i] = (C->lattice->f[opp] - feq_b_opp[i])
//                            - (C->lattice->f[i]   - feq_f[i]);
//             }

//            // --- collision + immersed boundary in one update ---

//             // Pré-calcula o termo de colisão MRT (vetorial) uma única vez
//             Vector9r df_coll = Vector9r::Zero();
//             if (S.collision_operator != "BGK") {
//                 // df_coll = M^{-1} S M (f - feq_f)
//                 df_coll = C->lattice->m_inv
//                         * C->lattice->relaxation_matrix
//                         * C->lattice->m
//                         * (C->lattice->f - feq_f);
//             }

//             for (int i = 0; i < C->lattice->number_of_nodes; ++i) {
//                 const double fi = C->lattice->f[i];
//                 double fi_new = fi;

//                 if (S.collision_operator == "BGK") {
//                     const double coll_i = (1.0 / tau) * (fi - feq_f[i]);   // (fi - feq)/tau
//                     fi_new = fi - (1.0 - Bn) * coll_i + Bn * omega_s[i];
//                 } else {
//                     // MRT: usa o componente i do vetor df_coll
//                     fi_new = fi - (1.0 - Bn) * df_coll[i] + Bn * omega_s[i];
//                 }

//                 // Safety
//                 if (!std::isfinite(fi_new) || fi_new < -10.0 || fi_new > 10.0) {
//                     fi_new = feq_f[i];
//                 }

//                 C->lattice->f[i] = fi_new;
//             }

//             // Noble & Torczynski force calculation (Eq. 31):
//             // F_p = Σ_n B_n × (Σ_i Ω_i^s × c_i)
//             // Apply only to interface cells (avoid pure solid)
//             if (chi < 0.99) {
//                 Vector3r cell_force = Vector3r::Zero();
//                 for (int i = 0; i < C->lattice->number_of_nodes; ++i) {
//                     cell_force += omega_s[i] * C->lattice->discrete_velocity[i];
//                 }
//                 // Reaction force on the solid
//                 B->state->hydro_force -= Bn * cell_force;
                
//                 // Noble & Torczynski torque calculation (Eq. 32):
//                 Vector3r r = C->state->pos - B->state->pos;  // x_n - X_p
//                 // 2D cross product: r × F = r_x * F_y - r_y * F_x
//                 B->state->hydro_torque -= Bn * (r[0] * cell_force[1] - r[1] * cell_force[0]);
//             }
//         }
//     }
// }


void ImbBoundary::action()
{
    Scene& S = Scene::get_Scene();
    if (S.bodies.empty()) return;

    const double tau = S.relaxation_time;
    const int Q = 9; // D2Q9

    // Tunables (good defaults)
    const double chi_eps    = 1e-12;  // treat as pure fluid
    const double chi_full_0 = 0.98;   // start blending toward solid
    const double chi_full_1 = 0.999;  // almost solid
    const double f_clamp_lo = -10.0;  // safety clamp
    const double f_clamp_hi =  10.0;

    // Optional temporal smoothing of chi (requires a "prev_chi" stored per cell)
    // If you don't have it, set use_chi_smoothing=false and it will use raw chi.
    const bool   use_chi_smoothing = false;
    const double beta_chi = 0.3; // 0..1 (smaller = more smoothing)

    auto smoothstep = [](double x) {
        x = std::clamp(x, 0.0, 1.0);
        return x*x*(3.0 - 2.0*x);
    };

    // Reset hydrodynamic accumulators
    for (auto& B : S.bodies) {
        B->state->hydro_force  = Vector3r::Zero();
        B->state->hydro_torque = 0.0;
    }

    // Loop bodies and their interacting cells
    for (auto& B : S.bodies)
    {
        const Vector3r ub = B->state->vel;

        for (auto& C : B->fluid_interactions)
        {
            if (C->is_wall) continue; // never IMB walls

            // --- chi (solid fraction) ---
            double chi = C->solid_fraction;
            chi = std::clamp(chi, 0.0, 1.0);

            // If you have stored previous chi per cell, you can smooth it here.
            // Example: C->prev_solid_fraction (not shown in your classes)
            if (use_chi_smoothing) {
                // chi = (1-beta_chi)*C->prev_solid_fraction + beta_chi*chi;
                // C->prev_solid_fraction = chi;
            }

            // Pure fluid cell: do nothing here (let plain collision engine handle it)
            // Since you're removing FluidCollision, we must still collide pure fluid cells
            // elsewhere (either keep FluidCollision OR include all fluid cells in this engine).
            // In your current architecture, B->fluid_interactions shouldn't include pure fluid cells,
            // so this guard is correct.
            if (chi <= chi_eps) continue;

            // If lattice search flags as solid, skip (no collision update here)
            // IMPORTANT: we still want to "sanitize" near-solid cells BEFORE skipping,
            // because they may later become fluid again.
            // We'll handle near-solid first, then skip solids.
            // -------------------------------------------------------------------------

            // --- safety for macros ---
            double   rho = C->material->density;
            Vector3r uf  = C->state->vel;

            if (!std::isfinite(rho) || rho < 0.1 || !std::isfinite(uf[0]) || !std::isfinite(uf[1])) {
                rho = S.initial_density;
                uf  = ub;
                for (int i = 0; i < Q; ++i) {
                    C->lattice->f[i] = C->CalculateEqFunction(rho, uf, i);
                }
                C->material->density = rho;
                C->state->vel        = uf;
            }

            // --- compute Bn (Noble & Torczynski weighting) ---
            // Bn = chi*(tau-0.5) / ((1-chi) + (tau-0.5))
            const double denom = (1.0 - chi) + (tau - 0.5);
            const double Bn = (denom > 0.0) ? (chi * (tau - 0.5) / denom) : 0.0;

            // --- near-solid stabilization: blend f -> feq(rho, ub) smoothly as chi -> 1 ---
            // This avoids "memory" in cells that toggle between solid/fluid and reduces spikes.
            if (chi >= chi_full_0) {
                double x = (chi - chi_full_0) / (chi_full_1 - chi_full_0);
                const double gamma = smoothstep(x); // 0..1
                for (int i = 0; i < Q; ++i) {
                    const double feq_bi = C->CalculateEqFunction(rho, ub, i);
                    C->lattice->f[i] = (1.0 - gamma) * C->lattice->f[i] + gamma * feq_bi;
                }
                C->state->vel = ub; // optional: keep macro velocity consistent near solid
            }

            // Now if the cell is marked solid, skip further work
            if (C->is_solid) continue;

            // --- compute feq_f and omega_s (Noble & Torczynski) ---
            Vector9r feq_f   = Vector9r::Zero();
            Vector9r omega_s = Vector9r::Zero();

            for (int i = 0; i < Q; ++i)
            {
                const int opp = C->lattice->opposite_node[i];

                const double feq_f_i     = C->CalculateEqFunction(rho, uf, i);
                const double feq_b_opp   = C->CalculateEqFunction(rho, ub, opp);

                feq_f[i] = feq_f_i;

                // Ω_i^s = [f_opp - f^eq(ρ,u_b,opp)] - [f_i - f^eq(ρ,u_f,i)]
                omega_s[i] = (C->lattice->f[opp] - feq_b_opp)
                           - (C->lattice->f[i]   - feq_f_i);
            }

            // --- collision term (precompute MRT vector once) ---
            Vector9r df_coll = Vector9r::Zero();
            const bool isBGK = (S.collision_operator == "BGK");

            if (!isBGK) {
                // df_coll = M^{-1} S M (f - feq_f)
                df_coll = C->lattice->m_inv
                        * C->lattice->relaxation_matrix
                        * C->lattice->m
                        * (C->lattice->f - feq_f);
            }

            // --- update distributions: f_new = f - (1-Bn)*Collision + Bn*omega_s ---
            for (int i = 0; i < Q; ++i)
            {
                const double fi = C->lattice->f[i];
                double fi_new   = fi;

                if (isBGK) {
                    const double coll_i = (1.0 / tau) * (fi - feq_f[i]);     // BGK
                    fi_new = fi - (1.0 - Bn) * coll_i + Bn * omega_s[i];
                } else {
                    // MRT: df_coll is already the collision operator in velocity space
                    fi_new = fi - (1.0 - Bn) * df_coll[i] + Bn * omega_s[i];
                }

                // Safety clamp (keeps run from exploding; still prefer fixing physics)
                if (!std::isfinite(fi_new) || fi_new < f_clamp_lo || fi_new > f_clamp_hi) {
                    fi_new = feq_f[i];
                }

                C->lattice->f[i] = fi_new;
            }

            // --- hydrodynamic force/torque (Noble & Torczynski Eq. 31/32 style) ---
            // F_cell = Bn * Σ_i omega_s[i] * c_i   (reaction on solid is negative)
            Vector3r cell_force = Vector3r::Zero();
            for (int i = 0; i < Q; ++i) {
                cell_force += omega_s[i] * C->lattice->discrete_velocity[i];
            }
            cell_force *= Bn;

            B->state->hydro_force -= cell_force;

            // 2D torque about body center (z-component)
            const Vector3r r = C->state->pos - B->state->pos;
            B->state->hydro_torque -= (r[0] * cell_force[1] - r[1] * cell_force[0]);
        }
    }
}


void ContactResolution::action(){
    Scene& S = Scene::get_Scene();

    int bodySize = (int)S.bodies.size();
	ASSERT(bodySize > 0);
	//Brute force method:
	for (int i = 0; i < bodySize - 1; ++i) {
		for (int j = i + 1; j < bodySize; ++j) {
			if (S.bodies[i]->CheckInteraction(S.bodies[j]->id))	continue;
			if ((S.bodies[i]->state->pos - S.bodies[j]->state->pos).norm() < S.bodies[i]->shape->radius + S.bodies[j]->shape->radius) {
				S.interactions.push_back(std::make_shared<Interaction>(S.bodies[i], S.bodies[j]));
				S.bodies[i]->inter.push_back(S.interactions[S.interactions.size() - 1]);
				S.bodies[j]->inter.push_back(S.interactions[S.interactions.size() - 1]);
			}
		}
	}
}

void InteractionLoop::action(){
    Scene& S = Scene::get_Scene();
    
    for (auto& I : S.interactions) {
		I->CalculateUnitVectorandContact();
		I->CalculateForceAndShearIncrements();
		I->ApplyFrictionLaw();
	}
}

void BodyLoop::action(){
    Scene& S = Scene::get_Scene();

	//Force calculation:
	for (auto& B : S.bodies) {
		//Force Resetter:
		B->state->force = Vector3r::Zero();
		B->state->moment = 0.0;

		//Contact force
		for (auto& I : B->inter) {
			auto Il = I.lock();
			auto iBody1 = Il->body1.lock();
			auto iBody2 = Il->body2.lock();
			if (B->id == iBody1->id) {
				B->state->force += -Il->normal_force * Il->unit_normal_vector;
				B->state->force += -Il->shear_force * Il->unit_shear_vector;
				B->state->moment += -Il->shear_force * iBody1->shape->radius;
			}
			if (B->id == iBody2->id) {
				B->state->force += Il->normal_force * Il->unit_normal_vector;
				B->state->force += Il->shear_force * Il->unit_shear_vector;
				B->state->moment += -Il->shear_force * iBody2->shape->radius;
			}
		}

        Vector3r unitUp = { 0, 1, 0 };
        Vector3r unitDown = { 0, -1, 0 };
        Vector3r unitRight = { 1, 0, 0 };
        Vector3r unitLeft = { -1, 0, 0 };
        Vector3r force = Vector3r::Zero();
        if ((B->state->pos[0] - B->shape->radius) < S.model_min_corner[0])  force += S.border_stiffness * abs(B->shape->radius - B->state->pos[0]) * unitRight;
        if ((B->state->pos[0] + B->shape->radius) > S.model_max_corner[0])  force += S.border_stiffness * abs(B->shape->radius - (S.model_max_corner[0] - B->state->pos[0])) * unitLeft;
        if ((B->state->pos[1] - B->shape->radius) < S.model_min_corner[1])  force += S.border_stiffness * abs(B->shape->radius - B->state->pos[1]) * unitUp;
        if ((B->state->pos[1] + B->shape->radius) > S.model_max_corner[1])  force += S.border_stiffness * abs(B->shape->radius - (S.model_max_corner[1] - B->state->pos[1])) * unitDown;
    
		//Body, fluid and border force:
        B->state->force += force;
		B->state->force += B->material->mass * S.gravity;
		
		// For IMB coupling
		if (S.dem_coupling == "IMB") {
           
            // Buoyancy force
			double V_particle = M_PI * B->shape->radius * B->shape->radius;
			B->state->force -= S.initial_density * V_particle * S.gravity;
			
			// Hydrodynamic force and torque
			B->state->force += B->state->hydro_force;
            B->state->moment += B->state->hydro_torque;
		}
	}
}

void Integrator::action(){
    Scene& S = Scene::get_Scene();

    for (auto& B : S.bodies) {

		//Calculate accelaration from forces:
		Vector3r  linAccel = Vector3r::Zero();
		Vector3r  f = B->state->force;
		double    m = B->state->moment;
		double    rotAccel = 0.0;
		int       signV = 0;
		int       signM = 0;

		for (int i = 0; i < linAccel.size(); ++i) {
			B->state->vel[i] > 0 ? signV = 1 : signV = -1;
			linAccel[i] = ((f[i] - S.local_damping * abs(f[i]) * signV) / B->material->mass) * B->blockedDOFs[i];
		}

		B->state->rotVel > 0 ? signM = 1 : signM = -1;
		rotAccel = ((m - S.local_damping * abs(m) * signM) / B->shape->inertia_moment) * B->blockedMDOFs;

		//Update velocity and position (LeapFrog method):
		if (S.nIter == 0) {
			B->state->vel += linAccel * S.dt_dem * 0.5;
			B->state->rotVel += rotAccel * S.dt_dem * 0.5;
		}
		else {
			B->state->vel += linAccel * S.dt_dem;
			B->state->rotVel += rotAccel * S.dt_dem;
		}
		B->state->pos += B->state->vel * S.dt_dem;
		B->state->rot += B->state->rotVel * S.dt_dem;
	}
	++S.nIter;
	ASSERT(S.nIter > 0);

}

void UpdateContact::action(){
    Scene& S = Scene::get_Scene();

    S.interactions.erase(std::remove_if(std::begin(S.interactions), std::end(S.interactions), [](std::shared_ptr<Interaction> I) {return !I->CheckContact(); }), std::end(S.interactions));
	for (auto& B : S.bodies) {
		B->inter.erase(std::remove_if(std::begin(B->inter), std::end(B->inter), [](std::weak_ptr<Interaction> I) {return I.expired(); }), std::end(B->inter));
	}
}



