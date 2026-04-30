#include <fstream>
#include "scene/scene.h"
#include "dem/body.h"
#include "property/state.h"
#include "property/shape.h"
#include "scene/output.h"

void Output::DisplaySimulationInfo(){
	Scene& S = Scene::get_Scene();
	int ignore = std::system("clear");
	std::cout << "----------------------- 2D LBM/DEM Simulation  ---------------------------------" << "\n";
	std::cout << "Iteration Number: " << S.iter                                                     << "\n";
	std::cout << "Number of Cells : " << S.cells.size()                                             << "\n";
	std::cout << "Number of Bodies: " << S.bodies.size()                                            << "\n";
	std::cout << "----------------------- LBM Parameters ----------------------------------------"  << "\n";
	std::cout << "Time Step       : " << S.dt_lbm                                                   << "\n";
	std::cout << "Lattice Spacing : " << S.lattice_spacing                                          << "\n";
	std::cout << "Relaxation Time : " << S.relaxation_time                                          << "\n";
    std::cout << "Guo Forcing     : " << S.GUO_fluid_forcing[0] << ", " << S.GUO_fluid_forcing[1]   << "\n";
	std::cout << "-------------------------- DEM Parameters --------------------------------------" << "\n";
    std::cout << "Number of bodies: " << S.bodies.size()                                            << "\n";
	std::cout << "Time Step       : " << S.dt_dem                                                   << "\n";
	std::cout << "Friction Angle  : " << S.friction_angle                                           << "\n";
	std::cout << "Shear Stiffness : " << S.shear_stiffness                                          << "\n";
	std::cout << "Normal Stiffness: " << S.normal_stiffness                                         << "\n";
}

void Output::ExportFluidVtk(std::string file_prefix) {
    Scene& S = Scene::get_Scene();

    const int Nx = static_cast<int>(S.model_max_corner[0]);
    const int Ny = static_cast<int>(S.model_max_corner[1]);
    const int Nz = static_cast<int>(S.model_max_corner[2]);
    const int Npts = Nx * Ny * Nz;

    const double cs2 = 1.0 / 3.0;

    std::ofstream out(file_prefix + std::to_string(S.iter) + ".vtk");
    out << std::fixed << std::setprecision(6);

    // --- Header ---
    out << "# vtk DataFile Version 3.0\n";
    out << "Fluid state\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << Npts << "\n";

    // --- Geometry mask: 1=fluid, 0=solid/wall ---
    out << "SCALARS 1.Geometry float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (auto& Cptr : S.cells) {
        const auto& C = *Cptr;
        float geometry = (C.is_solid || C.is_wall) ? 0.0f : 1.0f;
        out << geometry << "\n";
    }

    // --- Density ---
    out << "SCALARS 2.Density float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (auto& Cptr : S.cells) {
        const auto& C = *Cptr;
        out << static_cast<float>(C.material->density) << "\n";
    }

    // --- Pressure = cs^2 * rho ---
    out << "SCALARS 3.Pressure float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (auto& Cptr : S.cells) {
        const auto& C = *Cptr;
        const float p = static_cast<float>(cs2 * C.material->density);
        out << p << "\n";
    }

    // --- Velocity vector ---
    out << "VECTORS 4.Velocity float\n";

    for (auto& Cptr : S.cells) {
        const auto& C = *Cptr;
        if (C.is_solid || C.is_wall || C.solid_fraction == 1.0) {
            out << "0 0 0\n"; 
        } else {
            out << static_cast<float>(C.state->vel[0]) << " "
                << static_cast<float>(C.state->vel[1]) << " "
                << static_cast<float>(C.state->vel[2]) << "\n";
        }
    }

    // --- Geometry mask: 1=fluid, 0=solid/wall ---
    out << "SCALARS 5.SolidFraction float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (auto& Cptr : S.cells) {
        const auto& C = *Cptr;
        out << C.solid_fraction << "\n";
    }
    out.close();
}

void Output::ExportParticleVtk(std::string _file_name) {
    Scene& S = Scene::get_Scene();

    // create output file name: e.g., dem_0001.vtk
    std::ostringstream fname;
    fname << _file_name << std::setfill('0') << std::setw(4) << S.iter << ".vtk";

    std::ofstream out(fname.str());
    if (!out.is_open()) {
        std::cerr << "Error: cannot open " << fname.str() << " for writing\n";
        return;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "DEM particle output\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    const size_t N = S.bodies.size();
    out << "POINTS " << N << " double\n";
    for (const auto& B : S.bodies) {
        const auto& p = B->state->pos;
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    out << "\nPOINT_DATA " << N << "\n";

    // --- Radius ---
    out << "SCALARS Radius double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& B : S.bodies)
        out << B->shape->radius << "\n";

    // --- Velocity ---
    out << "\nVECTORS Velocity double\n";
    for (const auto& B : S.bodies) {
        const auto& v = B->state->vel;
        out << v.x() << " " << v.y() << " " << v.z() << "\n";
    }

    // --- Force ---
    out << "\nVECTORS Force double\n";
    for (const auto& B : S.bodies) {
        const auto& f = B->state->force;
        out << f.x() << " " << f.y() << " " << f.z() << "\n";
    }

    out.close();
}

void Output::ExportFluidCsv(Vector3r _pos, std::string _file_name){
	Scene& S = Scene::get_Scene();
	int id = S.CalculateCellId(_pos[0], _pos[1]);
	std::ofstream out;
	out.open(_file_name + std::to_string(S.iter) + ".csv", std::ios_base::app);
	out << S.cells[id]->state->pos[0]     << ", " << S.cells[id]->state->pos[1]     << ", " << S.cells[id]->state->pos[2] << ", " 
	    << S.cells[id]->state->vel[0]     << ", " << S.cells[id]->state->vel[1]     << ", " << S.cells[id]->state->vel[2] << ", " 
		<< S.cells[id]->state->vel.norm() << ", " << S.cells[id]->material->density << "\n";
	out.close();
}

void Output::ExportVelocityProfileCsv(const std::string& file_name) {
    Scene& S = Scene::get_Scene();

    int Nx = static_cast<int>(S.model_max_corner[0]);
    int Ny = static_cast<int>(S.model_max_corner[1]);

    std::ofstream out(file_name + "_" + std::to_string(S.iter) + ".csv");
    out << "y, ux_avg\n";

    for (int j = 0; j < Ny; ++j) {
        double ux_sum = 0.0;
        int count = 0;
        for (int i = 0; i < Nx; ++i) {
            auto& C = *S.cells[S.CalculateCellId(i, j)];
            if (C.is_solid || C.is_wall) continue;
            ux_sum += C.state->vel[0];
            ++count;
        }
        if (count > 0) {
            double ux_avg = ux_sum / count;
            out << j << "," << ux_avg << "\n";
        }
    }
    out.close();
}


void Output::ExportParticleCsv(std::string _file_name){
	Scene& S = Scene::get_Scene();
	std::ofstream out;
	out.open(_file_name, std::ios_base::app);
	if (S.bodies.size() == 1){
		out << S.bodies[0]->state->pos[1]<< ", " << S.bodies[0]->material->density * S.bodies[0]->state->pos[1] * (9.81) << ", " << S.bodies[0]->material->density * (S.bodies[0]->state->vel[1] * S.bodies[0]->state->vel[1]) / 2 << "\n";
	} else {
		out << S.bodies[0]->state->pos[0] << ", " << S.bodies[0]->material->density * (S.bodies[0]->state->vel[0] * S.bodies[0]->state->vel[0]) / 2 << ", " 
		    << S.bodies[1]->state->pos[0] << ", " << S.bodies[1]->material->density * (S.bodies[1]->state->vel[0] * S.bodies[1]->state->vel[0]) / 2 << "\n";
	}
	out.close();
}