#include "scene/Scene.h"


static Scene scene;

Scene& Scene::get_Scene(){ return scene;}

void Scene::AddRectangularCanal(){
    Scene& S = Scene::get_Scene();
    int id;
    
     //Add cells:
    for (int j = 0; j < model_max_corner[1]; ++j)
    for (int i = 0; i < model_max_corner[0]; ++i){
        id = CalculateCellId(i, j);
        cells.push_back(std::make_shared<Cell>(id, std::make_shared<Lattice>(), std::make_shared<Material>(), std::make_shared<State>()));
        cells[id]->state->pos = Vector3r(i, j, 0);
        cells[id]->set_neighbor_node();
        double h = lattice_spacing * 0.5;

        //Setting boost geometry polygon
        Vector3r A = Vector3r(i,j,0) + h * Vector3r(-1,-1, 0);
        Vector3r B = Vector3r(i,j,0) + h * Vector3r(-1, 1, 0);
        Vector3r C = Vector3r(i,j,0) + h * Vector3r( 1, 1, 0);
        Vector3r D = Vector3r(i,j,0) + h * Vector3r( 1,-1, 0);

        cells[id]->grid.outer().push_back(point(A[0], A[1]));
        cells[id]->grid.outer().push_back(point(B[0], B[1]));
        cells[id]->grid.outer().push_back(point(C[0], C[1]));
        cells[id]->grid.outer().push_back(point(D[0], D[1]));
        cells[id]->grid.outer().push_back(point(A[0], A[1]));

        boost::geometry::correct(cells[id]->grid);
        }

        // Lattice search initialization if body present
        if (S.bodies.size() > 0 && S.dem_coupling == "IMB") {
            S.engines.push_back(std::make_shared<LatticeSearch>());
            for (auto& E : S.engines){
                E->action();
            }
        }

        //Define walls:
        //Top:
        for (int i = 0; i < model_max_corner[0]; ++i){
            id = CalculateCellId(i, model_max_corner[1]-1);
            cells[id]->is_wall = true;
            if (S.dem_coupling == "IMB") cells[id]->solid_fraction = 1.0;
        }
        //Bot:
        for (int i = 0; i < model_max_corner[0]; ++i){
            id = CalculateCellId(i, 0);
            cells[id]->is_wall = true;
            if (S.dem_coupling == "IMB") cells[id]->solid_fraction = 1.0;
        }
        
        for (auto& C : cells){
            C->set_initial_condition();
            if (S.collision_operator == "MRT")  C->set_mrt_parameters();
        }
}

void Scene::addParticle(std::string _fileName, double _density, Vector3r _velocity) {
	int id = bodies.size();
	bodies.push_back(std::make_shared<Body>(id, std::make_shared<Shape>(), std::make_shared<State>(), std::make_shared<Material>()));
    bodies[id]->shape->generate_particle(_fileName);
    Vector3r sum = Vector3r::Zero();
    for(auto c : bodies[id]->shape->coordinates){
        sum += c;
    }
    bodies[id]->state->pos = sum / bodies[id]->shape->coordinates.size();
    bodies[id]->shape->radius = (bodies[id]->state->pos - bodies[id]->shape->coordinates[0]).norm();
    bodies[id]->material->density = _density;
    bodies[id]->state->vel = _velocity;
    bodies[id]->material->mass = _density * M_PI * bodies[id]->shape->radius * bodies[id]->shape->radius;
    bodies[id]->shape->inertia_moment = bodies[id]->material->density * bodies[id]->shape->radius * bodies[id]->shape->radius * 0.5;
}

void Scene::AddDisk(Vector3r _center, double _radius, double _density){
    int id = bodies.size();
    bodies.push_back(std::make_shared<Body>(id, std::make_shared<Shape>(), std::make_shared<State>(), std::make_shared<Material>()));
    
    // Build cloud polygon and coordinates vector
    bodies[id]->shape->cloud.clear();
    bodies[id]->shape->coordinates.clear();
    
    for (double theta = 0.0; theta < 2.0 * M_PI; theta += 0.1) {
        double x = _center[0] + _radius * std::cos(theta);
        double y = _center[1] + _radius * std::sin(theta);
        boost::geometry::append(bodies[id]->shape->cloud.outer(), point(x, y));
        bodies[id]->shape->coordinates.emplace_back(Vector3r(x, y, 0.0));
    }
    boost::geometry::correct(bodies[id]->shape->cloud);
    
    // Set body properties
    bodies[id]->state->pos = _center;
    bodies[id]->shape->radius = _radius;
    
    // Material properties (2D disk with unit thickness)
    bodies[id]->material->density = _density;
    bodies[id]->material->mass = _density * M_PI * _radius * _radius;  // mass = ρ × π × R²
    bodies[id]->shape->inertia_moment = 0.5 * bodies[id]->material->mass * _radius * _radius;  // I = 0.5 × m × R²
}

void Scene::MoveToNextTimeStep(){
    for(auto E : engines){
        E->action();
    }
    time += dt_lbm;
    ++iter;
}

void Scene::SetRelaxationParamters(double _s1, double _s2, double _s3, double _s4, double _s5, double _s6, double _s7, double _s8, double _s9){
    relaxation_parameters  << _s1, _s2, _s3, _s4, _s5, _s6, _s7, _s8, _s9;
}

void Scene::addDiskPacking(double _Rmin, double _Rmax, double _wall_gap, double _inlet_buf, double _outlet_buf, double _non_overlap, 
                            int _max_particles, int _wall_particles, double _rho_p, double _Rw_min, double _Rw_max){

    // Get the number of cells in the x and y directions
    int Nx = static_cast<int>(model_max_corner[0]);
    int Ny = static_cast<int>(model_max_corner[1]);
    
     // Particle sizes and number of particles
     double Rmin = _Rmin;
     double Rmax = _Rmax;
     int max_particles = _max_particles;
 
     // Buffers and gaps
     double wall_gap   = _wall_gap;
     double inlet_buf  = _inlet_buf;
     double outlet_buf = _outlet_buf;
     double x_min_p = inlet_buf  + Rmax + wall_gap;
     double x_max_p = Nx - outlet_buf - Rmax - wall_gap;
     double y_min_p = Rmax + wall_gap;
     double y_max_p = Ny - Rmax - wall_gap;
     double non_overlap = _non_overlap;

     //Check overlap
     auto overlap_ok = [&](const Vector3r& pos, double Rnew){
        for (auto& B : bodies){
            const double dist = (pos - B->state->pos).head<2>().norm();
            if (dist < non_overlap * (Rnew + B->shape->radius))
                return false;
        }
        return true;
    };

    // Random packing
    std::mt19937 gen(123);
    std::uniform_real_distribution<double> urx(x_min_p, x_max_p);
    std::uniform_real_distribution<double> ury(y_min_p, y_max_p);
    std::uniform_real_distribution<double> urr(Rmin, Rmax);

    for (int n = 0; n < max_particles; ++n){
        bool placed = false;
        for (int trial = 0; trial < 800 && !placed; ++trial){
            const double R = urr(gen);
            Vector3r pos(urx(gen), ury(gen), 0.0);
            if (overlap_ok(pos, R)) {
                AddDisk(pos, R, _rho_p);
                placed = true;
            }
        }
    }

    // Extra small disks near top/bottom to reduce preferential channels
    double Rw_min = _Rw_min;
    double Rw_max = _Rw_max;

    std::uniform_real_distribution<double> urx_w(x_min_p, x_max_p);
    std::uniform_real_distribution<double> ury_bottom(Rw_max + 0.1, 3.0 * Rw_max);
    std::uniform_real_distribution<double> ury_top(Ny - 3.0 * Rw_max, Ny - (Rw_max + 0.1));
    std::uniform_real_distribution<double> urr_w(Rw_min, Rw_max);

    for (int n = 0; n < _wall_particles; ++n){
        bool placed = false;
        for (int trial = 0; trial < 800 && !placed; ++trial){
            const double R = urr_w(gen);
            const bool top = (n % 2 == 0);
            const double y = top ? ury_top(gen) : ury_bottom(gen);
            Vector3r pos(urx_w(gen), y, 0.0);
            if (overlap_ok(pos, R)) {
                AddDisk(pos, R, _rho_p);
                placed = true;
            }
        }
    }
}