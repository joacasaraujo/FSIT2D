#include "Shape.h"

void Shape::generate_particle(std::string _file_name){
    //Atribuição dos pontos das coordenadas
    std::ifstream file(_file_name);
    double x,y,z;
    while(file >> x >> y >> z){
        coordinates.emplace_back(Vector3r(x,y,z));
    }

    for(int i = 0; i < coordinates.size(); ++i){
        cloud.outer().push_back(point(coordinates[i][0],coordinates[i][1]));
        if(i == coordinates.size()-1){
            cloud.outer().push_back(point(coordinates[0][0],coordinates[0][1]));
        }
    }

    //Geração das conectividades
    for(int i = 0; i < coordinates.size(); ++i){
        if(i == coordinates.size()-1){
            conectivities.emplace_back(Vector3r(coordinates.size()-1,0,0));
        } else {
            conectivities.emplace_back(Vector3r(i,i+1,0));
        }
    }
    boost::geometry::correct(cloud);
    file.close();
    
}