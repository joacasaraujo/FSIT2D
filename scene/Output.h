#ifndef OUTPUT_H
#define OUTPUT_H

#include "maths/Math.h"


class Scene;

class Output {
    public:
    void DisplaySimulationInfo();
    
    void ExportFluidVtk(std::string _file_name);
    void ExportParticleVtk(std::string _file_name);

    void ExportFluidCsv(Vector3r _pos, std::string _file_name);
    void ExportVelocityProfileCsv(const std::string& file_name);
    void ExportParticleCsv(std::string _file_name);
};

#endif