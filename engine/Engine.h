#ifndef ENGINE_H
#define ENGINE_H

#include "maths/Math.h"
class Scene;

class Engine {
    public:
        virtual void action() = 0;
};

class ZouAndHeBC : public Engine {
    public:
        virtual void action() override;
};

class FluidCollision : public Engine {
    public:
        virtual void action() override;
};

class FluidStreaming : public Engine {
    public:
        virtual void action() override;
};

class ApplyFluidForcing : public Engine {
    public:
        virtual void action() override;
};

class FluidParticleForce : public Engine {
    public:
        virtual void action() override;
};

class ImbBoundary : public Engine {
    public:
        virtual void action() override;
};

class LatticeSearch : public Engine {
    public:
        virtual void action() override;
};

class NeighborSearch : public Engine {
    public:
    virtual void action() override;
};

class ContactResolution : public Engine {
    public:
    virtual void action() override;
};

class InteractionLoop : public Engine {
    public:
    virtual void action() override;
};

class BodyLoop : public Engine {
    public:
    virtual void action() override;
};

class Integrator : public Engine {
    public:
    virtual void action() override;
};

class UpdateContact : public Engine {
    public:
    virtual void action() override;
};

#endif