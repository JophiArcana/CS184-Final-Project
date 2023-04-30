#ifndef POINTMASS_H
#define POINTMASS_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

// Forward declarations
class Halfedge;

struct PointMass {
    PointMass(Vector3D position, Vector3D velocity, bool pinned)
            : pinned(pinned), start_position(position), position(position), tentative_position(position),
              velocity(velocity), acceleration(0), stage(true) {}

    Vector3D normal();

//    Vector3D velocity(double delta_t) {
//        return (position - last_position) / delta_t;
//    }

    // static values
    bool pinned;
    Vector3D start_position;

    // dynamic values
    Vector3D position, tentative_position;
    Vector3D velocity;
    Vector3D forces;
    Vector3D acceleration; // TODO: verify if unnecessary?

    // for cell updates
    bool stage;

    // mesh reference
    Halfedge *halfedge;
};

#endif /* POINTMASS_H */
