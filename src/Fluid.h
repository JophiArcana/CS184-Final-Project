//
// Created by Wentinn Liao on 4/12/23.
//

#ifndef CLOTHSIM_FLUID_H
#define CLOTHSIM_FLUID_H

#include "pointMass.h"
#include "./collision/collisionObject.h"
#include "CGL/CGL.h"


class Fluid {
    Fluid(int length, int width, int height, int nParticles);

    std::vector<PointMass> &get_position(int x, int y, int z);

    void simulate(double frames_per_sec, double simulation_steps,
                  std::vector<Vector3D> external_accelerations);

    std::vector<PointMass> *grid;

    int LENGTH, WIDTH, HEIGHT;
    int NUM_PARTICLES;
    double PARTICLE_RADIUS = 0.05;

    std::vector<CollisionObject *> collisionObjects;

    void buildFluidMesh();

};


#endif //CLOTHSIM_FLUID_H


