//
// Created by Wentinn Liao on 4/12/23.
//

#ifndef CLOTHSIM_FLUID_H
#define CLOTHSIM_FLUID_H

#include <random>
#include "pointMass.h"
#include "./collision/collisionObject.h"
#include "FluidMesh.h"
#include "CGL/CGL.h"

#define uniform(a, b)   (a + (double) (b - a) * std::rand() / RAND_MAX)


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
    FluidMesh *mesh;

    void buildFluidMesh();

};


#endif //CLOTHSIM_FLUID_H


