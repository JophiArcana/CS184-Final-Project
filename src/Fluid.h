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

struct FluidParameters {
    FluidParameters() {};

    FluidParameters(double density,
                    double rms_velocity,
                    double average_distance,
                    double molar_mass)
            : density(density), rms_velocity(rms_velocity), average_distance(average_distance), molar_mass(molar_mass) {}

    ~FluidParameters() {}

    double density;             // kg/m^3
    double rms_velocity;        // m/s
    double average_distance;    // m
    double molar_mass;          // kg/mol
};

class Fluid {
public:
    static const FluidParameters WATER;

    Fluid(double length, double width, double height, int nParticles, FluidParameters params);

    std::vector<PointMass> &get_position(const Vector3D &pos);

    void simulate(double frames_per_sec, double simulation_steps,
                  const std::vector<Vector3D> &external_accelerations);

    std::vector<PointMass> *grid;

    const double LENGTH, WIDTH, HEIGHT;
    int G_LENGTH, G_WIDTH, G_HEIGHT;
    const int NUM_PARTICLES;
    const FluidParameters PARAMS;
    double SMOOTHING_RADIUS;
    double PARTICLE_MASS;

    std::vector<CollisionObject *> collisionObjects;
    FluidMesh *mesh;

    double W(const PointMass &p, const PointMass &q);
    std::vector<double> batch_density(int index);
    std::vector<Vector3D> batch_kernel_gradient(int index);

    std::vector<Vector3D> pressure_gradient(int index);

    std::vector<Vector3D> velocity_laplacian(int index);

    void buildFluidMesh();

private:
    double KERNEL_COEFF;
};


#endif //CLOTHSIM_FLUID_H


