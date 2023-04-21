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
    FluidParameters() = default;

    FluidParameters(double density,
                    double rms_velocity,
                    double average_distance,
                    double molar_mass)
            : density(density), rms_velocity(rms_velocity), average_distance(average_distance), molar_mass(molar_mass) {}

    ~FluidParameters() = default;

    double density;             // kg/m^3
    double rms_velocity;        // m/s
    double average_distance;    // m
    double molar_mass;          // kg/mol
};

class Fluid {
public:
    static const FluidParameters WATER;

    Fluid(double length, double width, double height, int nParticles, FluidParameters params);

    std::vector<PointMass> &get_position(const Vector3D &pos) const;

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

    double W(const PointMass &pi, const PointMass &pj) const;
    Vector3D grad_W(const PointMass &pi, const PointMass &pj) const;

    std::vector<double> batch_density(int index) const;
    std::vector<Vector3D> batch_grad_pressure(int index) const;
    std::vector<Vector3D> batch_laplacian_velocity(int index) const;

    void buildFluidMesh();

private:
    double KERNEL_COEFF;
    double SELF_KERNEL;
};


#endif //CLOTHSIM_FLUID_H


