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

// #define uniform(a, b)   (a + (double) (b - a) * std::rand() / RAND_MAX)

struct FluidParameters {
    FluidParameters() = default;

    FluidParameters(double density,
                    double rms_velocity,
                    double average_distance,
                    double molar_mass,
                    double kinematic_viscosity,
                    double tait_coefficient)
            : density(density), rms_velocity(rms_velocity), average_distance(average_distance), molar_mass(molar_mass), kinematic_viscosity(kinematic_viscosity), tait_coefficient(tait_coefficient) {}

    ~FluidParameters() = default;

    double density;                 // kg/m^3
    double rms_velocity;            // m/s
    double average_distance;        // m
    double molar_mass;              // kg/mol
    double kinematic_viscosity;     // m^2/s
    double tait_coefficient;        // kg/ms^2
};

class Fluid {
public:
    static const FluidParameters WATER;

    Fluid(double length, double width, double height, int nParticles, FluidParameters params);

    std::vector<PointMass *> &get_position(const Vector3D &pos) const;

    void simulate(double frames_per_sec, double simulation_steps,
                  const std::vector<Vector3D> &external_accelerations);

    std::vector<PointMass *> *grid;

    const double LENGTH, WIDTH, HEIGHT;
    int G_LENGTH, G_WIDTH, G_HEIGHT;
    const int NUM_PARTICLES;
    const FluidParameters PARAMS;
    double SMOOTHING_RADIUS;
    double PARTICLE_MASS;

    std::vector<CollisionObject *> collisionObjects;

    FluidMesh *mesh;
    void buildFluidMesh();

    std::vector<double> timestamps;

    double W(PointMass *pi, PointMass *pj) const;
    Vector3D unnormalized_grad_W(PointMass *pi, PointMass *pj) const;

    std::vector<std::vector<double>> batch_W(int index) const;
    std::vector<std::vector<Vector3D>> batch_unnormalized_grad_W(int index) const;

    std::vector<double> batch_density(const std::vector<std::vector<double>> &W) const;
    std::vector<double> batch_pressure(const std::vector<double> &density) const;
    std::vector<Vector3D> batch_scaled_grad_pressure(const std::vector<double> &pressure, const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const;
    std::vector<Vector3D> batch_scaled_laplacian_velocity(int index, const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const;

private:
    double KERNEL_COEFF;
    double SELF_KERNEL;
};


#endif //CLOTHSIM_FLUID_H


