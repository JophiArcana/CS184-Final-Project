//
// Created by Wentinn Liao on 4/12/23.
//

#ifndef CLOTHSIM_FLUID_H
#define CLOTHSIM_FLUID_H

#define MULTITHREAD true

#include <random>
#include <deque>
#include "pointMass.h"
#include "./collision/collisionObject.h"
#include "FluidMesh.h"
#include "CGL/CGL.h"

// #define uniform(a, b)   (a + (double) (b - a) * std::rand() / RAND_MAX)

enum e_orientation {
    HORIZONTAL = 0, VERTICAL = 1
};

struct FluidParameters {
    FluidParameters() = default;

    FluidParameters(double density,
                    double rms_velocity,
                    double average_distance,
                    double molar_mass,
                    double kinematic_viscosity,
                    double tait_coefficient,
                    double tait_gamma)
            : density(density), rms_velocity(rms_velocity), average_distance(average_distance), molar_mass(molar_mass),
              kinematic_viscosity(kinematic_viscosity), tait_coefficient(tait_coefficient), tait_gamma(tait_gamma) {}

    ~FluidParameters() = default;

    double density;                 // kg/m^3
    double rms_velocity;            // m/s
    double average_distance;        // m
    double molar_mass;              // kg/mol
    double kinematic_viscosity;     // m^2/s
    double tait_coefficient;        // kg/ms^2
    double tait_gamma;              // 1
};

class Fluid {
public:
    static const FluidParameters WATER;

    Fluid();

    Fluid(FluidParameters params);

    Fluid(double length, double width, double height, int nParticles, FluidParameters params);

    std::vector<PointMass *> list;
    std::vector<PointMass *> *grid() const;
    int get_index(const Vector3D &pos) const;
    std::vector<PointMass *> &get_cell(const Vector3D &pos) const;
    std::vector<std::vector<int>> global_neighbor_indices;

    void simulate(double frames_per_sec, double simulation_steps,
                  const std::vector<Vector3D> &external_accelerations);

    void forward_movement(const std::vector<Vector3D> &external_accelerations, double delta_t);
    void incompressibility_adjustment(int n_iter, double step_size);

    void collision_update(double delta_t);
    void constrain_update();
    void cell_update();

    std::vector<PointMass *> *grid1, *grid2;
    bool grid_toggle;

    // although length, etc. are constant through the simulation, not "const" to assign during creation
    double LENGTH, WIDTH, HEIGHT;
    double CELL_SIZE;
    int G_LENGTH, G_WIDTH, G_HEIGHT, G_SIZE;
    int NUM_PARTICLES;
    FluidParameters PARAMS;
    double VOLUME, SMOOTHING_RADIUS, PARTICLE_MASS;


    std::vector<CollisionObject *> collisionObjects;

    FluidMesh *mesh;

    void buildFluidMesh();
    void debugFluidMesh();

    std::vector<double> timestamps;
    e_orientation orientation;

    double W(PointMass *pi, PointMass *pj) const;
    Vector3D scaled_grad_W(PointMass *pi, PointMass *pj) const;

    std::vector<std::vector<double>> batch_W(int index) const;
    std::vector<std::vector<Vector3D>> batch_scaled_grad_W(int index) const;

    std::vector<double> batch_density(const std::vector<std::vector<double>> &W) const;

    std::vector<std::vector<Vector3D>>
    global_lambda_displacement(const std::vector<std::vector<double>> &global_density,
                               const std::vector<std::vector<std::vector<double>>> &global_W,
                               const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const;
    std::vector<std::vector<double>> batch_scorr(const std::vector<std::vector<double>> &W) const;

    std::vector<std::vector<Vector3D>>
    global_curl_velocity(const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const;
    std::vector<std::vector<Vector3D>>
    global_normalized_grad_norm_curl_velocity(const std::vector<std::vector<double>> &global_normalized_curl_velocity_norm,
                                              const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const;

private:
    double KERNEL_COEFF, GRAD_KERNEL_COEFF;
    double SELF_KERNEL;
    double SCORR_COEFF;
    double SELF_SCORR;
};


#endif //CLOTHSIM_FLUID_H


