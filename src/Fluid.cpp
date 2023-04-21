//
// Created by Wentinn Liao on 4/12/23.
//

#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9);

/** TODO: add fluid density maybe */
Fluid::Fluid(int length, int width, int height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params) {
    double volume = length * width * height;
    double true_num_particles = volume * params.density * 6.022E23 / params.molar_mass;
    this->SMOOTHING_RADIUS = params.average_distance * std::cbrt(true_num_particles / nParticles);

    grid = new std::vector<PointMass>[LENGTH * WIDTH * HEIGHT];

    // make planes
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(1, 0, 0), 0.3));
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(0, 1, 0), 0.3));
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(0, 0, 1), 0.3));
    this->collisionObjects.push_back(new Plane(Vector3D(LENGTH, 0, 0), Vector3D(-1, 0, 0), 0.3));
    this->collisionObjects.push_back(new Plane(Vector3D(0, WIDTH, 0), Vector3D(0, -1, 0), 0.3));

    // making particles
    // using code from https://stackoverflow.com/questions/38244877/how-to-use-stdnormal-distribution
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device random_device;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(random_device());
    // std of water velocity is approx 642.50364 m/s
    std::normal_distribution<float> norm_dist_gen(0, 642.503643481);
    for (int i = 0; i < nParticles; i += 1) {
        Vector3D position = Vector3D(uniform(0, LENGTH), uniform(0, WIDTH), uniform(0, 0.5 * HEIGHT));
        Vector3D velocity = Vector3D(norm_dist_gen(gen), norm_dist_gen(gen), norm_dist_gen(gen));
        PointMass pt = PointMass(position, velocity, false);

        this->get_position(position).push_back(pt);
    }
}

std::vector<PointMass> &Fluid::get_position(const Vector3D &pos) {
    return grid[(int) pos[2] + HEIGHT * ((int) pos[1] + WIDTH * (int) pos[0])];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;
    double mass = (LENGTH * WIDTH * HEIGHT * density) / NUM_PARTICLES;

    Vector3D total_external_force(0);
    for (const Vector3D &acc: external_accelerations) {
        total_external_force += acc;
    total_external_force *= mass;

    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            pt.forces += total_external_force;
        }
    }

    /** Take out this part, should be handled in Verlet integration */
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            pt.velocity += total_acc;
        }
    }

    /** handle collisions **/
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            for (CollisionObject *co: this->collisionObjects) {
                co->collide(pt);
            }
        }
    }

    /** TODO: handle self collisions */

    /** update position **/
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            pt.position = pt.position + pt.velocity * delta_t;
        }
    }
}

void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
}
