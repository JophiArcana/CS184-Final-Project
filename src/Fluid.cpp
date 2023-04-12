//
// Created by Wentinn Liao on 4/12/23.
//

#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"



Fluid::Fluid(int length, int width, int height, int nParticles) {
    this->LENGTH = length;
    this->WIDTH = width;
    this->HEIGHT = height;
    this->NUM_PARTICLES = nParticles;
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
    std::normal_distribution<float> norm_dist_gen(0, 2);
    for (int i = 0; i < nParticles; i += 1) {
        Vector3D position = Vector3D(uniform(0, LENGTH), uniform(0, WIDTH), uniform(0, 0.5 * HEIGHT));
        Vector3D velocity = Vector3D(norm_dist_gen(gen), norm_dist_gen(gen), norm_dist_gen(gen));
        PointMass pt = PointMass(position, velocity, false);

        grid[int(position[2]) + HEIGHT * int(position[1]) + WIDTH * HEIGHT * int(position[0])].push_back(pt);
    }
}

std::vector<PointMass> &Fluid::get_position(int x, int y, int z) {
    return grid[z + HEIGHT * y + WIDTH * HEIGHT * x];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, std::vector<Vector3D> external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass pt: this->grid[index]) {
            for (Vector3D acc: external_accelerations) {
                pt.velocity += (acc * delta_t);
            }
        }
    }

    /** handle collisions **/
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass pt: this->grid[index]) {
            for (CollisionObject *co: this->collisionObjects) {
                co->collide(pt);
            }
        }
    }

    /** TODO: handle self collisions */

    /** update position **/
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass pt: this->grid[index]) {
            pt.position = pt.position + pt.velocity * delta_t;
        }
    }
}

void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
}
