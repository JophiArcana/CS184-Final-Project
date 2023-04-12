//
// Created by Wentinn Liao on 4/12/23.
//

#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

#define uniform(a, b)   (a + (double) (b - a) * std::rand() / RAND_MAX)

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
    for (int i = 0; i < nParticles; i += 1) {
        Vector3D position = Vector3D();
        position[0] = uniform(0, LENGTH);
        position[1] = uniform(0, WIDTH);
        position[2] = uniform(0, 0.5 * HEIGHT);
        PointMass pt = PointMass(position, false);

        grid[int(position[2]) + HEIGHT * int(position[1]) + WIDTH * HEIGHT * int(position[0])].push_back(pt);
    }
}

std::vector<PointMass> &Fluid::get_position(int x, int y, int z) {
    return grid[z + HEIGHT * y + WIDTH * HEIGHT * x];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, std::vector<Vector3D> external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    /** TODO update position based on velocity **/
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass pt: this->grid[index]) {
            pt.position = pt.position; // + pt.velocity(delta_t);
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

    /** TODO: handle self collisions **/

    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass pt: this->grid[index]) {
            pt.last_position = pt.position;
        }
    }
}