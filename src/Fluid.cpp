//
// Created by Wentinn Liao on 4/12/23.
//

#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528);

/** TODO: add fluid density maybe */
Fluid::Fluid(double length, double width, double height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params) {
    double volume = length * width * height;
    double true_num_particles = volume * params.density * 6.022E23 / params.molar_mass;
    double ratio = true_num_particles / nParticles;

    this->SMOOTHING_RADIUS = params.average_distance * std::cbrt(ratio);
    this->PARTICLE_MASS = volume * params.density / nParticles;
    this->SELF_KERNEL = 1 / (PI * std::pow(this->SMOOTHING_RADIUS, 3));
    this->KERNEL_COEFF = 1.5 * this->SELF_KERNEL;

    this->G_LENGTH = (int) (LENGTH / (2 * this->SMOOTHING_RADIUS));
    this->G_WIDTH = (int) (WIDTH / (2 * this->SMOOTHING_RADIUS));
    this->G_HEIGHT = (int) (HEIGHT / (2 * this->SMOOTHING_RADIUS));

    grid = new std::vector<PointMass>[G_LENGTH * G_WIDTH * G_HEIGHT];

    // make planes
    // TODO: change LENGTH and WIDTH to G_LENGTH and G_WIDTH?
    double wall_friction = 0.3;
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(1, 0, 0), wall_friction));
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(0, 1, 0), wall_friction));
    this->collisionObjects.push_back(new Plane(Vector3D(0, 0, 0), Vector3D(0, 0, 1), wall_friction));
    this->collisionObjects.push_back(new Plane(Vector3D(LENGTH, 0, 0), Vector3D(-1, 0, 0), wall_friction));
    this->collisionObjects.push_back(new Plane(Vector3D(0, WIDTH, 0), Vector3D(0, -1, 0), wall_friction));

    // making particles
    // using code from https://stackoverflow.com/questions/38244877/how-to-use-stdnormal-distribution
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device random_device;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(random_device());
    // std of water velocity is approx 642.50364 m/s, average velocity over #ratio particles
    // negligible, may remove random sampling completely
    std::normal_distribution<double> norm_dist_gen(0, params.rms_velocity / ratio);
    for (int i = 0; i < nParticles; i += 1) {
        // TODO: change LENGTH and WIDTH to G_LENGTH and G_WIDTH? (and HEIGHT to G_HEIGHT)?
        // No this is fine, we use LENGTH for position in space, G_LENGTH is number of grid cells
        Vector3D position = Vector3D(uniform(0, LENGTH), uniform(0, WIDTH), uniform(0, 0.5 * HEIGHT));
        Vector3D velocity = Vector3D(norm_dist_gen(gen), norm_dist_gen(gen), norm_dist_gen(gen));
        PointMass pt = PointMass(position, velocity, false);

        this->get_position(position).push_back(pt);
    }
}

std::vector<PointMass> &Fluid::get_position(const Vector3D &pos) const {
    Vector3D indices = pos / (2 * this->SMOOTHING_RADIUS);
    return grid[(int) indices[2] + G_HEIGHT * ((int) indices[1] + G_WIDTH * (int) indices[0])];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    /** TODO: use lambda * h / vmax for timestep, lambda = 0.4, compute running vmax */
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    Vector3D total_external_force(0);
    for (const Vector3D &acc: external_accelerations)
        total_external_force += acc;
    total_external_force *= PARTICLE_MASS;

    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            pt.forces += total_external_force;
        }
    }

    /** intermolecular forces. this (commented-out) code
     * doesn't assume that neighboring boxes of size SMOOTHING_RADIUS do not affect each other **/
    /*
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            for (int i = max((int)ceil(-SMOOTHING_RADIUS) + index / WIDTH * HEIGHT, 0); i <= SMOOTHING_RADIUS; i++) {
                for (int j = max((int)ceil(-SMOOTHING_RADIUS) + (index / HEIGHT) % WIDTH, 0); j <= SMOOTHING_RADIUS; j++) {
                    for (int k = max((int)ceil(-SMOOTHING_RADIUS) + index % HEIGHT,0); k <= SMOOTHING_RADIUS; k++) {

                    }
                }
            }
        }
    }
    */

    /** Take out this part, should be handled in Verlet integration */
    for (int index = 0; index < LENGTH * WIDTH * HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            // pt.velocity += total_acc;
        }
    }

    /** handle collisions with objects **/
    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            for (CollisionObject *co: this->collisionObjects) {
                co->collide(pt);

            }
        }
    }

    /** TODO: handle self collisions */

    /** update position **/
    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index += 1) {
        for (PointMass &pt: this->grid[index]) {
            pt.position = pt.position + pt.velocity * delta_t;
        }
    }
}

/** Kernel function **/
// If use Gaussian, I can definitely kernelize this for O(n) batch computation
double Fluid::W(const PointMass &pi, const PointMass &pj) const {
    double q = (pi.position - pj.position).norm() / this->SMOOTHING_RADIUS;
    double c = 0;
    if (q < 1) {
        c = 2. / 3 + q * q * (q / 2 - 1);
    } else if (q <= 2) {
        c = std::pow(2 - q, 3) / 6;
    }
    return c * this->KERNEL_COEFF;
}

/** gradient of Kernel **/
Vector3D Fluid::grad_W(const PointMass &pi, const PointMass &pj) const {
    Vector3D xji = pj.position - pi.position;
    double q = xji.norm() / this->SMOOTHING_RADIUS;
    double c = 0;
    if (q < 1) {
        c = -2 + 1.5 * q;
    } else if (q <= 2) {
        c = -0.5 * q + 2 - 2 / q;
    }
    return (c * this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS)) * xji;
}

std::vector<double> Fluid::batch_density(int index) const {
    const std::vector<PointMass> &cell = this->grid[index];
    size_t n = cell.size();
    std::vector<double> result(n);

    for (int i = 0; i < n; i++) {
        result[i] += this->SELF_KERNEL / 2;
        for (int j = i + 1; j < n; j++) {
            double w = this->W(cell[i], cell[j]);
            result[i] += w;
            result[j] += w;
        }
        result[i] *= (2 * this->PARTICLE_MASS);
    }
    return result;
}

std::vector<Vector3D> Fluid::batch_grad_pressure(int index) const {
    std::vector<Vector3D> batch;
    for (PointMass &pt: this->grid[index]) {
        Vector3D v = Vector3D(0);
        for (int i = 0; i < )
    }

}

void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
}
