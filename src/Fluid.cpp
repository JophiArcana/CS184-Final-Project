//
// Created by Wentinn Liao on 4/12/23.
//

#include <numeric>
#include <utility>
#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528, 1E-6, 3.0714285E8);


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

    this->grid = new std::vector<PointMass *>[G_LENGTH * G_WIDTH * G_HEIGHT];

    // make planes
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
        Vector3D position = Vector3D(LENGTH * std::rand() / RAND_MAX, WIDTH * std::rand() / RAND_MAX, 0.5 * HEIGHT * std::rand() / RAND_MAX);
        Vector3D velocity = Vector3D(norm_dist_gen(gen), norm_dist_gen(gen), norm_dist_gen(gen));

        this->get_position(position).push_back(new PointMass(position, velocity, false));
    }
}

std::vector<PointMass *> &Fluid::get_position(const Vector3D &pos) const {
    Vector3D indices = pos / (2 * this->SMOOTHING_RADIUS);
    return grid[(int) indices[2] + G_HEIGHT * ((int) indices[1] + G_WIDTH * (int) indices[0])];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    Vector3D total_external_acceleration(0);
    for (const Vector3D &acc: external_accelerations)
        total_external_acceleration += acc;


    /** Acceleration computation */
    double vmax = 0;
    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index += 1) {
        std::vector<std::vector<double>> W = this->batch_W(index);
        std::vector<std::vector<Vector3D>> unnormalized_grad_W = this->batch_unnormalized_grad_W(index);

        std::vector<double> density = this->batch_density(W);
        std::vector<double> pressure = this->batch_pressure(density);

        std::vector<Vector3D> scaled_grad_pressure = this->batch_scaled_grad_pressure(pressure, density, unnormalized_grad_W);
        std::vector<Vector3D> scaled_laplacian_velocity = this->batch_scaled_laplacian_velocity(index, density, unnormalized_grad_W);

        size_t n = this->grid[index].size();
        for (int i = 0; i < n; i++) {
            PointMass *pt = this->grid[index][i];
            pt->acceleration = -scaled_grad_pressure[i] + scaled_laplacian_velocity[i] + total_external_acceleration;

            vmax = std::max(vmax, pt->velocity.norm());
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

    /** Adaptive step integration **/
    double delta_t = 0.4 * this->SMOOTHING_RADIUS / vmax;
    this->timestamps.push_back((this->timestamps.empty()) ? 0 : (this->timestamps.back() + delta_t));

    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index += 1) {
        for (PointMass *pt: this->grid[index]) {
            pt->velocity += pt->acceleration * delta_t;
            pt->position += pt->velocity * delta_t;
        }
    }

    /** handle collisions with objects **/
    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index += 1) {
        for (PointMass *pt: this->grid[index]) {
            for (CollisionObject *co: this->collisionObjects) {
                co->collide(*pt);

            }
        }
    }
}


void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
}


/** Kernel function **/
double Fluid::W(PointMass *pi, PointMass *pj) const {
    double q = (pi->position - pj->position).norm() / this->SMOOTHING_RADIUS;
    double c = 0;
    if (q < 1) {
        c = 2. / 3 + q * q * (q / 2 - 1);
    } else if (q <= 2) {
        c = std::pow(2 - q, 3) / 6;
    }
    return c * this->KERNEL_COEFF;
}


/** Note: all grad_W are scaled by a constant factor, so this scaling is saved for once all
 * vectors have been accumulated for efficiency. */
Vector3D Fluid::unnormalized_grad_W(PointMass *pi, PointMass *pj) const {
    Vector3D xji = pj->position - pi->position;
    double q = xji.norm() / this->SMOOTHING_RADIUS;
    if (q == 0) {
        return {0};
    } else {
        double c = 0;
        if (q < 1) {
            c = -2 + 1.5 * q;
        } else if (q <= 2) {
            c = -0.5 * q + 2 - 2 / q;
        }
        return c * xji;
        // Below is the normalized version. Constant factor taken out for efficiency
        // return (c * this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS)) * xji;
    }
}


std::vector<std::vector<double>> Fluid::batch_W(int index) const {
    const std::vector<PointMass *> &cell = this->grid[index];

    size_t n = cell.size();
    std::vector<std::vector<double>> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            result[i][j] = result[j][i];
        }
        for (int j = i + 1; j < n; j++) {
            result[i][j] = this->W(cell[i], cell[j]);
        }
        result[i][i] = this->SELF_KERNEL;
    }
    return result;
}


std::vector<std::vector<Vector3D>> Fluid::batch_unnormalized_grad_W(int index) const {
    const std::vector<PointMass *> &cell = this->grid[index];

    size_t n = cell.size();
    std::vector<std::vector<Vector3D>> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            result[i][j] = -result[j][i];
        }
        for (int j = i; j < n; j++) {
            result[i][j] = this->unnormalized_grad_W(cell[i], cell[j]);
        }
        result[i][i] = {0};
    }
    return result;
}


std::vector<double> Fluid::batch_density(const std::vector<std::vector<double>> &W) const {
    size_t n = W.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++)
        result[i] = this->PARTICLE_MASS * std::accumulate(W[i].begin(), W[i].end(), 0.);
    return result;
}


std::vector<double> Fluid::batch_pressure(const std::vector<double> &density) const {
    size_t n = density.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++)
        result[i] = PARAMS.tait_coefficient * (std::pow(density[i] / PARAMS.density, 7) - 1);
    return result;
}


/** Note: gradient computation involves multiplying by density, while the force calculation involves dividing
 * so the operation is ignored for efficiency. */
std::vector<Vector3D> Fluid::batch_scaled_grad_pressure(const std::vector<double> &pressure, const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const {
    size_t n = pressure.size();
    double coeff = (this->KERNEL_COEFF * this->PARTICLE_MASS) / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);

    double normalized_pressure[n];
    for (int i = 0; i < n; i++)
        normalized_pressure[i] = pressure[i] / (density[i] * density[i]);

    std::vector<Vector3D> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            result[i] += (normalized_pressure[i] + normalized_pressure[j]) * unnormalized_grad_W[i][j];
        result[i] *= coeff;
    }
    return result;
}


/** All resulting are scaled by the viscosity when making force calculations, so the scaling merged
 * to the coefficient for efficiency. */
std::vector<Vector3D> Fluid::batch_scaled_laplacian_velocity(int index, const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const {
    const std::vector<PointMass *> &cell = this->grid[index];
    size_t n = cell.size();
    double coeff = (2 * PARAMS.kinematic_viscosity * this->KERNEL_COEFF * this->PARTICLE_MASS) / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);
    double eps = 0.01 * this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS;

    std::vector<Vector3D> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Vector3D xij = cell[i]->position - cell[j]->position;
            Vector3D vij = cell[i]->velocity - cell[j]->velocity;
            result[i] += (dot(xij, unnormalized_grad_W[i][j]) / (density[j] * (xij.norm2() + eps))) * vij;
        }
        result[i] *= coeff;
    }
    return result;
}
