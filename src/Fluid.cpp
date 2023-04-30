//
// Created by Wentinn Liao on 4/12/23.
//

#include <numeric>
#include <utility>
#include <chrono>
#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528, 1E-6, 2.1510E9, 7);

Fluid::Fluid(double length, double width, double height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params) {
    double volume = length * width * height;
    double true_num_particles = volume * PARAMS.density * 6.022E23 / PARAMS.molar_mass;
    double ratio = true_num_particles / nParticles;

    this->SMOOTHING_RADIUS = PARAMS.average_distance * std::cbrt(ratio);
    this->PARTICLE_MASS = volume * PARAMS.density / nParticles;
    this->SELF_KERNEL = 1. / (PI * std::pow(this->SMOOTHING_RADIUS, 3));
    this->KERNEL_COEFF = 1.5 * this->SELF_KERNEL;

    this->G_LENGTH = (int) (LENGTH / (2 * this->SMOOTHING_RADIUS)) + 1;
    this->G_WIDTH = (int) (WIDTH / (2 * this->SMOOTHING_RADIUS)) + 1;
    this->G_HEIGHT = (int) (HEIGHT / (2 * this->SMOOTHING_RADIUS)) + 1;
    this->G_BUFFER = max(max(this->G_LENGTH, this->G_WIDTH), this->G_HEIGHT) / 4;

    cout << "Here " << PARAMS.molar_mass / (6.022E23 * PARAMS.density * std::pow(PARAMS.average_distance, 4)) << endl;

    this->grid = new std::deque<PointMass *>[(G_LENGTH + 2 * G_BUFFER) * (G_WIDTH + 2 * G_BUFFER) * (G_HEIGHT + 2 * G_BUFFER)];


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

std::deque<PointMass *> &Fluid::get_position(const Vector3D &pos) const {
    Vector3D indices = pos / (2 * this->SMOOTHING_RADIUS);
    indices += Vector3D(G_BUFFER);
    return grid[(int) indices[2] + (G_HEIGHT + 2 * G_BUFFER) * ((int) indices[1] + (G_WIDTH + 2 * G_BUFFER) * (int) indices[0])];
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    double start_t = (double) chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count() / 1000;

    Vector3D total_external_acceleration(0);
    for (const Vector3D &acc: external_accelerations)
        total_external_acceleration += acc;

    /** Predicted movement */
    for (int h = G_BUFFER; h < G_HEIGHT + G_BUFFER; h++) {
        for (int w = G_BUFFER; w < G_WIDTH + G_BUFFER; w++) {
            for (int l = G_BUFFER; l < G_LENGTH + G_BUFFER; l++) {
                std::deque<PointMass *> cell = this->grid[h + (G_HEIGHT + 2 * G_BUFFER) * (w + (G_WIDTH + 2 * G_BUFFER) * l)];
                PointMass *pm;
                while ((pm = cell.front())->stage) {
                    cell.pop_front();
                    pm->velocity += delta_t * total_external_acceleration;
                    pm->tentative_position += delta_t * pm->velocity;
                    pm->stage = false;
                    this->get_position(pm->tentative_position).push_back(pm);
                }
            }
        }
    }



    /** Acceleration computation */

    double vmax = 0;
    for (int index = 0; index < G_LENGTH * G_WIDTH * G_HEIGHT; index++) {
        size_t n = this->grid[index].size();
        if (n == 0) {
            continue;
        }
        std::vector<std::vector<double>> W = this->batch_W(index);
        std::vector<std::vector<Vector3D>> unnormalized_grad_W = this->batch_unnormalized_grad_W(index);

        std::vector<double> density = this->batch_density(W);
        std::vector<double> pressure = this->batch_pressure(density);


        std::vector<Vector3D> scaled_grad_pressure = this->batch_scaled_grad_pressure(pressure, density, unnormalized_grad_W);
        std::vector<Vector3D> scaled_laplacian_velocity = this->batch_scaled_laplacian_velocity(index, density, unnormalized_grad_W);


        for (int i = 0; i < n; i++) {
            PointMass *pt = this->grid[index][i];
            pt->acceleration = -scaled_grad_pressure[i] + scaled_laplacian_velocity[i] + total_external_acceleration;
            // cout << pt->acceleration << endl;
            cout << -scaled_grad_pressure[i] << " " << scaled_laplacian_velocity[i] << " " << total_external_acceleration << endl;
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
    delta_t = min(delta_t, 0.001);
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
                co->collide(*pt, delta_t);
            }
        }
    }

    /** std::vector<PointMass *> *ngrid = new std::vector<PointMass *>[G_LENGTH * G_WIDTH * G_HEIGHT];

    for (int i = 0; i < G_LENGTH * G_WIDTH * G_HEIGHT; i += 1) {
        for (PointMass *pm : grid[i]) {
            Vector3D indices = pm->position / (2 * this->SMOOTHING_RADIUS);
            cout << pm->position << " " << indices << endl;
            int index = (int) indices[2] + G_HEIGHT * ((int) indices[1] + G_WIDTH * (int) indices[0]);
            ngrid[index].push_back(pm);
        }
    }

    delete[] grid;
    grid = ngrid; */

    // cout << "max v: " << vmax << endl;
    double end_t = (double) chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count() / 1000;
    // cout << end_t - start_t << endl;
}


void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
}


/** Kernel function **/
double Fluid::W(PointMass *pi, PointMass *pj) const {
    double q = (pi->tentative_position - pj->tentative_position).norm() / this->SMOOTHING_RADIUS;
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
    Vector3D xij = pi->tentative_position - pj->tentative_position;
    double q = xij.norm() / this->SMOOTHING_RADIUS;
    if (q == 0) {
        return {0};
    } else {
        double c = 0;
        if (q < 1) {
            c = -2 + 1.5 * q;
        } else if (q <= 2) {
            c = -0.5 * q + 2 - 2 / q;
        }
        // Below is the normalized version. Constant factor taken out for efficiency
        // return (c * this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS)) * xij;
        return c * xij;
    }
}


std::vector<std::vector<double>> Fluid::batch_W(int index) const {
    const std::deque<PointMass *> &cell = this->grid[index];

    size_t n = cell.size();
    std::vector<std::vector<double>> result(n);
    int i = 0;
    for (auto iter_i = cell.begin(); iter_i != cell.end(); iter_i++) {
        result[i] = std::vector<double>(n);
        for (int j = 0; j < i; j++) {
            result[i][j] = result[j][i];
        }
        int j = i + 1;
        for (auto iter_j = std::next(iter_i); iter_j != cell.end(); iter_j++) {
            result[i][j] = this->W(*iter_i, *iter_j);
            j++;
        }
        result[i][i] = this->SELF_KERNEL;
        i++;
    }
    return result;
}


std::vector<std::vector<Vector3D>> Fluid::batch_unnormalized_grad_W(int index) const {
    const std::deque<PointMass *> &cell = this->grid[index];

    size_t n = cell.size();
    std::vector<std::vector<Vector3D>> result(n);
    int i = 0;
    for (auto iter_i = cell.begin(); iter_i != cell.end(); iter_i++) {
        result[i] = std::vector<Vector3D>(n);
        for (int j = 0; j < i; j++) {
            result[i][j] = -result[j][i];
        }
        int j = i + 1;
        for (auto iter_j = std::next(iter_i); iter_j < cell.end(); iter_j++) {
            result[i][j] = this->unnormalized_grad_W(*iter_i, *iter_j);
            j++;
        }
        result[i][i] = {0};
        i++;
    }
    return result;
}


std::vector<double> Fluid::batch_lambda(const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const {
    size_t n = density.size();
    double coeff = (PARAMS.density * this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS) / (this->PARTICLE_MASS * this->KERNEL_COEFF);
    coeff *= coeff;

    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double norm2 = unnormalized_grad_W[i][j].norm2();
            result[i] += norm2;
            result[j] += norm2;
        }
        result[i] += std::accumulate(unnormalized_grad_W[i].begin(), unnormalized_grad_W[i].end(), Vector3D(0)).norm2();
        result[i] = (coeff * (1 - density[i] / PARAMS.density)) / result[i];
    }
    return result;
}


std::vector<double> Fluid::batch_density(const std::vector<std::vector<double>> &W) const {
    size_t n = W.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = this->PARTICLE_MASS * std::accumulate(W[i].begin(), W[i].end(), 0.);
    }
    return result;
}


std::vector<double> Fluid::batch_pressure(const std::vector<double> &density) const {
    size_t n = density.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = (PARAMS.tait_coefficient / PARAMS.tait_gamma) * (std::pow(density[i] / PARAMS.density, PARAMS.tait_gamma) - 1);
        // cout << result[i] << endl;
    }
    return result;
}


/** Note: gradient computation involves multiplying by density, while the force calculation involves dividing
 * so the operation is ignored for efficiency. */
std::vector<Vector3D> Fluid::batch_scaled_grad_pressure(const std::vector<double> &pressure, const std::vector<double> &density, const std::vector<std::vector<Vector3D>> &unnormalized_grad_W) const {
    size_t n = pressure.size();
    double coeff = (this->KERNEL_COEFF * this->PARTICLE_MASS) / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);
    double normalized_pressure[n];
    for (int i = 0; i < n; i++) {
        normalized_pressure[i] = pressure[i] / (density[i] * density[i]);
    }

    std::vector<Vector3D> result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += (normalized_pressure[i] + normalized_pressure[j]) * unnormalized_grad_W[i][j];
        }
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
