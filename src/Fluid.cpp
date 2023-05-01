//
// Created by Wentinn Liao on 4/12/23.
//

#include <numeric>
#include <utility>
#include <chrono>
#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

#define RELAXATION_EPS 0.001
#define VORTICITY_EPS 0.001

int PointMass::ID_NUM = 0;

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528, 1E-6, 2.1510E9, 7);

Fluid::Fluid(double length, double width, double height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params), grid_toggle(true) {
    this->VOLUME = length * width * height;
    double true_num_particles = VOLUME * PARAMS.density * 6.022E23 / PARAMS.molar_mass;
    double ratio = true_num_particles / nParticles;

    this->SMOOTHING_RADIUS = PARAMS.average_distance * std::cbrt(ratio);
    this->PARTICLE_MASS = VOLUME * PARAMS.density / nParticles;

    this->SELF_KERNEL = 1. / (PI * std::pow(this->SMOOTHING_RADIUS, 3));
    this->KERNEL_COEFF = 1.5 * this->SELF_KERNEL;
    this->SCORR_COEFF = -0.001 / std::pow(0.65 * this->KERNEL_COEFF, 4);
    this->SELF_SCORR = this->SCORR_COEFF * std::pow(this->SELF_KERNEL, 4);

    this->CELL_SIZE = 2. * this->SMOOTHING_RADIUS;
    this->G_LENGTH = (int) (LENGTH / CELL_SIZE) + 1;
    this->G_WIDTH = (int) (WIDTH / CELL_SIZE) + 1;
    this->G_HEIGHT = (int) (3 * HEIGHT / CELL_SIZE);
    this->G_SIZE = G_LENGTH * G_WIDTH * G_HEIGHT;

    this->grid1 = new std::vector<PointMass *>[G_SIZE];
    this->grid2 = new std::vector<PointMass *>[G_SIZE];
    cout << "Grid size " << G_SIZE << endl;
    if (this->grid1 == nullptr || this->grid2 == nullptr)
        throw std::bad_alloc();


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
        Vector3D position = Vector3D(LENGTH * std::rand() / RAND_MAX, WIDTH * std::rand() / RAND_MAX,
                                     0.5 * HEIGHT * std::rand() / RAND_MAX);
        Vector3D velocity = Vector3D(norm_dist_gen(gen), norm_dist_gen(gen), norm_dist_gen(gen));
        PointMass *pm = new PointMass(position, velocity, false);
        this->list.push_back(pm);
        this->get_cell(position).push_back(pm);
    }

}

inline std::vector<PointMass *> *Fluid::grid() const {
    return (this->grid_toggle) ? grid1 : grid2;
}

inline int Fluid::get_index(const Vector3D &pos) const {
    Vector3D indices = pos / CELL_SIZE;
    if (pos[0] < 0 || G_LENGTH < pos[0] || pos[1] < 0 || G_WIDTH < pos[1] || pos[2] < 0 || G_HEIGHT < pos[2]) {
        cout << pos << " out of bounds" << endl;
        throw std::runtime_error("particle out of bounds");
    }
    return (int) indices[2] + G_HEIGHT * ((int) indices[1] + G_WIDTH * (int) indices[0]);
}

inline std::vector<PointMass *> &Fluid::get_cell(const Vector3D &pos) const {
    return this->grid()[this->get_index(pos)];
}

std::vector<int> Fluid::neighbor_indices(int index) const {
    const std::vector<PointMass *> &cell = this->grid()[index];

    int l = index / (G_WIDTH * G_HEIGHT), w = (index / G_HEIGHT) % G_WIDTH, h = index % G_HEIGHT;
    std::vector<int> result(0);
    for (int tl = max(0, l - 1); tl < min(G_LENGTH, l + 2); tl++) {
        for (int tw = max(0, w - 1); tw < min(G_WIDTH, w + 2); tw++) {
            for (int th = max(0, h - 1); th < min(G_HEIGHT, h + 2); th++) {
                result.push_back(th + G_HEIGHT * (tw + G_WIDTH * tl));
            }
        }
    }
    return result;
}

void
Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    double start_t = (double) chrono::duration_cast<chrono::nanoseconds>(
            chrono::system_clock::now().time_since_epoch()).count();

    Vector3D total_external_acceleration(0);
    for (const Vector3D &acc: external_accelerations)
        total_external_acceleration += acc;
    // cout << "Acceleration computation done" << endl;

    /** Predicted movement */
    // Vector3D vmax(0);
    for (int index = 0; index < G_SIZE; index++) {
        for (PointMass *pm: this->grid()[index]) {
            pm->velocity += delta_t * total_external_acceleration;
            pm->tentative_position += delta_t * pm->velocity;
            pm->collided = false;

//            if (pm->velocity.norm2() > vmax.norm2())
//                vmax = pm->velocity;
            // cout << "Movement: " << pm->tentative_position << " " << pm->velocity << endl;
//            if (pm->velocity[2] > 10)
//                throw std::runtime_error("Movement prediction resulted in high velocity");
            if (isnan(pm->tentative_position[0]) || isnan(pm->velocity[0]))
                throw std::runtime_error("Movement prediction resulted in NaN position or velocity");
            this->collision_update(pm, delta_t);
        }
    }
    // cout << "Movement vmax " << vmax << endl;
    // cout << "Movement prediction done" << endl;
    this->cell_update();
    // cout << "Movement prediction done" << endl;

    /** Position updates */
    this->incompressibility_adjustment(3, delta_t);
    // cout << "Position updates done" << endl;

    /** Velocity updates */
    std::vector<std::vector<int>> global_neighbor_indices(G_SIZE);
    std::vector<std::vector<std::vector<double>>> global_W(G_SIZE);
    std::vector<std::vector<std::vector<Vector3D>>> global_scaled_grad_W(G_SIZE);
    std::vector<std::vector<double>> global_density(G_SIZE);
    for (int index = 0; index < G_SIZE; index++) {
        for (PointMass *pm: this->grid()[index]) {
            Vector3D prev = pm->velocity;
            if (!pm->collided)
                pm->velocity = (pm->tentative_position - pm->position) / delta_t;
            pm->position = pm->tentative_position;
//            if (pm->velocity[2] > 10) {
//                cout << prev << endl;
//                throw std::runtime_error("Position adjustment/velocity update resulted in high velocity");
//            }
            if (isnan(pm->tentative_position[0]) || isnan(pm->velocity[0]))
                throw std::runtime_error("Velocity update resulted in NaN position or velocity");
            // cout << "Final: " << pm->position << " " << pm->velocity << endl;

            global_neighbor_indices[index] = this->neighbor_indices(index);
            global_W[index] = this->batch_W(index, global_neighbor_indices[index]);
            global_scaled_grad_W[index] = this->batch_scaled_grad_W(index, global_neighbor_indices[index]);
            global_density[index] = this->batch_density(global_W[index]);
        }
    }
    std::vector<std::vector<Vector3D>> global_curl_velocity = this->global_curl_velocity(global_neighbor_indices,
                                                                                         global_scaled_grad_W);
    std::vector<std::vector<double>> global_normalized_curl_velocity_norm(G_SIZE);
    for (int index = 0; index < G_SIZE; index++) {
        size_t n = this->grid()[index].size();
        global_normalized_curl_velocity_norm[index] = std::vector<double>(n);
        for (int i = 0; i < n; i++) {
            global_normalized_curl_velocity_norm[index][i] =
                    global_curl_velocity[index][i].norm() / (global_density[index][i] * global_density[index][i]);
        }
    }
    std::vector<std::vector<Vector3D>> global_normalized_grad_norm_curl_velocity = this->global_normalized_grad_norm_curl_velocity(
            global_neighbor_indices, global_normalized_curl_velocity_norm, global_scaled_grad_W);

    double max_vorticity_acc = 0;
    for (int index = 0; index < G_SIZE; index++) {
        const std::vector<PointMass *> &cell = this->grid()[index];
        size_t n = cell.size();
        for (int i = 0; i < n; i++) {
            Vector3D vorticity_acc = VORTICITY_EPS * cross(global_normalized_grad_norm_curl_velocity[index][i],
                                                           global_curl_velocity[index][i]);
            max_vorticity_acc = max(max_vorticity_acc, vorticity_acc.norm());
            // cout << vorticity_acc << endl;
            cell[i]->velocity += delta_t * vorticity_acc;
//            if (cell[i]->velocity[2] > 10)
//                throw std::runtime_error("Vorticity resulted in high velocity");
        }
    }
    // cout << max_vorticity_acc << endl;

    for (int index = 0; index < G_SIZE; index++) {
        const std::vector<PointMass *> &cell = this->grid()[index];
        size_t n = cell.size();

        const std::vector<int> &neighbor_indices = global_neighbor_indices[index];

        for (int i = 0; i < n; i++) {
            Vector3D dv(0);
            int j = 0;
            for (int neighbor_index: neighbor_indices) {
                std::vector<PointMass *> neighbors = this->grid()[neighbor_index];
                for (int dj = 0; dj < neighbors.size(); dj++) {
                    dv += (neighbors[dj]->velocity - cell[i]->velocity) * global_W[index][i][j + dj];
                }
                j += neighbors.size();
            }
            // cout << (0.1 * VOLUME / NUM_PARTICLES) * dv << endl;
            cell[i]->velocity += (0.1 * VOLUME / NUM_PARTICLES) * dv;
//            if (cell[i]->velocity[2] > 10)
//                throw std::runtime_error("Viscosity resulted in high velocity");
            // cout << "Final: " << cell[i]->position << " " << cell[i]->velocity << endl;
        }
    }

    // cout << "Velocity update vmax " << vmax << endl;
    // cout << "Velocity updates done" << endl;

    double end_t = (double) chrono::duration_cast<chrono::nanoseconds>(
            chrono::system_clock::now().time_since_epoch()).count();
    cout << "Cycle done in " << (end_t - start_t) * 1E-6 << "ms" << endl;
}


void Fluid::incompressibility_adjustment(int n_iter, double delta_t) {
    /** Position updates */
    double coeff = 0.2 * this->KERNEL_COEFF / (PARAMS.density * this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);
    for (int _ = 0; _ < n_iter; _++) {
        for (int index = 0; index < G_SIZE; index++) {
            // cout << index << " of " << G_SIZE << endl;
            const std::vector<PointMass *> &cell = this->grid()[index];
            if (cell.size() > 1) {
                size_t n = cell.size();

                std::vector<int> neighbor_indices = this->neighbor_indices(index);

                std::vector<std::vector<double>> W = this->batch_W(index, neighbor_indices);
                std::vector<std::vector<Vector3D>> scaled_grad_W = this->batch_scaled_grad_W(index, neighbor_indices);
                std::vector<double> density = this->batch_density(W);

                std::vector<double> lambda = this->batch_lambda(density, scaled_grad_W);
                std::vector<std::vector<double>> scorr = this->batch_scorr(W);

                for (int i = 0; i < n; i++) {
                    PointMass *pm = cell[i];
                    Vector3D dp(0);
                    for (int j = 0; j < n; j++) {
                        dp += (this->PARTICLE_MASS * (lambda[i] + lambda[j]) + scorr[i][j]) * scaled_grad_W[i][j];
                    }
                    pm->tentative_position += (coeff * dp);
                    if (isnan(pm->tentative_position[0]) || isnan(pm->velocity[0])) {
                        throw std::runtime_error("Position adjustment resulted in NaN position or velocity");
                    }
//                    cout << "Adjustment: " << pm->tentative_position << " " << pm->velocity << endl;
                    this->collision_update(pm, delta_t);
                }
            }
        }
//        cout << "Max adjustment " << coeff * max_position_change << endl;
//        cout << this->PARTICLE_MASS << endl;
//        cout << "Adjustment vmax " << vmax << endl;
        this->cell_update();
    }
}


void Fluid::collision_update(PointMass *pm, double delta_t) {
    /** TODO: update to use tentative position */
    for (CollisionObject *co: this->collisionObjects) {
        co->collide(*pm, delta_t);
    }
}


void Fluid::cell_update() {
    // std::vector<PointMass *> *new_grid = new std::vector<PointMass *>[G_SIZE];
    // std::vector<PointMass *> new_grid[G_SIZE];
    // cout << "cell update start" << endl;
    std::vector<PointMass *> *old_grid, *new_grid;
    if (this->grid_toggle) {
        old_grid = this->grid1;
        new_grid = this->grid2;
    } else {
        old_grid = this->grid2;
        new_grid = this->grid1;
    }

    // cout << "clearing" << endl;
    for (int index = 0; index < G_SIZE; index++) {
        new_grid[index].clear();
    }
    // cout << "assignment" << endl;
    for (int index = 0; index < G_SIZE; index++) {
        for (PointMass *pm: old_grid[index]) {
            new_grid[this->get_index(pm->tentative_position)].push_back(pm);
        }
    }
    this->grid_toggle = ~this->grid_toggle;
    // delete[] this->grid;
    // this->grid = new_grid;
    // cout << "cell update end" << endl;
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
Vector3D Fluid::scaled_grad_W(PointMass *pi, PointMass *pj) const {
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
        if (isnan((c * xij)[0])) {
            cout << pi->tentative_position << " " << pj->tentative_position << endl;
            throw std::runtime_error("scaled_grad_W produced NaN value");
        }
        return c * xij;
    }
}


std::vector<std::vector<double>> Fluid::batch_W(int index, const std::vector<int> &neighbor_indices) const {
    const std::vector<PointMass *> &cell = this->grid()[index];
    size_t n = cell.size();

    int l = index / (G_WIDTH * G_HEIGHT), w = (index / G_HEIGHT) % G_WIDTH, h = index % G_HEIGHT;
    std::vector<std::vector<double>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<double>(0);
        for (int neighbor_index: neighbor_indices) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->W(cell[i], pm));
            }
        }
    }
//    for (int i = 0; i < n; i++) {
//        result[i] = std::vector<double>(n);
//        for (int j = 0; j < i; j++) {
//            result[i][j] = result[j][i];
//        }
//        for (int j = i + 1; j < n; j++) {
//            result[i][j] = this->W(cell[i], cell[j]);
//        }
//        result[i][i] = this->SELF_KERNEL;
//    }
    return result;
}


std::vector<std::vector<Vector3D>>
Fluid::batch_scaled_grad_W(int index, const std::vector<int> &neighbor_indices) const {
    const std::vector<PointMass *> &cell = this->grid()[index];
    size_t n = cell.size();

    int l = index / (G_WIDTH * G_HEIGHT), w = (index / G_HEIGHT) % G_WIDTH, h = index % G_HEIGHT;
    std::vector<std::vector<Vector3D>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<Vector3D>(0);
        for (int neighbor_index: neighbor_indices) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->scaled_grad_W(cell[i], pm));
            }
        }
    }
//    for (int i = 0; i < n; i++) {
//        result[i] = std::vector<Vector3D>(n);
//        for (int j = 0; j < i; j++) {
//            result[i][j] = -result[j][i];
//        }
//        for (int j = i + 1; j < n; j++) {
//            result[i][j] = this->scaled_grad_W(cell[i], cell[j]);
//        }
//        result[i][i] = {0};
//    }
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


std::vector<double> Fluid::batch_lambda(const std::vector<double> &density,
                                        const std::vector<std::vector<Vector3D>> &scaled_grad_W) const {
    size_t n = density.size();
    double coeff = (PARAMS.density * this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS) /
                   (this->PARTICLE_MASS * this->KERNEL_COEFF);
    coeff *= coeff;

    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        double denominator = 0;
        Vector3D sum(0);
        for (const Vector3D &grad_W: scaled_grad_W[i]) {
            denominator += grad_W.norm2();
            sum += grad_W;
        }
        // cout << "Result " << i << " before large term: " << result[i] << endl;
        // cout << scaled_grad_W[0][1] << endl;
        denominator += sum.norm2();
        result[i] = (coeff * (1 - density[i] / PARAMS.density)) / (denominator + RELAXATION_EPS);
        if (isnan(result[i])) {
//            for (const Vector3D &grad_W: scaled_grad_W[i]) {
//                cout << grad_W << endl;
//            }
            cout << denominator << endl;
            cout << (coeff * (1 - density[i] / PARAMS.density)) << endl;
            throw std::runtime_error("lambda produced NaN value");
        }
    }
    return result;
}


std::vector<std::vector<double>> Fluid::batch_scorr(const std::vector<std::vector<double>> &W) const {
    size_t n = W.size();
    std::vector<std::vector<double>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<double>(0);
        for (double w: W[i]) {
            result[i].push_back(this->SCORR_COEFF * std::pow(w, 4));
//            if (isnan(result[i][j]))
//                throw std::runtime_error("scorr produced NaN value");
        }
//        result[i][i] = this->SELF_SCORR;
    }
    return result;
}

std::vector<std::vector<Vector3D>>
Fluid::global_curl_velocity(const std::vector<std::vector<int>> &global_neighbor_indices,
                            const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const {
    double coeff = (this->KERNEL_COEFF * this->VOLUME) /
                   (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS * this->NUM_PARTICLES);
    // double coeff = this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);

    std::vector<std::vector<Vector3D>> result(G_SIZE);
    for (int index = 0; index < G_SIZE; index++) {
        const std::vector<PointMass *> &cell = this->grid()[index];
        size_t n = cell.size();
        result[index] = std::vector<Vector3D>(n);

        const std::vector<int> &neighbor_indices = global_neighbor_indices[index];
        const std::vector<std::vector<Vector3D>> &scaled_grad_W = global_scaled_grad_W[index];

        for (int i = 0; i < n; i++) {
            int j = 0;
            for (int neighbor_index: neighbor_indices) {
                std::vector<PointMass *> neighbors = this->grid()[neighbor_index];
                for (int dj = 0; dj < neighbors.size(); dj++) {
                    result[index][i] += cross(neighbors[dj]->velocity - cell[i]->velocity, scaled_grad_W[i][j + dj]);
                }
                j += neighbors.size();
            }
            result[index][i] *= coeff;
            // cout << result[index][i] << endl;
        }
    }
    return result;
}

std::vector<std::vector<Vector3D>>
Fluid::global_normalized_grad_norm_curl_velocity(const std::vector<std::vector<int>> &global_neighbor_indices,
                                                 const std::vector<std::vector<double>> &global_normalized_curl_velocity_norm,
                                                 const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const {
    // double coeff = this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);

    std::vector<std::vector<Vector3D>> result(G_SIZE);
    for (int index = 0; index < G_SIZE; index++) {
        const std::vector<PointMass *> &cell = this->grid()[index];
        size_t n = cell.size();
        result[index] = std::vector<Vector3D>(n);

        const std::vector<int> &neighbor_indices = global_neighbor_indices[index];
        const std::vector<double> &normalized_curl_velocity_norm = global_normalized_curl_velocity_norm[index];
        const std::vector<std::vector<Vector3D>> &scaled_grad_W = global_scaled_grad_W[index];

        for (int i = 0; i < n; i++) {
            int j = 0;
            for (int neighbor_index: neighbor_indices) {
                std::vector<PointMass *> neighbors = this->grid()[neighbor_index];
                for (int dj = 0; dj < neighbors.size(); dj++) {
                    result[index][i] += (normalized_curl_velocity_norm[i] +
                                         global_normalized_curl_velocity_norm[neighbor_index][dj]) *
                                        scaled_grad_W[i][j + dj];
                    // cout << (cell[i]->tentative_position - this->grid()[neighbor_index][dj]->tentative_position).norm() / this->SMOOTHING_RADIUS << endl;
                }
                j += neighbors.size();
            }
            if (!(result[index][i] == 0))
                result[index][i].normalize();
            if (isnan(result[index][i][0])) {
                cout << i << endl;
                throw std::runtime_error("normalized grad norm curl velocity produced NaN value");
            }
        }
    }
    return result;
}

