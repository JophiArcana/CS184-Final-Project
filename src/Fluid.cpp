//
// Created by Wentinn Liao on 4/12/23.
//

#include <numeric>
#include <utility>
#include <chrono>
#include <thread>
#include "Fluid.h"
#include "./collision/collisionObject.h"
#include "./collision/plane.h"

#ifdef MULTITHREAD
    #define N_THREADS 8
#endif

#define RELAXATION_EPS 0.001
#define VORTICITY_EPS 0.001
#define XSPH_EPS 0.1


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
    this->GRAD_KERNEL_COEFF = this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);
    this->SCORR_COEFF = -0.001 / std::pow(0.65 * this->KERNEL_COEFF, 4);
    this->SELF_SCORR = this->SCORR_COEFF * std::pow(this->SELF_KERNEL, 4);

    this->CELL_SIZE = 2. * this->SMOOTHING_RADIUS;
    this->G_LENGTH = (int) (LENGTH / CELL_SIZE) + 1;
    this->G_WIDTH = (int) (WIDTH / CELL_SIZE) + 1;
    this->G_HEIGHT = (int) (3 * HEIGHT / CELL_SIZE);
    this->G_SIZE = G_LENGTH * G_WIDTH * G_HEIGHT;

    this->grid1 = new std::vector<PointMass *>[G_SIZE];
    this->grid2 = new std::vector<PointMass *>[G_SIZE];
    this->list = std::vector<PointMass *>();
    this->list.reserve(NUM_PARTICLES);
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

    double start_t = (double) chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch()).count();

    /** Predicted movement */
    this->forward_movement(external_accelerations, delta_t);
    // cout << "Movement prediction done" << endl;

    /** Position updates */
    this->incompressibility_adjustment(3, delta_t);
    // cout << "Position updates done" << endl;

    /** Velocity updates */
    std::vector<std::vector<int>> global_neighbor_indices(G_SIZE);
    std::vector<std::vector<std::vector<double>>> global_W(G_SIZE);
    std::vector<std::vector<std::vector<Vector3D>>> global_scaled_grad_W(G_SIZE);
    std::vector<std::vector<double>> global_density(G_SIZE);
#ifdef MULTITHREAD
    auto velocity_update = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
            for (PointMass *pm: this->grid()[index]) {
                if (!pm->collided)
                    pm->velocity = (pm->tentative_position - pm->position) / delta_t;
                pm->position = pm->tentative_position;
            }
            global_neighbor_indices[index] = this->neighbor_indices(index);
            global_W[index] = this->batch_W(index, global_neighbor_indices[index]);
            global_scaled_grad_W[index] = this->batch_scaled_grad_W(index, global_neighbor_indices[index]);
            global_density[index] = this->batch_density(global_W[index]);
        }
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(velocity_update, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif
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

#ifdef MULTITHREAD
    auto vorticity_update = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
            const std::vector<PointMass *> &cell = this->grid()[index];
            size_t n = cell.size();
            for (int i = 0; i < n; i++) {
                Vector3D vorticity_acc = VORTICITY_EPS * cross(global_normalized_grad_norm_curl_velocity[index][i],
                                                               global_curl_velocity[index][i]);
                cell[i]->velocity += delta_t * vorticity_acc;
            }
        }
#ifdef MULTITHREAD
    };
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(vorticity_update, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif

    // cout << max_vorticity_acc << endl;

    double coeff = XSPH_EPS * VOLUME / NUM_PARTICLES;
#ifdef MULTITHREAD
    auto viscosity_update = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
            const std::vector<PointMass *> &cell = this->grid()[index];
            size_t n = cell.size();

            for (int i = 0; i < n; i++) {
                Vector3D dv(0);
                int j = 0;
                for (int neighbor_index: global_neighbor_indices[index]) {
                    std::vector<PointMass *> neighbors = this->grid()[neighbor_index];
                    for (int dj = 0; dj < neighbors.size(); dj++) {
                        dv += (neighbors[dj]->velocity - cell[i]->velocity) * global_W[index][i][j + dj];
                    }
                    j += neighbors.size();
                }
                cell[i]->velocity += coeff * dv;
            }
        }
#ifdef MULTITHREAD
    };
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(viscosity_update, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif


    // cout << "Velocity update vmax " << vmax << endl;
    // cout << "Velocity updates done" << endl;

    double end_t = (double) chrono::duration_cast<chrono::nanoseconds>(
            chrono::system_clock::now().time_since_epoch()).count();
    cout << "Cycle done in " << (end_t - start_t) * 1E-6 << "ms" << endl;
}


void Fluid::forward_movement(const std::vector<Vector3D> &external_accelerations, double delta_t) {
    Vector3D total_external_acceleration(0);
    for (const Vector3D &acc: external_accelerations)
        total_external_acceleration += acc;

#ifdef MULTITHREAD
    auto movement_update = [&](int thread_num) {
        for (int i = thread_num; i < NUM_PARTICLES; i += N_THREADS) {
#else
        for (int i = 0; i < NUM_PARTICLES; i++) {
#endif
            PointMass *pm = this->list[i];
            pm->velocity += delta_t * total_external_acceleration;
            pm->tentative_position += delta_t * pm->velocity;
            pm->collided = false;
        }
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(movement_update, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif
    this->collision_update(delta_t);
    this->cell_update();
}


void Fluid::incompressibility_adjustment(int n_iter, double delta_t) {
    /** Position updates */
    double coeff = 0.2 * GRAD_KERNEL_COEFF / PARAMS.density;
    for (int _ = 0; _ < n_iter; _++) {
#ifdef MULTITHREAD
        auto position_update = [&](int thread_num) {
            for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
            for (int index = 0; index < G_SIZE; index++) {
#endif
                const std::vector<PointMass *> &cell = this->grid()[index];
                if (cell.size() > 1) {
                    size_t n = cell.size();

                    std::vector<int> neighbor_indices = this->neighbor_indices(index);

                    std::vector<std::vector<double>> W = this->batch_W(index, neighbor_indices);
                    std::vector<std::vector<Vector3D>> scaled_grad_W = this->batch_scaled_grad_W(index,
                                                                                                 neighbor_indices);
                    std::vector<double> density = this->batch_density(W);

                    std::vector<double> lambda = this->batch_lambda(density, scaled_grad_W);
                    std::vector<std::vector<double>> scorr = this->batch_scorr(W);

                    for (int i = 0; i < n; i++) {
                        PointMass *pm = cell[i];
                        Vector3D dp(0);
                        for (int j = 0; j < n; j++) {
                            dp += (PARTICLE_MASS * (lambda[i] + lambda[j]) + scorr[i][j]) * scaled_grad_W[i][j];
                        }
                        pm->tentative_position += (coeff * dp);
                    }
                }
            }
#ifdef MULTITHREAD
        };
        std::thread threads[N_THREADS];
        for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
            threads[thread_num] = std::thread(position_update, thread_num);
        }
        for (auto &thread : threads) {
            thread.join();
        }
#endif
        this->collision_update(delta_t);
        this->cell_update();
    }
}


void Fluid::collision_update(double delta_t) {
    /** TODO: update to use tentative position */
    for (int i = 0; i < NUM_PARTICLES; i++) {
        PointMass *pm = this->list[i];
        for (CollisionObject *co: this->collisionObjects) {
            co->collide(*pm, delta_t);
        }
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
    vector<double> pressures = vector<double>((G_LENGTH + 1) * (G_WIDTH + 1) * (G_HEIGHT + 1));

    for (int i = 0; i <= G_LENGTH; i += 1) {
        for (int j = 0; j <= G_WIDTH; j += 1) {
            for (int k = 0; k <= G_HEIGHT; k += 1) {
                int index = k + (G_WIDTH + 1) * (j + (G_LENGTH + 1) * i);
                pressures[index] = 0; // TODO implement
            }
        }
    }

    vector<int> dx = {0, 0, 0, 0, 1, 1, 1, 1};
    vector<int> dy = {0, 0, 1, 1, 0, 0, 1, 1};
    vector<int> dz = {0, 1, 0, 1, 0, 1, 0, 1};

    // probably don't need this
    vector<int> cubeBitmap = vector<int>(G_LENGTH * G_WIDTH * G_LENGTH);

    FluidMesh mesh;
    vector<int> neighbors = {100, 010, 001};

    for (int i = 0; i < G_LENGTH; i += 1) {
        for (int j = 0; j < G_WIDTH; j += 1) {
            for (int k = 0; k < G_HEIGHT; k += 1) {
                int index = k + G_WIDTH * (j + G_LENGTH * i);
                for (int q = 0; q < 8; q += 1) {
                    int index2 = index + dz[q] + G_WIDTH * (dy[q] + G_LENGTH * dx[q]);
                    if (pressures[index2] < 0.0) { // TODO implement with threshold
                        cubeBitmap[index] += 1 << q;
                    }
                }

                if (cubeBitmap[index] == 0 || cubeBitmap[index] == 255) { // all on one side of the surface
                    continue;
                }

                // TODO convert into triangular mesh
                int visited = 0;

                for (int q = 0; q < 8; q += 1) {
                    if (visited & (1 << q) != 0) { // already visited
                        continue;
                    }
                    int index2 = index + dz[q] + G_WIDTH * (dy[q] + G_LENGTH * dx[q]);
                    if (pressures[index2] > 0.0) { // TODO implement with threshold
                        continue;
                    }

                    vector<int> queue = vector<int>(); // stores 010
                    vector<Vector3D> vertices = vector<Vector3D>(); //
                    queue.push_back(q);
                    while (queue.size() > 0) {
                        int ind = 1; //queue.pop_back();
                        if (visited & (1 << ind) != 0) {
                            continue;
                        }
                        visited &= 1 << ind;

                        for (int toXor: neighbors) {
                            int nextInd = ind ^ toXor;
                            // buggy code!
                            int index3 = index2 + toXor & 1 + G_WIDTH * (toXor & 2 + G_LENGTH * toXor & 4);
                            if (pressures[index3] > 0.0) { // TODO implement with threshold
                                // vertex on edge
                                // edges.push_back(ind << 3 | nextInd);
                                double pressure1 = pressures[index2];
                                double pressure2 = pressures[index3];
                                double t = pressure2 / (pressure2 - pressure1);
                                Vector3D point = Vector3D(i, j, k);
                                point[0] += t * ((ind & 4) >> 2) + (1 - t) * ((nextInd & 4) >> 2);
                                point[1] += t * ((ind & 2) >> 1) + (1 - t) * ((nextInd & 2) >> 1);
                                point[2] += t * ((ind & 1) >> 0) + (1 - t) * ((nextInd & 1) >> 0);
                                vertices.push_back(point);
                            } else {
                                queue.push_back(nextInd);
                            }
                        }
                    }
                    // convert vertices into a mesh

                    if (vertices.size() < 3) {
                        cout << "ERROR" << endl;
                    }

                    for (int i = 0; i <= vertices.size() - 3; i += 1) { // TODO change uv values of triangle
                        // i, i + 1, i + 2
                        Triangle triangle = Triangle(&vertices[i], &vertices[i + 1], &vertices[i + 2], vertices[i], vertices[i + 1], vertices[i + 2]);
                        mesh.triangles.push_back(&triangle);
                    }

                }
            }
        }
    }


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
        // return (c * GRAD_KERNEL_COEFF) * xij;
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

    std::vector<std::vector<double>> result(n);
// #ifdef MULTITHREAD
//    auto compute_W = [&](int thread_num) {
//        for (int i = thread_num; i < n; i += N_THREADS) {
//            result[i] = std::vector<double>(0);
//            for (int neighbor_index: neighbor_indices) {
//                for (PointMass *pm: this->grid()[neighbor_index]) {
//                    result[i].push_back(this->W(cell[i], pm));
//                }
//            }
//        }
//    };
//    std::thread threads[N_THREADS];
//    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
//        threads[thread_num] = std::thread(compute_W, thread_num);
//    }
//    for (auto &thread : threads) {
//        thread.join();
//    }
// #else
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<double>(0);
        for (int neighbor_index: neighbor_indices) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->W(cell[i], pm));
            }
        }
    }
// #endif
    return result;
}


std::vector<std::vector<Vector3D>>
Fluid::batch_scaled_grad_W(int index, const std::vector<int> &neighbor_indices) const {
    const std::vector<PointMass *> &cell = this->grid()[index];
    size_t n = cell.size();

    std::vector<std::vector<Vector3D>> result(n);
//#ifdef MULTITHREAD
//    auto compute_grad_W = [&](int thread_num) {
//        for (int i = thread_num; i < n; i += N_THREADS) {
//            result[i] = std::vector<Vector3D>(0);
//            for (int neighbor_index: neighbor_indices) {
//                for (PointMass *pm: this->grid()[neighbor_index]) {
//                    result[i].push_back(this->scaled_grad_W(cell[i], pm));
//                }
//            }
//        }
//    };
//    std::thread threads[N_THREADS];
//    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
//        threads[thread_num] = std::thread(compute_grad_W, thread_num);
//    }
//    for (auto &thread : threads) {
//        thread.join();
//    }
//#else
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<Vector3D>(0);
        for (int neighbor_index: neighbor_indices) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->scaled_grad_W(cell[i], pm));
            }
        }
    }
// #endif
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
    double coeff = PARAMS.density / (PARTICLE_MASS * GRAD_KERNEL_COEFF);
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
    double coeff = (GRAD_KERNEL_COEFF * VOLUME) / NUM_PARTICLES;

    std::vector<std::vector<Vector3D>> result(G_SIZE);
#ifdef MULTITHREAD
    auto f = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
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
            }
        }
#ifdef  MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(f, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif
    return result;
}

std::vector<std::vector<Vector3D>>
Fluid::global_normalized_grad_norm_curl_velocity(const std::vector<std::vector<int>> &global_neighbor_indices,
                                                 const std::vector<std::vector<double>> &global_normalized_curl_velocity_norm,
                                                 const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const {
    std::vector<std::vector<Vector3D>> result(G_SIZE);

#ifdef MULTITHREAD
    auto f = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
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
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(f, thread_num);
    }
    for (auto &thread : threads) {
        thread.join();
    }
#endif
    return result;
}

