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

#define RELAXATION_EPS 0.00001
#define VORTICITY_EPS 0.0005
#define XSPH_EPS 0.1


int PointMass::ID_NUM = 0;

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528, 1E-6, 2.1510E9, 7);

Fluid::Fluid(double length, double width, double height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params), grid_toggle(true) {
    this->VOLUME = length * width * height;
    double true_num_particles = VOLUME * PARAMS.density * 6.022E23 / PARAMS.molar_mass;
    double ratio = true_num_particles / nParticles;

    this->SMOOTHING_RADIUS = 1.8 * PARAMS.average_distance * std::cbrt(ratio);
    this->PARTICLE_MASS = VOLUME * PARAMS.density / nParticles;

    cout << this->SMOOTHING_RADIUS << endl;

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
    // this->collisionObjects.push_back(new Plane(Vector3D(0, 0, HEIGHT), Vector3D(0, 0, -1), wall_friction));


    // making particles
    // using code from https://stackoverflow.com/questions/38244877/how-to-use-stdnormal-distribution
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device random_device;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(random_device());
    // std of water velocity is approx 642.50364 m/s, average velocity over #ratio particles
    // negligible, may remove random sampling completely
    std::normal_distribution<double> norm_dist(0, params.rms_velocity / ratio);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    for (int i = 0; i < nParticles; i += 1) {
        Vector3D position = Vector3D(LENGTH * (0.5 * (1 - std::sqrt(0.5)) + std::sqrt(0.5) * uniform_dist(gen)), WIDTH * (0.5 * (1 - std::sqrt(0.5)) + std::sqrt(0.5) * uniform_dist(gen)), HEIGHT * (2 * uniform_dist(gen)));
        Vector3D velocity = Vector3D(norm_dist(gen), norm_dist(gen), norm_dist(gen));
        PointMass *pm = new PointMass(position, velocity, false);
        this->list.push_back(pm);
        this->get_cell(position).push_back(pm);
    }
    FluidMesh *fmesh = new FluidMesh();
    this->mesh = fmesh;

    this->global_neighbor_indices = std::vector<std::vector<int>>(G_SIZE);
    for (int index = 0; index < G_SIZE; index++) {
        int l = index / (G_WIDTH * G_HEIGHT), w = (index / G_HEIGHT) % G_WIDTH, h = index % G_HEIGHT;
        this->global_neighbor_indices[index] = std::vector<int>(0);
        for (int tl = max(0, l - 1); tl < min(G_LENGTH, l + 2); tl++) {
            for (int tw = max(0, w - 1); tw < min(G_WIDTH, w + 2); tw++) {
                for (int th = max(0, h - 1); th < min(G_HEIGHT, h + 2); th++) {
                    this->global_neighbor_indices[index].push_back(th + G_HEIGHT * (tw + G_WIDTH * tl));
                }
            }
        }
    }

    double average_density = 0;
    for (int index = 0; index < G_SIZE; index++) {
        std::vector<double> density = this->batch_density(this->batch_W(index));
        average_density = std::accumulate(density.begin(), density.end(), average_density);
    }
    average_density /= NUM_PARTICLES;
    cout << "Average starting density " << average_density << endl;
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

void
Fluid::simulate(double frames_per_sec, double simulation_steps, const std::vector<Vector3D> &external_accelerations) {
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    double start_t = (double) chrono::duration_cast<chrono::nanoseconds>(
            chrono::system_clock::now().time_since_epoch()).count();

    /** Predicted movement */
    this->forward_movement(external_accelerations, delta_t);

    /** Position updates */
    this->incompressibility_adjustment(1, 0.15);

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
            global_W[index] = this->batch_W(index);
            global_scaled_grad_W[index] = this->batch_scaled_grad_W(index);
            global_density[index] = this->batch_density(global_W[index]);
        }
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(velocity_update, thread_num);
    }
    for (auto &thread: threads) {
        thread.join();
    }
#endif
    std::vector<std::vector<Vector3D>> global_curl_velocity = this->global_curl_velocity(global_scaled_grad_W);
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
            global_normalized_curl_velocity_norm, global_scaled_grad_W);

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
    for (auto &thread: threads) {
        thread.join();
    }
#endif

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
    for (auto &thread: threads) {
        thread.join();
    }
#endif

    // this->buildFluidMesh();
    this->debugFluidMesh();

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
    for (auto &thread: threads) {
        thread.join();
    }
#endif
    this->collision_update(delta_t);
    this->cell_update();
}


void Fluid::incompressibility_adjustment(int n_iter, double step_size) {
    /** Position updates */
    double scorr_coeff = GRAD_KERNEL_COEFF / PARAMS.density;
    for (int _ = 0; _ < n_iter; _++) {
        std::vector<std::vector<int>> global_neighbor_indices(G_SIZE);
        std::vector<std::vector<std::vector<double>>> global_W(G_SIZE);
        std::vector<std::vector<std::vector<Vector3D>>> global_scaled_grad_W(G_SIZE);

        std::vector<std::vector<double>> global_density(G_SIZE);
        std::vector<std::vector<std::vector<double>>> global_scorr(G_SIZE);

#ifdef MULTITHREAD
        auto compute_lambda_scorr = [&](int thread_num) {
            for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
            for (int index = 0; index < G_SIZE; index++) {
#endif
                global_W[index] = this->batch_W(index);
                global_scaled_grad_W[index] = this->batch_scaled_grad_W(index);

                global_density[index] = this->batch_density(global_W[index]);
                global_scorr[index] = this->batch_scorr(global_W[index]);
            }
#ifdef MULTITHREAD
        };
        std::thread threads[N_THREADS];
        for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
            threads[thread_num] = std::thread(compute_lambda_scorr, thread_num);
        }
        for (auto &thread: threads) {
            thread.join();
        }
#endif
        std::vector<std::vector<Vector3D>> global_lambda_displacement = this->global_lambda_displacement(
                global_density,
                global_W,
                global_scaled_grad_W
        );

#ifdef MULTITHREAD
        auto position_update = [&](int thread_num) {
            for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
            for (int index = 0; index < G_SIZE; index++) {
#endif
                const std::vector<PointMass *> &cell = this->grid()[index];
                size_t n = cell.size();

                const std::vector<int> &neighbor_indices = global_neighbor_indices[index];
                const std::vector<std::vector<Vector3D>> &scaled_grad_W = global_scaled_grad_W[index];

                for (int i = 0; i < n; i++) {
                    Vector3D scorr_displacement(0);
                    int j = 0;
                    for (int neighbor_index: neighbor_indices) {
                        const std::vector<PointMass *> &neighbors = this->grid()[neighbor_index];
                        for (int dj = 0; dj < this->grid()[neighbor_index].size(); dj++) {
                            scorr_displacement += global_scorr[index][i][j] * scaled_grad_W[i][j];
                        }
                        j += neighbors.size();
                    }
                    cell[i]->tentative_position += step_size * (global_lambda_displacement[index][i] + scorr_coeff * scorr_displacement);
                }
            }
#ifdef MULTITHREAD
        };
        for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
            threads[thread_num] = std::thread(position_update, thread_num);
        }
        for (auto &thread: threads) {
            thread.join();
        }
#endif
        this->constrain_update();
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

void Fluid::constrain_update() {
    /** TODO: update to use tentative position */
    for (int i = 0; i < NUM_PARTICLES; i++) {
        PointMass *pm = this->list[i];
        for (CollisionObject *co: this->collisionObjects) {
            co->constrain(*pm);
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
    this->grid_toggle = !this->grid_toggle;
    // delete[] this->grid;
    // this->grid = new_grid;
    // cout << "cell update end" << endl;
}


void Fluid::buildFluidMesh() {
    /** TODO: implement mesh construction */
    mesh->triangles.clear();

    vector<double> pressures = vector<double>((G_LENGTH + 1) * (G_WIDTH + 1) * (G_HEIGHT + 1));

    for (int i = 0; i <= G_LENGTH; i += 1) {
        for (int j = 0; j <= G_WIDTH; j += 1) {
            for (int k = 0; k <= G_HEIGHT; k += 1) {
                int index = k + (G_WIDTH + 1) * (j + (G_LENGTH + 1) * i);
                pressures[index] = 0.0;
            }
        }
    }

    vector<int> dx = {0, 0, 0, 0, 1, 1, 1, 1};
    vector<int> dy = {0, 0, 1, 1, 0, 0, 1, 1};
    vector<int> dz = {0, 1, 0, 1, 0, 1, 0, 1};

    for (int i = 0; i < G_LENGTH; i += 1) {
        for (int j = 0; j < G_WIDTH; j += 1) {
            for (int k = 0; k < G_HEIGHT; k += 1) {
                int index = k + G_WIDTH * (j + G_LENGTH * i);
                vector<PointMass *> points = this->grid()[index];
                for (int q = 0; q < 8; q += 1) {
                    Vector3D pos = Vector3D(i + dx[q], j + dy[q], k + dz[q]);
                    int index2 = index + dz[q] + G_WIDTH * (dy[q] + G_LENGTH * dx[q]);
                    for (PointMass *pm: points) {
                        pressures[index2] += 1 / min(0.01, (pos - pm->position).norm());
                    }
                }
            }
        }
    }

    // probably don't need this
    vector<int> cubeBitmap = vector<int>(G_SIZE);

    vector<int> neighbors = {4, 2, 1};

    // TODO change value
    double THRESHOLD = 0.1;

    for (int i = 0; i < G_LENGTH; i += 1) {
        for (int j = 0; j < G_WIDTH; j += 1) {
            for (int k = 0; k < G_HEIGHT; k += 1) {
//                cout << i << " " << j << " " << k << endl;
                int index = k + G_WIDTH * (j + G_LENGTH * i);
                cubeBitmap[index] = 0;
                for (int q = 0; q < 8; q += 1) {
                    int index2 = index + dz[q] + G_WIDTH * (dy[q] + G_LENGTH * dx[q]);
                    if (pressures[index2] < THRESHOLD) {
                        cubeBitmap[index] |= (1 << q);
                    }
                }

                if (cubeBitmap[index] == 0 || cubeBitmap[index] == 255) { // all on one side of the surface
                    continue;
                }

                int visited = 0;

//                cout << cubeBitmap[index] << endl;

                for (int q = 0; q < 8; q += 1) {
                    if ((visited & (1 << q)) != 0) { // already visited
                        continue;
                    }
                    int index2 = index + dz[q] + G_WIDTH * (dy[q] + G_LENGTH * dx[q]);
                    if (pressures[index2] > THRESHOLD) {
                        continue;
                    }

                    vector<int> queue = vector<int>(); // stores 010
                    vector<Vector3D> vertices = vector<Vector3D>(); //
                    queue.push_back(q);
                    while (queue.size() > 0) {
                        int ind = queue.back();
                        queue.pop_back();
                        // cout << ind << " " << visited << endl;
                        if ((visited & (1 << ind)) != 0) {
                            continue;
                        }

                        // cout << queue.size() << endl;

                        visited |= (1 << ind);

                        for (int toXor: neighbors) {
                            int nextInd = ind ^ toXor;
//                            if (nextInd == 3) {
//                                cout << ind << " " << toXor << " " << nextInd << endl;
//                            }
                            int index3 = index + (nextInd & 1) +
                                         G_WIDTH * (((nextInd & 2) >> 1) + G_LENGTH * ((nextInd & 4) >> 2));
                            if (pressures[index3] > THRESHOLD) {
                                // vertex on edge
                                // edges.push_back(ind << 3 | nextInd);
                                double pressure1 = pressures[index2];
                                double pressure2 = pressures[index3];
                                double t = pressure2 / (pressure2 - pressure1);
                                Vector3D point = Vector3D(i, j, k);
                                point[0] += t * ((ind & 4) >> 2) + (1 - t) * ((nextInd & 4) >> 2);
                                point[1] += t * ((ind & 2) >> 1) + (1 - t) * ((nextInd & 2) >> 1);
                                point[2] += t * ((ind & 1) >> 0) + (1 - t) * ((nextInd & 1) >> 0);
                                point[0] *= LENGTH / G_LENGTH;
                                point[1] *= WIDTH / G_WIDTH;
                                point[2] *= HEIGHT / G_HEIGHT;
                                vertices.push_back(point);
                            } else {
                                queue.push_back(nextInd);
                            }
                        }
                    }
                    // convert vertices into a mesh

                    if (vertices.size() < 3) {
//                        cout << "ERROR " << vertices.size() << endl;
                        continue;
                    } else {
//                        cout << "SUCCESS" << endl;
                    }

                    for (int ii = 0; ii <= vertices.size() - 3; ii += 1) { // TODO change uv values of triangle
                        // i, i + 1, i + 2
                        Triangle triangle(vertices[0], vertices[ii + 1], vertices[ii + 2], vertices[0],
                                          vertices[ii + 1], vertices[ii + 2]);
                        mesh->triangles.push_back(triangle);
                    }

                }
            }
        }
    }
}


void Fluid::debugFluidMesh() {
    mesh->triangles.clear();
    Vector3D offset[4] = {
            {0.,          0.,          0.05},
            {0.,          0.04714045,  -0.01666667},
            {-0.04082483, -0.02357023, -0.01666667},
            {0.04082483,  -0.02357023, -0.01666667}
    };
    for (PointMass *pm: this->list) {
        mesh->triangles.emplace_back(
                pm->position + offset[0],
                pm->position + offset[1],
                pm->position + offset[2],
                pm->position + offset[0],
                pm->position + offset[1],
                pm->position + offset[2]
        );
        mesh->triangles.emplace_back(
                pm->position + offset[0],
                pm->position + offset[2],
                pm->position + offset[3],
                pm->position + offset[0],
                pm->position + offset[2],
                pm->position + offset[3]
        );
        mesh->triangles.emplace_back(
                pm->position + offset[0],
                pm->position + offset[3],
                pm->position + offset[1],
                pm->position + offset[0],
                pm->position + offset[3],
                pm->position + offset[1]
        );
        mesh->triangles.emplace_back(
                pm->position + offset[3],
                pm->position + offset[2],
                pm->position + offset[1],
                pm->position + offset[3],
                pm->position + offset[2],
                pm->position + offset[1]
        );
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
    if (isnan(c)) {
        cout << pi->tentative_position << " " << pj->tentative_position << endl;
        throw std::runtime_error("W produced NaN value");
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


std::vector<std::vector<double>> Fluid::batch_W(int index) const {
    const std::vector<PointMass *> &cell = this->grid()[index];
    size_t n = cell.size();

    std::vector<std::vector<double>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<double>(0);
        for (int neighbor_index: global_neighbor_indices[index]) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->W(cell[i], pm));
            }
        }
    }
    return result;
}


std::vector<std::vector<Vector3D>>
Fluid::batch_scaled_grad_W(int index) const {
    const std::vector<PointMass *> &cell = this->grid()[index];
    size_t n = cell.size();

    std::vector<std::vector<Vector3D>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<Vector3D>(0);
        for (int neighbor_index: global_neighbor_indices[index]) {
            for (PointMass *pm: this->grid()[neighbor_index]) {
                result[i].push_back(this->scaled_grad_W(cell[i], pm));
            }
        }
    }
    return result;
}


std::vector<double> Fluid::batch_density(const std::vector<std::vector<double>> &W) const {
    size_t n = W.size();
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = this->PARTICLE_MASS * std::accumulate(W[i].begin(), W[i].end(), 0.);
        if (isnan(result[i])) {
            throw std::runtime_error("density produced NaN value");
        }
    }
    return result;
}


std::vector<std::vector<Vector3D>> Fluid::global_lambda_displacement(const std::vector<std::vector<double>> &global_density,
                                                                     const std::vector<std::vector<std::vector<double>>> &global_W,
                                                                     const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const {
    std::vector<std::vector<Vector3D>> result(G_SIZE);
#ifdef MULTITHREAD
    auto compute_lambda = [&](int thread_num) {
        for (int index = thread_num; index < G_SIZE; index += N_THREADS) {
#else
        for (int index = 0; index < G_SIZE; index++) {
#endif
            const std::vector<PointMass *> &cell = this->grid()[index];
            size_t n = cell.size();
            result[index] = std::vector<Vector3D>(n);

            const std::vector<int> &neighbor_indices = global_neighbor_indices[index];
            const std::vector<std::vector<double>> &W = global_W[index];
            const std::vector<std::vector<Vector3D>> &scaled_grad_W = global_scaled_grad_W[index];

            for (int i = 0; i < n; i++) {
                double weight = 0;
                Vector3D denominator(0);

                int j = 0;
                for (int neighbor_index: neighbor_indices) {
                    for (int dj = 0; dj < this->grid()[neighbor_index].size(); dj++) {
                        double w = W[i][j + dj];
                        if (global_density[neighbor_index][dj] > PARAMS.density) {
                            denominator += ((w + SELF_KERNEL) * scaled_grad_W[i][j + dj]);
                        }
                        weight += w;
                    }
                    j += this->grid()[neighbor_index].size();
                }
                denominator /= weight;
                result[index][i] = -denominator;
                if (isnan(result[index][i][0])) {
                    throw std::runtime_error("lambda produced NaN value");
                } else if (result[index][i].norm() > 1) {
                    cout << denominator << endl;
                    throw std::runtime_error("lambda produced large value");
                }
            }
        }
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(compute_lambda, thread_num);
    }
    for (auto &thread: threads) {
        thread.join();
    }
#endif
    return result;
}


std::vector<std::vector<double>> Fluid::batch_scorr(const std::vector<std::vector<double>> &W) const {
    size_t n = W.size();
    std::vector<std::vector<double>> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = std::vector<double>(0);
        for (double w: W[i]) {
            result[i].push_back(this->SCORR_COEFF * std::pow(w, 4));
            if (isnan(result[i][result[i].size() - 1])) {
                throw std::runtime_error("scorr produced NaN value");
            }

        }
    }
    return result;
}

std::vector<std::vector<Vector3D>>
Fluid::global_curl_velocity(const std::vector<std::vector<std::vector<Vector3D>>> &global_scaled_grad_W) const {
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
                    const std::vector<PointMass *> &neighbors = this->grid()[neighbor_index];
                    for (int dj = 0; dj < neighbors.size(); dj++) {
                        result[index][i] += cross(neighbors[dj]->velocity - cell[i]->velocity,
                                                  scaled_grad_W[i][j + dj]);
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
    for (auto &thread: threads) {
        thread.join();
    }
#endif
    return result;
}

std::vector<std::vector<Vector3D>>
Fluid::global_normalized_grad_norm_curl_velocity(const std::vector<std::vector<double>> &global_normalized_curl_velocity_norm,
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
                    const std::vector<PointMass *> &neighbors = this->grid()[neighbor_index];
                    for (int dj = 0; dj < neighbors.size(); dj++) {
                        result[index][i] += (normalized_curl_velocity_norm[i] +
                                             global_normalized_curl_velocity_norm[neighbor_index][dj]) *
                                            scaled_grad_W[i][j + dj];
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
    for (auto &thread: threads) {
        thread.join();
    }
#endif
    return result;
}

