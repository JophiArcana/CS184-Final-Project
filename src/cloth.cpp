
#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

#define uniform(a, b)   (a + (double) (b - a) * std::rand() / RAND_MAX)

using namespace std;


Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
    this->width = width;
    this->height = height;
    this->num_width_points = num_width_points;
    this->num_height_points = num_height_points;
    this->thickness = thickness;

    buildGrid();
    buildClothMesh();
}

Cloth::~Cloth() {
    point_masses.clear();
    springs.clear();

    if (clothMesh) {
        delete clothMesh;
    }
}

void Cloth::buildGrid() {
    // TODO (Part 1): Build a grid of masses and springs.
    double x_inc = (double) this->width / (this->num_width_points - 1);
    double y_inc = (double) this->height / (this->num_height_points - 1);

    for (int yi = 0; yi < this->num_height_points; yi++)
        for (int xi = 0; xi < this->num_width_points; xi++)
            if (this->orientation == HORIZONTAL)
                this->point_masses.emplace_back(Vector3D(xi * x_inc, 1, yi * y_inc), false);
            else
                this->point_masses.emplace_back(Vector3D(xi * x_inc, yi * y_inc, uniform(-0.001, 0.001)), false);

    for (vector<int> &p: this->pinned)
        this->point_masses[p[1] * this->num_width_points + p[0]].pinned = true;

    for (int yi = 0; yi < this->num_height_points; yi++) {
        for (int xi = 0; xi < this->num_width_points; xi++) {
            int index = yi * this->num_width_points + xi;
            if (xi < this->num_width_points - 1)
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + 1],
                        STRUCTURAL
                );
            if (yi < this->num_height_points - 1) {
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + this->num_width_points],
                        STRUCTURAL
                );
            }
            if (xi < this->num_width_points - 1 && yi < this->num_height_points - 1)
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + this->num_width_points + 1],
                        SHEARING
                );
            if (xi > 0 && yi < this->num_height_points - 1)
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + this->num_width_points - 1],
                        SHEARING
                );
            if (xi < this->num_width_points - 2)
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + 2],
                        BENDING
                );
            if (yi < this->num_height_points - 2)
                this->springs.emplace_back(
                        &this->point_masses[index],
                        &this->point_masses[index + 2 * this->num_width_points],
                        BENDING
                );
        }
    }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
    double mass = width * height * cp->density / num_width_points / num_height_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    // TODO (Part 2): Compute total force acting on each point mass.
    Vector3D total_external_force(0);
    for (Vector3D &acc: external_accelerations)
        total_external_force += acc;
    total_external_force *= mass;

    for (PointMass &p: this->point_masses)
        p.forces = total_external_force;    // External force

    for (Spring &s: this->springs) {
        Vector3D ab = s.pm_b->position - s.pm_a->position;
        double ab_norm = ab.norm();

        double spring_force_magnitude = 0;
        if ((s.spring_type == CGL::STRUCTURAL && cp->enable_structural_constraints) ||
            (s.spring_type == CGL::SHEARING && cp->enable_shearing_constraints) ||
            (s.spring_type == CGL::BENDING && cp->enable_bending_constraints)) {
            spring_force_magnitude = cp->ks * (ab_norm - s.rest_length);
        }
        if (s.spring_type == CGL::BENDING)
            spring_force_magnitude *= 0.2;

        Vector3D F = ab * (spring_force_magnitude / ab_norm);
        s.pm_a->forces += F;
        s.pm_b->forces -= F;
    }


    // TODO (Part 2): Use Verlet integration to compute new point mass positions
    for (PointMass &p: this->point_masses) {
        if (!p.pinned) {
            Vector3D new_position = p.position + (1 - 0.01 * cp->damping) * (p.position - p.last_position) +
                                    p.forces * (delta_t * delta_t / mass);
            p.last_position = p.position;
            p.position = new_position;
        }
    }


    // TODO (Part 4): Handle self-collisions.
    this->build_spatial_map();
    for (PointMass &pm: this->point_masses) {
        this->self_collide(pm, simulation_steps);
    }


    // TODO (Part 3): Handle collisions with other primitives.
    for (PointMass &pm: this->point_masses) {
        for (CollisionObject *co: *collision_objects) {
            co->collide(pm);
        }
    }


    // TODO (Part 2): Constrain the changes to be such that the spring does not change
    // in length more than 10% per timestep [Provot 1995].
    for (Spring &s: this->springs) {
        Vector3D last_ab = s.pm_b->last_position - s.pm_a->last_position, ab = s.pm_b->position - s.pm_a->position;
        double last_ab_norm = last_ab.norm(), ab_norm = ab.norm();
        if (ab_norm > 1.1 * last_ab_norm) {
            Vector3D adjustment = (1 - 1.1 * last_ab_norm / ab_norm) * ab;
            if (s.pm_a->pinned) {
                s.pm_b->position -= adjustment;
            } else if (s.pm_b->pinned) {
                s.pm_a->position += adjustment;
            } else {
                s.pm_a->position += adjustment / 2;
                s.pm_b->position -= adjustment / 2;
            }
        }
    }

}

void Cloth::build_spatial_map() {
    for (const auto &entry: map) {
        delete (entry.second);
    }
    map.clear();

    // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (PointMass &pm: this->point_masses) {
        float h = this->hash_position(pm.position);
        if (this->map.find(h) == this->map.end()) {
            this->map[h] = new vector<PointMass *>();
        }
        this->map[h]->push_back(&pm);
    }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
    Vector3D correction(0);
    int num_corrections = 0;
    for (const PointMass *p: *this->map[this->hash_position(pm.position)]) {
        if (&pm != p) {
            Vector3D displacement = pm.position - p->position;
            double d = displacement.norm();
            if (d < 2 * this->thickness) {
                correction += displacement * ((2 * this->thickness) / d - 1);
                num_corrections++;
            }
        }
    }
    if (num_corrections > 0)
        pm.position += correction / (num_corrections * simulation_steps);
}

float Cloth::hash_position(Vector3D pos) {
    // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    double w = 3. * this->width / this->num_width_points,
            h = 3. * this->height / this->num_height_points,
            t = max(w, h);
    double w_trunc = pos[0] - std::fmod(pos[0], w),
            h_trunc = pos[1] - std::fmod(pos[1], h),
            t_trunc = pos[2] - std::fmod(pos[2], t);
    double p = 1e9 + 7;
    return (float) (p * (p * w_trunc + h_trunc) + t_trunc);
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
    PointMass *pm = &point_masses[0];
    for (int i = 0; i < point_masses.size(); i++) {
        pm->position = pm->start_position;
        pm->last_position = pm->start_position;
        pm++;
    }
}

void Cloth::buildClothMesh() {
    if (point_masses.size() == 0) return;

    ClothMesh *clothMesh = new ClothMesh();
    vector<Triangle *> triangles;

    // Create vector of triangles
    for (int y = 0; y < num_height_points - 1; y++) {
        for (int x = 0; x < num_width_points - 1; x++) {
            PointMass *pm = &point_masses[y * num_width_points + x];
            // Get neighboring point masses:
            /*                      *
             * pm_A -------- pm_B   *
             *             /        *
             *  |         /   |     *
             *  |        /    |     *
             *  |       /     |     *
             *  |      /      |     *
             *  |     /       |     *
             *  |    /        |     *
             *      /               *
             * pm_C -------- pm_D   *
             *                      *
             */

            float u_min = x;
            u_min /= num_width_points - 1;
            float u_max = x + 1;
            u_max /= num_width_points - 1;
            float v_min = y;
            v_min /= num_height_points - 1;
            float v_max = y + 1;
            v_max /= num_height_points - 1;

            PointMass *pm_A = pm;
            PointMass *pm_B = pm + 1;
            PointMass *pm_C = pm + num_width_points;
            PointMass *pm_D = pm + num_width_points + 1;

            Vector3D uv_A = Vector3D(u_min, v_min, 0);
            Vector3D uv_B = Vector3D(u_max, v_min, 0);
            Vector3D uv_C = Vector3D(u_min, v_max, 0);
            Vector3D uv_D = Vector3D(u_max, v_max, 0);


            // Both triangles defined by vertices in counter-clockwise orientation
            triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                                             uv_A, uv_C, uv_B));
            triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                                             uv_B, uv_C, uv_D));
        }
    }

    // For each triangle in row-order, create 3 edges and 3 internal halfedges
    for (int i = 0; i < triangles.size(); i++) {
        Triangle *t = triangles[i];

        // Allocate new halfedges on heap
        Halfedge *h1 = new Halfedge();
        Halfedge *h2 = new Halfedge();
        Halfedge *h3 = new Halfedge();

        // Allocate new edges on heap
        Edge *e1 = new Edge();
        Edge *e2 = new Edge();
        Edge *e3 = new Edge();

        // Assign a halfedge pointer to the triangle
        t->halfedge = h1;

        // Assign halfedge pointers to point masses
        t->pm1->halfedge = h1;
        t->pm2->halfedge = h2;
        t->pm3->halfedge = h3;

        // Update all halfedge pointers
        h1->edge = e1;
        h1->next = h2;
        h1->pm = t->pm1;
        h1->triangle = t;

        h2->edge = e2;
        h2->next = h3;
        h2->pm = t->pm2;
        h2->triangle = t;

        h3->edge = e3;
        h3->next = h1;
        h3->pm = t->pm3;
        h3->triangle = t;
    }

    // Go back through the cloth mesh and link triangles together using halfedge
    // twin pointers

    // Convenient variables for math
    int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points - 1) * 2;

    bool topLeft = true;
    for (int i = 0; i < triangles.size(); i++) {
        Triangle *t = triangles[i];

        if (topLeft) {
            // Get left triangle, if it exists
            if (i % num_width_tris != 0) { // Not a left-most triangle
                Triangle *temp = triangles[i - 1];
                t->pm1->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm1->halfedge->twin = nullptr;
            }

            // Get triangle above, if it exists
            if (i >= num_width_tris) { // Not a top-most triangle
                Triangle *temp = triangles[i - num_width_tris + 1];
                t->pm3->halfedge->twin = temp->pm2->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle to bottom right; guaranteed to exist
            Triangle *temp = triangles[i + 1];
            t->pm2->halfedge->twin = temp->pm1->halfedge;
        } else {
            // Get right triangle, if it exists
            if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
                Triangle *temp = triangles[i + 1];
                t->pm3->halfedge->twin = temp->pm1->halfedge;
            } else {
                t->pm3->halfedge->twin = nullptr;
            }

            // Get triangle below, if it exists
            if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
                Triangle *temp = triangles[i + num_width_tris - 1];
                t->pm2->halfedge->twin = temp->pm3->halfedge;
            } else {
                t->pm2->halfedge->twin = nullptr;
            }

            // Get triangle to top left; guaranteed to exist
            Triangle *temp = triangles[i - 1];
            t->pm1->halfedge->twin = temp->pm2->halfedge;
        }

        topLeft = !topLeft;
    }

    clothMesh->triangles = triangles;
    this->clothMesh = clothMesh;
}

