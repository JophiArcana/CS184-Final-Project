#include "iostream"
#include <nanogui/nanogui.h>

// #include "../clothMesh.h"
#include "../FluidMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001
#define DAMPING 0.5

void Plane::collide(PointMass &pm, double delta_t) {
    // TODO (Part 3): Handle collisions with planes.
    double d = dot(pm.tentative_position - this->point, this->normal);
    if (d <= 0) {
        pm.velocity = (pm.tentative_position - pm.position) / delta_t;
        // cout << "Plane normal " << this->normal << endl;
        // cout << "Pre-collision position " << pm.tentative_position << endl;
        double vd = dot(pm.velocity, this->normal);

        double crossover_t = d / vd;
        Vector3D collision_position = pm.tentative_position - crossover_t * pm.velocity;

        // cout << "Pre-collision velocity " << pm.velocity << endl;
        Vector3D old_v = pm.velocity;
        pm.velocity = DAMPING * ((-2 * vd) * this->normal + pm.velocity);
//        if (pm.velocity[2] > 10) {
//            cout << "old velocity: " << old_v << endl;
//            throw std::runtime_error("Collision out velocity exceeded old velocity");
//        }
        // cout << "\tPost-collision velocity " << pm.velocity << endl;
        pm.tentative_position = collision_position + SURFACE_OFFSET * this->normal;
        pm.collided = true;

        // cout << "\tPost-collision position " << pm.tentative_position << endl;

        /** Vector3D vec = pm.velocity * delta_t;
        // dot(pm.position - t * vec, normal) = 0
        // dot product is linear
        // dot(pm.position, normal) - t * dot(vec, normal) = 0
        double t = dot(pm.position, normal) / dot(vec, normal);

        // calculating offset
        double currDistance = dot(pm.position - this->point, this->normal);
        double ratio = (currDistance - SURFACE_OFFSET) / currDistance;

        Vector3D correctionVector = t * vec * ratio;

        pm.position = pm.position - correctionVector * (1 - friction);

        pm.velocity = pm.velocity - 2 * this->normal * pm.velocity; */
    }
}

void Plane::render(GLShader &shader) {
    nanogui::Color color(0.5f, 1.f, 1.f, 0.4f);

    Vector3f sPoint(point.x, point.y, point.z);
    Vector3f sNormal(normal.x, normal.y, normal.z);
    Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                       normal.x - normal.y);
    sParallel.normalize();
    Vector3f sCross = sNormal.cross(sParallel);

    MatrixXf positions(3, 4);
    MatrixXf normals(3, 4);

    positions.col(0) << sPoint + 2 * (sCross + sParallel);
    positions.col(1) << sPoint + 2 * (sCross - sParallel);
    positions.col(2) << sPoint + 2 * (-sCross + sParallel);
    positions.col(3) << sPoint + 2 * (-sCross - sParallel);

    normals.col(0) << sNormal;
    normals.col(1) << sNormal;
    normals.col(2) << sNormal;
    normals.col(3) << sNormal;

    if (shader.uniform("u_color", false) != -1) {
        shader.setUniform("u_color", color);
    }
    shader.uploadAttrib("in_position", positions);
    if (shader.attrib("in_normal", false) != -1) {
        shader.uploadAttrib("in_normal", normals);
    }

    shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
}
