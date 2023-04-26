#include "iostream"
#include <nanogui/nanogui.h>

// #include "../clothMesh.h"
#include "../FluidMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass &pm, double delta_t) {
    // TODO (Part 3): Handle collisions with planes.
    if (dot(pm.position - pm.velocity * delta_t - this->point, this->normal) * dot(pm.position - this->point, this->normal) >= 0) {
        return; // no collision
    }

    Vector3D vec = pm.velocity * delta_t;
    // dot(pm.position - t * vec, normal) = 0
    // dot product is linear
    // dot(pm.position, normal) - t * dot(vec, normal) = 0
    double t = dot(pm.position, normal) / dot(vec, normal);

    // calculating offset
    double currDistance = dot(pm.position - this->point, this->normal);
    double ratio = (currDistance - SURFACE_OFFSET) / currDistance;

    Vector3D correctionVector = t * vec * ratio;

    pm.position = pm.position - correctionVector * (1 - friction);

    pm.velocity = pm.velocity - 2 * this->normal * pm.velocity;
}

void Plane::render(GLShader &shader) {
    nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

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
