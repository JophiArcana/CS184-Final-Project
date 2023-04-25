#include <nanogui/nanogui.h>

// #include "../clothMesh.h"
#include "../FluidMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
    // TODO (Part 3): Handle collisions with spheres.

    Vector3D vec = pm.position - this->origin;
    double dist = vec.norm();
    if (dist >= this->radius) {
        return;
    }

    Vector3D tangentPointSphere = vec * this->radius / dist;

    Vector3D correctionVec = tangentPointSphere + this->origin - pm.position;
    pm.position = pm.position + (1 - this->friction) * correctionVec;
}

void Sphere::render(GLShader &shader) {
    // We decrease the radius here so flat triangles don't behave strangely
    // and intersect with the sphere when rendered
    m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
