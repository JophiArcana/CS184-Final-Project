//
// Created by Wentinn Liao on 4/12/23.
//

#ifndef CLOTHSIM_FLUIDMESH_H
#define CLOTHSIM_FLUIDMESH_H

#include <vector>
#include <mutex>

#include "CGL/CGL.h"
#include "pointMass.h"

using namespace CGL;
using namespace std;

class Triangle {
public:
    Triangle(Vector3D pm1, Vector3D pm2, Vector3D pm3, Vector3D uv1, Vector3D uv2, Vector3D uv3)
            : pm1(pm1), pm2(pm2), pm3(pm3), uv1(uv1), uv2(uv2), uv3(uv3) {
        normal = -cross(pm2 - pm1, pm3 - pm2);
        normal.normalize();
    }

    // Static references to constituent mesh objects
    Vector3D pm1;
    Vector3D pm2;
    Vector3D pm3;

    Vector3D normal;

    // UV values for each of the points.
    // Uses Vector3D for convenience. This means that the z dimension
    // is not used, and xy corresponds to uv.
    Vector3D uv1;
    Vector3D uv2;
    Vector3D uv3;

    Halfedge *halfedge;
}; // struct Triangle

class Edge {
public:
    Halfedge *halfedge;
}; // struct Edge

class Halfedge {
public:
    Edge *edge;
    Halfedge *next;
    Halfedge *twin;
    Triangle *triangle;
    PointMass *pm;
}; // struct Halfedge

class FluidMesh {
public:
    ~FluidMesh() {}

    std::vector<Triangle> triangles;
    std::mutex mtx;

}; // struct FluidMesh


#endif //CLOTHSIM_FLUIDMESH_H
