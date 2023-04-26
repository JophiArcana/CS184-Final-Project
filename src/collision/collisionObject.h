#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

#include <nanogui/nanogui.h>

// #include "../clothMesh.h"
#include "../FluidMesh.h"

using namespace CGL;
using namespace std;
using namespace nanogui;

class CollisionObject {
public:
    virtual void render(GLShader &shader) = 0;

    virtual void collide(PointMass &pm, double delta_t) = 0;

private:
    double friction;
};

#endif /* COLLISIONOBJECT */
