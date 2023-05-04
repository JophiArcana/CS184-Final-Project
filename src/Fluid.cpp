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

Vector3D edgeTable[12][2] = {{Vector3D(0,0,0), Vector3D(1,0,0)},
                            {Vector3D(1,0,0), Vector3D(1,1,0)},
                             {Vector3D(1,1,0), Vector3D(0,1,0)},
                             {Vector3D(0,1,0), Vector3D(0,0,0)},
                             {Vector3D(0,0,1), Vector3D(1,0,1)},
                             {Vector3D(1,0,1), Vector3D(1,1,1)},
                             {Vector3D(1,1,1), Vector3D(0,1,1)},
                             {Vector3D(0,1,1), Vector3D(0,0,1)},
                             {Vector3D(0,0,0), Vector3D(0,0,1)},
                             {Vector3D(1,0,0), Vector3D(1,0,1)},
                             {Vector3D(1,1,0), Vector3D(1,1,1)},
                             {Vector3D(0,1,0), Vector3D(0,1,1)}};

Vector3D midpointTable[12] = {Vector3D(0.5,0,0),
                             Vector3D(1,0.5,0),
                             Vector3D(0.5,1,0),
                             Vector3D(0,0.5,0),
                             Vector3D(0.5,0,1),
                             Vector3D(1,0.5,1),
                             Vector3D(0.5,1,1),
                             Vector3D(0,0.5,1),
                             Vector3D(0,0,0.5),
                             Vector3D(1,0,0.5),
                             Vector3D(1,1,0.5),
                            Vector3D(0,1,0.5)};

// taken from https://graphics.stanford.edu/~mdfisher/Code/MarchingCubes/MarchingCubes.cpp
int trianglePosTable[256][16] =
        {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
         {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
         {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
         {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
         {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
         {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
         {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
         {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
         {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
         {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
         {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
         {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
         {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
         {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
         {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
         {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
         {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
         {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
         {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
         {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
         {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
         {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
         {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
         {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
         {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
         {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
         {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
         {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
         {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
         {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
         {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
         {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
         {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
         {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
         {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
         {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
         {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
         {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
         {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
         {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
         {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
         {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
         {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
         {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
         {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
         {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
         {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
         {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
         {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
         {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
         {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
         {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
         {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
         {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
         {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
         {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
         {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
         {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
         {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
         {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
         {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
         {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
         {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
         {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
         {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
         {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
         {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
         {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
         {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
         {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
         {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
         {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
         {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
         {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
         {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
         {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
         {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
         {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
         {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
         {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
         {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
         {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
         {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
         {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
         {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
         {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
         {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
         {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
         {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
         {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
         {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
         {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
         {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
         {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
         {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
         {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
         {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

const struct FluidParameters Fluid::WATER(997, 642.503643481, 0.31E-9, 0.01801528);

Fluid::Fluid(double length, double width, double height, int nParticles, FluidParameters params) :
        LENGTH(length), WIDTH(width), HEIGHT(height), NUM_PARTICLES(nParticles), PARAMS(params), grid_toggle(true) {
    this->VOLUME = length * width * height / 3;
    double true_num_particles = VOLUME * PARAMS.density * 6.022E23 / PARAMS.molar_mass;
    double ratio = true_num_particles / nParticles;

    this->SMOOTHING_RADIUS = 1.8 * PARAMS.average_distance * std::cbrt(ratio);
    this->PARTICLE_MASS = VOLUME * PARAMS.density / nParticles;

    this->SELF_KERNEL = 1. / (PI * std::pow(this->SMOOTHING_RADIUS, 3));
    this->KERNEL_COEFF = 1.5 * this->SELF_KERNEL;
    this->GRAD_KERNEL_COEFF = this->KERNEL_COEFF / (this->SMOOTHING_RADIUS * this->SMOOTHING_RADIUS);
    this->SCORR_COEFF = -0.001 / std::pow(0.65 * this->KERNEL_COEFF, 4);
    this->SELF_SCORR = this->SCORR_COEFF * std::pow(this->SELF_KERNEL, 4);

    this->CELL_SIZE = 2. * this->SMOOTHING_RADIUS;
    this->G_LENGTH = (int) (LENGTH / CELL_SIZE) + 1;
    this->G_WIDTH = (int) (WIDTH / CELL_SIZE) + 1;
    this->G_HEIGHT = (int) (HEIGHT / CELL_SIZE) + 1;
    this->G_SIZE = G_LENGTH * G_WIDTH * G_HEIGHT;

    this->grid1 = new std::vector<PointMass *>[G_SIZE];
    this->grid2 = new std::vector<PointMass *>[G_SIZE];
    this->list = std::vector<PointMass *>();
    this->list.reserve(NUM_PARTICLES);
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
        // Vector3D position = Vector3D(LENGTH * (0.5 * (1 - std::sqrt(0.5)) + std::sqrt(0.5) * uniform_dist(gen)), WIDTH * (0.5 * (1 - std::sqrt(0.5)) + std::sqrt(0.5) * uniform_dist(gen)), HEIGHT * (2 * uniform_dist(gen)));
        Vector3D position = Vector3D(LENGTH * uniform_dist(gen), WIDTH * uniform_dist(gen), HEIGHT * (1. / 3 + uniform_dist(gen) / 3));
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

//    double average_density = 0;
//    for (int index = 0; index < G_SIZE; index++) {
//        std::vector<double> density = this->batch_density(this->batch_W(index));
//        average_density = std::accumulate(density.begin(), density.end(), average_density);
//    }
//    average_density /= NUM_PARTICLES;
//    cout << "Average starting density " << average_density << endl;
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



#ifdef DEBUG
    this->buildDebugFluidMesh();
#else
    this->buildFluidMesh();
#endif


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
    std::vector<PointMass *> *old_grid, *new_grid;
    if (this->grid_toggle) {
        old_grid = this->grid1;
        new_grid = this->grid2;
    } else {
        old_grid = this->grid2;
        new_grid = this->grid1;
    }

    for (int index = 0; index < G_SIZE; index++) {
        new_grid[index].clear();
    }
    for (int index = 0; index < G_SIZE; index++) {
        for (PointMass *pm: old_grid[index]) {
            new_grid[this->get_index(pm->tentative_position)].push_back(pm);
        }
    }
    this->grid_toggle = !this->grid_toggle;
}

double Fluid::helper(Vector3D position, Vector3D corner) const {
    double distance = (position / CELL_SIZE - corner).norm();
    if (distance > 1) {
        return 0;
    }
    return 1 / pow(max(0.1, distance), 3);
}


void Fluid::buildFluidMesh() {
    mesh->triangles.clear();

    int SUBDIVISION = 4;

    vector<int> dx = {0, 0, 0, 0, 1, 1, 1, 1};
    vector<int> dy = {0, 0, 1, 1, 0, 0, 1, 1};
    vector<int> dz = {0, 1, 0, 1, 0, 1, 0, 1};

    vector<int> neighbors = {4, 2, 1};

    // TODO change value
    double THRESHOLD = 15;

    // all vertices outside mesh with adjacent vertex inside the mesh
    vector<Vector3D> surface_vertices = vector<Vector3D>();

#ifdef MULTITHREAD
    auto generate_mesh = [&](int thread_num) {
        for (int subdivision_index = thread_num; subdivision_index < G_SIZE * (int) std::pow(SUBDIVISION, 3); subdivision_index += N_THREADS) {
#else
        for (int subdivision_index = 0;
             subdivision_index < G_SIZE * (int) std::pow(SUBDIVISION, 3); subdivision_index++) {
#endif
            int i = subdivision_index / (G_WIDTH * G_HEIGHT * SUBDIVISION * SUBDIVISION);
            int j = (subdivision_index / (G_HEIGHT * SUBDIVISION)) % (G_WIDTH * SUBDIVISION);
            int k = subdivision_index % (G_HEIGHT * SUBDIVISION);

            int x = i / SUBDIVISION;
            int y = j / SUBDIVISION;
            int z = k / SUBDIVISION;
            int index = z + G_HEIGHT * (y + G_WIDTH * x);
            int bitmap = 0;
            vector<double> pressures = vector<double>(8);

            for (int q = 0; q < 8; q += 1) {
                Vector3D pos = Vector3D(i * 1.0 / SUBDIVISION + dx[q], j * 1.0 / SUBDIVISION + dy[q],
                                        k * 1.0 / SUBDIVISION + dz[q]);

                for (int qqq = 0; qqq < 27; qqq += 1) {
                    int deltaz = qqq % 3 - 1;
                    int deltay = (qqq / 3) % 3 - 1;
                    int deltax = (qqq / 9) % 3 - 1;
                    int nx = x + deltax;
                    int ny = y + deltay;
                    int nz = z + deltaz;
                    if (nx < 0 || ny < 0 || nz < 0) {
                        continue;
                    }
                    if (nx >= G_LENGTH || ny >= G_WIDTH || nz >= G_HEIGHT) {
                        continue;
                    }
                    int index2 = nz + G_HEIGHT * (ny + G_WIDTH * nx);
                    for (PointMass *pm: this->grid()[index2]) {
                        pressures[q] += helper(pm->position, pos);
                    }
                }
//                for (PointMass *pm: this->grid()[index]) {
//                    pressures[q] += helper(pm->position, pos);
//                }
//                // point masses from neighboring cells
//                if (z != 0) {
//                    for (PointMass *pm: this->grid()[index - 1]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }
//                if (z != G_HEIGHT - 1) {
//                    for (PointMass *pm: this->grid()[index + 1]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }
//
//                if (y != 0) {
//                    for (PointMass *pm: this->grid()[index - G_HEIGHT]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }
//                if (y != G_WIDTH - 1) {
//                    for (PointMass *pm: this->grid()[index + G_HEIGHT]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }
//
//                if (x != 0) {
//                    for (PointMass *pm: this->grid()[index - G_HEIGHT * G_WIDTH]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }
//                if (x != G_LENGTH - 1) {
//                    for (PointMass *pm: this->grid()[index + G_HEIGHT * G_WIDTH]) {
//                        pressures[q] += helper(pm->position, pos);
//                    }
//                }

            }
            for (int q = 0; q < 8; q += 1) {
                if (pressures[q] < THRESHOLD) {
                    bitmap |= (1 << q);
                }
            }

            if (bitmap == 0 || bitmap == 255) { // all on one side of the surface
                continue;
            }

            // cout << this->grid()[index].size() << endl;

            int visited = 0;

            for (int q = 0; q < 8; q += 1) {
                if ((visited & (1 << q)) != 0) { // already visited
                    continue;
                }
                if (pressures[q] > THRESHOLD) {
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
                        if (pressures[nextInd] > THRESHOLD) {
                            // vertex on edge
                            double pressure1 = pressures[ind] - THRESHOLD;
                            double pressure2 = pressures[nextInd] - THRESHOLD;
                            if (pressure1 > 0 || pressure2 < 0) {
                                throw std::runtime_error("pressures wrong");
                            }
                            double t = pressure2 / (pressure2 - pressure1);
                            Vector3D point = Vector3D(i, j, k);
                            point[0] += t * dx[ind] + (1 - t) * dx[nextInd];
                            point[1] += t * dy[ind] + (1 - t) * dy[nextInd];
                            point[2] += t * dz[ind] + (1 - t) * dz[nextInd];
                            point[0] *= LENGTH / G_LENGTH / SUBDIVISION;
                            point[1] *= WIDTH / G_WIDTH / SUBDIVISION;
                            point[2] *= HEIGHT / G_HEIGHT / SUBDIVISION;
                            vertices.push_back(point);

                            Vector3D original_point = Vector3D(i, j, k);
                            original_point[0] += dx[ind];
                            original_point[1] += dy[ind];
                            original_point[2] += dz[ind];
                            original_point[0] *= LENGTH / G_LENGTH / SUBDIVISION;
                            original_point[1] *= WIDTH / G_WIDTH / SUBDIVISION;
                            original_point[2] *= HEIGHT / G_HEIGHT / SUBDIVISION;
#ifdef MULTITHREAD
                            mesh->mtx.lock();
                            surface_vertices.push_back(original_point);
                            mesh->mtx.unlock();
#else
                            surface_vertices.push_back(original_point);
#endif
                        } else {
                            queue.push_back(nextInd);
                        }
                    }
                }
                // convert vertices into a mesh

                if (vertices.size() < 3) {
//                        cout << "ERROR " << vertices.size() << endl;
//                        continue;
                } else {
                    // cout << "SUCCESS" << endl;
                }
                // TODO: check if we need to run this 8 times
                // i think we do
                for (int ii = 0; ii <= vertices.size() - 3; ii += 1) { // TODO change uv values of triangle
                    // 0, ii + 1, ii + 2
                    Triangle triangle(vertices[0], vertices[ii + 1], vertices[ii + 2], vertices[0],
                                      vertices[ii + 1], vertices[ii + 2]);
#ifdef MULTITHREAD
                    mesh->mtx.lock();
                    mesh->triangles.emplace_back(triangle);
                    mesh->mtx.unlock();
#else
                    mesh->triangles.emplace_back(triangle);
#endif
                }

            }
        }
        // I think this can maybe be done in the "for (int toXor: neighbors)" for loop to save pressure recalculation
        for (int vt = 0; vt <= surface_vertices.size(); vt += 1) {
            for (int i = -1; i < 1; i++) {
                for (int j = -1; j < 1; j++) {
                    for (int k = -1; k < 1; k++) {

                    }
                }
            }
            // TODO: finish this out
            // 0, ii + 1, ii + 2
            // Triangle triangle(vertices[0], vertices[ii + 1], vertices[ii + 2], vertices[0],
            //                  vertices[ii + 1], vertices[ii + 2]);
            // mesh->triangles.emplace_back(triangle);
        }
#ifdef MULTITHREAD
    };
    std::thread threads[N_THREADS];
    for (int thread_num = 0; thread_num < N_THREADS; thread_num++) {
        threads[thread_num] = std::thread(generate_mesh, thread_num);
    }
    for (auto &thread: threads) {
        thread.join();
    }
#endif
}


void Fluid::buildDebugFluidMesh() {
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

