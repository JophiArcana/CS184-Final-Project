#include <iostream>
#include <fstream>
#include <nanogui/nanogui.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include "misc/getopt.h" // getopt for windows
#else

#include <getopt.h>
#include <unistd.h>

#endif

#include <unordered_set>
#include <stdlib.h> // atoi for getopt inputs

#include "CGL/CGL.h"
#include "collision/plane.h"
#include "collision/sphere.h"
// #include "cloth.h"
#include "clothSimulator.h"
#include "json.hpp"
#include "misc/file_utils.h"
#include "Fluid.h"

typedef uint32_t gid_t;

using namespace std;
using namespace nanogui;

using json = nlohmann::json;

#define msg(s) cerr << "[ClothSim] " << s << endl;

const string SPHERE = "sphere";
const string PLANE = "plane";
// const string CLOTH = "cloth";
const string FLUID = "fluid";

// const unordered_set<string> VALID_KEYS = {SPHERE, PLANE, CLOTH};
const unordered_set<string> VALID_KEYS = {SPHERE, PLANE, FLUID};

ClothSimulator *app = nullptr;
GLFWwindow *window = nullptr;
Screen *screen = nullptr;

void error_callback(int error, const char *description) {
    puts(description);
}

void createGLContexts() {
    if (!glfwInit()) {
        return;
    }

    glfwSetTime(0);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_SAMPLES, 0);
    glfwWindowHint(GLFW_RED_BITS, 8);
    glfwWindowHint(GLFW_GREEN_BITS, 8);
    glfwWindowHint(GLFW_BLUE_BITS, 8);
    glfwWindowHint(GLFW_ALPHA_BITS, 8);
    glfwWindowHint(GLFW_STENCIL_BITS, 8);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

    // Create a GLFWwindow object
    window = glfwCreateWindow(800, 800, "Cloth Simulator", nullptr, nullptr);
    if (window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        throw std::runtime_error("Could not initialize GLAD!");
    }
    glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM

    glClearColor(0.2f, 0.25f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Create a nanogui screen and pass the glfw pointer to initialize
    screen = new Screen();
    screen->initialize(window, true);

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    glfwSwapInterval(1);
    glfwSwapBuffers(window);
}

void setGLFWCallbacks() {
    glfwSetCursorPosCallback(window, [](GLFWwindow *, double x, double y) {
        if (!screen->cursorPosCallbackEvent(x, y)) {
            app->cursorPosCallbackEvent(x / screen->pixelRatio(),
                                        y / screen->pixelRatio());
        }
    });

    glfwSetMouseButtonCallback(
            window, [](GLFWwindow *, int button, int action, int modifiers) {
                if (!screen->mouseButtonCallbackEvent(button, action, modifiers) ||
                    action == GLFW_RELEASE) {
                    app->mouseButtonCallbackEvent(button, action, modifiers);
                }
            });

    glfwSetKeyCallback(
            window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
                if (!screen->keyCallbackEvent(key, scancode, action, mods)) {
                    app->keyCallbackEvent(key, scancode, action, mods);
                }
            });

    glfwSetCharCallback(window, [](GLFWwindow *, unsigned int codepoint) {
        screen->charCallbackEvent(codepoint);
    });

    glfwSetDropCallback(window,
                        [](GLFWwindow *, int count, const char **filenames) {
                            screen->dropCallbackEvent(count, filenames);
                            app->dropCallbackEvent(count, filenames);
                        });

    glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
        if (!screen->scrollCallbackEvent(x, y)) {
            app->scrollCallbackEvent(x, y);
        }
    });

    glfwSetFramebufferSizeCallback(window,
                                   [](GLFWwindow *, int width, int height) {
                                       screen->resizeCallbackEvent(width, height);
                                       app->resizeCallbackEvent(width, height);
                                   });
}

void usageError(const char *binaryName) {
    printf("Usage: %s [options]\n", binaryName);
    printf("Required program options:\n");
    printf("  -f     <STRING>    Filename of scene\n");
    printf("  -r     <STRING>    Project root.\n");
    printf("                     Should contain \"shaders/Default.vert\".\n");
    printf("                     Automatically searched for by default.\n");
    printf("  -a     <INT>       Sphere vertices latitude direction.\n");
    printf("  -o     <INT>       Sphere vertices longitude direction.\n");
    printf("\n");
    exit(-1);
}

void incompleteObjectError(const char *object, const char *attribute) {
    cout << "Incomplete " << object << " definition, missing " << attribute << endl;
    exit(-1);
}

//bool loadObjectsFromFile(string filename, Cloth *cloth, ClothParameters *cp, vector<CollisionObject *> *objects,
//                         int sphere_num_lat, int sphere_num_lon) {
bool loadObjectsFromFile(string filename, Fluid *fluid, FluidParameters *fp, vector<CollisionObject *> *objects,
                         int sphere_num_lat, int sphere_num_lon) {
    // Read JSON from file
    ifstream i(filename);
    if (!i.good()) {
        return false;
    }
    json j;
    i >> j;

    // Loop over objects in scene
    for (json::iterator it = j.begin(); it != j.end(); ++it) {
        string key = it.key();

        // Check that object is valid
        unordered_set<string>::const_iterator query = VALID_KEYS.find(key);
        if (query == VALID_KEYS.end()) {
            cout << "Invalid scene object found: " << key << endl;
            exit(-1);
        }

        // Retrieve object
        json object = it.value();

        // Parse object depending on type (cloth, sphere, or plane)
        // if (key == CLOTH) {
        if (key == FLUID) {
            // Fluid
            double length, width, height;
            // TODO update this
            int g_length, g_width, g_height;
            int num_particles;
            e_orientation orientation;
            // vector<vector<int>> pinned;

            auto it_length = object.find("length");
            if (it_length != object.end()) {
                length = *it_length;
            } else {
                incompleteObjectError("fluid", "length");
            }

            auto it_width = object.find("width");
            if (it_width != object.end()) {
                width = *it_width;
            } else {
                incompleteObjectError("fluid", "width");
            }


            auto it_height = object.find("height");
            if (it_height != object.end()) {
                height = *it_height;
            } else {
                incompleteObjectError("fluid", "height");
            }

//            auto it_g_length = object.find("g_length");
//            if (it_g_length != object.end()) {
//                g_length = *it_g_length;
//            } else {
//                incompleteObjectError("fluid", "g_length");
//            }
//
//            auto it_g_width = object.find("g_width");
//            if (it_g_width != object.end()) {
//                g_width = *it_g_width;
//            } else {
//                incompleteObjectError("fluid", "g_width");
//            }
//
//            auto it_g_height = object.find("g_length");
//            if (it_g_height != object.end()) {
//                g_height = *it_g_height;
//            } else {
//                incompleteObjectError("fluid", "g_height");
//            }
            /*
            auto it_num_particles = object.find("num_particles");
            if (it_num_particles != object.end()) {
                num_particles = *it_num_particles;
            } else {
                incompleteObjectError("fluid", "num_particles");
            }
            */
//            auto it_orientation = object.find("orientation");
//            if (it_orientation != object.end()) {
//                orientation = *it_orientation;
//            } else {
//                incompleteObjectError("fluid", "orientation");
//            }
/*
            auto it_pinned = object.find("pinned");
            if (it_pinned != object.end()) {
                vector<json> points = *it_pinned;
                for (auto pt: points) {
                    vector<int> point = pt;
                    pinned.push_back(point);
                }
            }
*/
//            fluid->LENGTH = width;
//            fluid->WIDTH = width;
//            fluid->HEIGHT = height;
//            fluid->G_HEIGHT = g_height;
//            fluid->G_LENGTH = g_length;
//            fluid->G_WIDTH = g_width;
            // fluid->num_height_points = num_height_points;
            // cloth->thickness = thickness;
            // fluid->orientation = orientation;
            // cloth->pinned = pinned;

            // Cloth parameters
            double density, rms_velocity, average_distance, molar_mass, kinematic_viscosity, tait_coefficient;

            auto it_density = object.find("density");
            if (it_density != object.end()) {
                density = *it_density;
            } else {
                incompleteObjectError("fluid", "density");
            }

            auto it_rms_velocity = object.find("rms_velocity");
            if (it_rms_velocity != object.end()) {
                rms_velocity = *it_rms_velocity;
            } else {
                incompleteObjectError("fluid", "rms_velocity");
            }

            auto it_average_distance = object.find("average_distance");
            if (it_average_distance != object.end()) {
                average_distance = *it_average_distance;
            } else {
                incompleteObjectError("fluid", "average_distance");
            }

            auto it_molar_mass = object.find("molar_mass");
            if (it_molar_mass != object.end()) {
                molar_mass = *it_molar_mass;
            } else {
                incompleteObjectError("fluid", "molar_mass");
            }

            auto it_kinematic_viscosity = object.find("kinematic_viscosity");
            if (it_kinematic_viscosity != object.end()) {
                kinematic_viscosity = *it_kinematic_viscosity;
            } else {
                incompleteObjectError("fluid", "kinematic_viscosity");
            }

            auto it_tait_coefficient = object.find("tait_coefficient");
            if (it_tait_coefficient != object.end()) {
                tait_coefficient = *it_tait_coefficient;
            } else {
                incompleteObjectError("fluid", "tait_coefficient");
            }

            fp->density = density;
            fp->rms_velocity = rms_velocity;
            fp->average_distance = average_distance;
            fp->molar_mass = molar_mass;
            fp->kinematic_viscosity = kinematic_viscosity;
            fp->tait_coefficient = tait_coefficient;
        } else if (key == SPHERE) {
            Vector3D origin;
            double radius, friction;

            auto it_origin = object.find("origin");
            if (it_origin != object.end()) {
                vector<double> vec_origin = *it_origin;
                origin = Vector3D(vec_origin[0], vec_origin[1], vec_origin[2]);
            } else {
                incompleteObjectError("sphere", "origin");
            }

            auto it_radius = object.find("radius");
            if (it_radius != object.end()) {
                radius = *it_radius;
            } else {
                incompleteObjectError("sphere", "radius");
            }

            auto it_friction = object.find("friction");
            if (it_friction != object.end()) {
                friction = *it_friction;
            } else {
                incompleteObjectError("sphere", "friction");
            }

            Sphere *s = new Sphere(origin, radius, friction, sphere_num_lat, sphere_num_lon);
            objects->push_back(s);
        } else { // PLANE
            Vector3D point, normal;
            double friction;

            auto it_point = object.find("point");
            if (it_point != object.end()) {
                vector<double> vec_point = *it_point;
                point = Vector3D(vec_point[0], vec_point[1], vec_point[2]);
            } else {
                incompleteObjectError("plane", "point");
            }

            auto it_normal = object.find("normal");
            if (it_normal != object.end()) {
                vector<double> vec_normal = *it_normal;
                normal = Vector3D(vec_normal[0], vec_normal[1], vec_normal[2]);
            } else {
                incompleteObjectError("plane", "normal");
            }

            auto it_friction = object.find("friction");
            if (it_friction != object.end()) {
                friction = *it_friction;
            } else {
                incompleteObjectError("plane", "friction");
            }

            Plane *p = new Plane(point, normal, friction);
            objects->push_back(p);
        }
    }

    i.close();

    return true;
}

bool is_valid_project_root(const std::string &search_path) {
    std::stringstream ss;
    ss << search_path;
    ss << "/";
    ss << "shaders/Default.vert";

    return FileUtils::file_exists(ss.str());
}

// Attempt to locate the project root automatically
bool find_project_root(const std::vector<std::string> &search_paths, std::string &retval) {

    for (std::string search_path: search_paths) {
        if (is_valid_project_root(search_path)) {
            retval = search_path;
            return true;
        }
    }
    return false;
}

int main(int argc, char **argv) {
    // Attempt to find project root
    std::vector<std::string> search_paths = {
            ".",
            "..",
            "../..",
            "../../.."
    };
    std::string project_root;
    bool found_project_root = find_project_root(search_paths, project_root);

//    Cloth cloth;
//    ClothParameters cp;


    int c;

    int sphere_num_lat = 40;
    int sphere_num_lon = 40;

    std::string file_to_load_from;
    bool file_specified = false;

    while ((c = getopt(argc, argv, "f:r:a:o:")) != -1) {
        switch (c) {
            case 'f': {
                file_to_load_from = optarg;
                file_specified = true;
                break;
            }
            case 'r': {
                project_root = optarg;
                if (!is_valid_project_root(project_root)) {
                    std::cout
                            << "Warn: Could not find required file \"shaders/Default.vert\" in specified project root: "
                            << project_root << std::endl;
                }
                found_project_root = true;
                break;
            }
            case 'a': {
                int arg_int = atoi(optarg);
                if (arg_int < 1) {
                    arg_int = 1;
                }
                sphere_num_lat = arg_int;
                break;
            }
            case 'o': {
                int arg_int = atoi(optarg);
                if (arg_int < 1) {
                    arg_int = 1;
                }
                sphere_num_lon = arg_int;
                break;
            }
            default: {
                usageError(argv[0]);
                break;
            }
        }
    }

    if (!found_project_root) {
        std::cout << "Error: Could not find required file \"shaders/Default.vert\" anywhere!" << std::endl;
        return -1;
    } else {
        std::cout << "Loading files starting from: " << project_root << std::endl;
    }

    // std::cout << "hello there" << std::endl;

    if (!file_specified) { // No arguments, default initialization
        std::stringstream def_fname;
        def_fname << project_root;
        def_fname << "/scene/pinned2.json";
        file_to_load_from = def_fname.str();
    }

    FluidParameters fp = Fluid::WATER;
    Fluid fluid(2, 2, 1, 2000, fp);
    vector<CollisionObject *> objects;

    bool success = loadObjectsFromFile(file_to_load_from, &fluid, &fp, &objects, sphere_num_lat, sphere_num_lon);
    if (!success) {
        std::cout << "Warn: Unable to load from file: " << file_to_load_from << std::endl;
    }

    glfwSetErrorCallback(error_callback);

    createGLContexts();

    // Initialize the Cloth object
    // fluid.buildFluidMesh();
    fluid.debugFluidMesh();

    // Initialize the ClothSimulator object
    app = new ClothSimulator(project_root, screen);
    // app->loadCloth(&cloth);
    // app->loadClothParameters(&cp);
    app->loadFluid(&fluid);
    // TODO implement
    app->loadFluidParameters(&fp);
    app->loadCollisionObjects(&fluid.collisionObjects);
    // app->loadCollisionObjects(&objects);
    app->init();

    // Call this after all the widgets have been defined

    screen->setVisible(true);
    screen->performLayout();

    // Attach callbacks to the GLFW window

    setGLFWCallbacks();

    cout << "drawing contents" << endl;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        app->drawContents();

        // Draw nanogui
        screen->drawContents();
        screen->drawWidgets();

        glfwSwapBuffers(window);

        if (!app->isAlive()) {
            glfwSetWindowShouldClose(window, 1);
        }
    }

    return 0;
}
