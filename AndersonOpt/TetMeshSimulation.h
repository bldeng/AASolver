#ifdef USE_GL_FOLDER
#include <GL/glut.h>
#else
#include <glut.h>
#endif
#include "SimulationOpt.h"
#include "MeshTypes.h"
#include "FileParser.h"
#include "TetMeshIO.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

#define NO_MOTION			0
#define ZOOM_MOTION			1
#define ROTATE_MOTION		2
#define TRANSLATE_MOTION	3

bool idle_run = false;

SimulationOpt sOpt;

Matrix3X problem_def_poisitions;
Eigen::Matrix4Xi tet_connectivity;
Eigen::Matrix3Xi tet_boundary;
std::vector<int> handle_indices;

int iter;
bool use_fast_svd = false;
bool is_initialized = false;

///////////////////////////////////////////////////////////////////////////////////////////
//  Simulation entry here...
///////////////////////////////////////////////////////////////////////////////////////////
void Update() {
  sOpt.simulate(iter);
}

///////////////////////////////////////////////////////////////////////////////////////////
//  Get handle index
///////////////////////////////////////////////////////////////////////////////////////////
bool read_handle_files(const char *index_filename,
                       std::vector<int> &handle_indices) {
  if (!read_sequence_file(index_filename, handle_indices)) {
    std::cerr << "Error reading coordinate files" << std::endl;
    return false;
  }

  if (handle_indices.empty()) {
    std::cerr << "Error: no indices read" << std::endl;
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
//  class OPENGL_DRIVER
///////////////////////////////////////////////////////////////////////////////////////////
class OPENGL_DRIVER {
 public:
  //3D display configuration.
  static int screen_width, screen_height;
  static int mesh_mode;
  static int render_mode;
  static double zoom, swing_angle, elevate_angle;
  static double center[3];
  static int motion_mode, mouse_x, mouse_y;

  OPENGL_DRIVER(int *argc, char **argv) {
    if (!is_initialized) {
      if (!TetMeshIO::load(argv[1], problem_def_poisitions, tet_connectivity,
                           tet_boundary)) {
        std::cerr
            << "Error: unable to read problem definition mesh from the file "
            << argv[1] << std::endl;
        return;
      }

      if (!read_handle_files(argv[2], handle_indices)) {
        return;
      }

      int n_vtx = problem_def_poisitions.cols();

      for (int i = 0; i < static_cast<int>(handle_indices.size()); ++i) {
        if (handle_indices[i] < 0 || handle_indices[i] >= n_vtx) {
          std::cerr << "Error: invalid handle index: " << handle_indices[i]
                    << std::endl;
          return;
        }
      }

      Parameters param;
      if (!param.load(argv[3])) {
        std::cerr << "Error: unable to load option file " << argv[1]
                  << std::endl;
        return;
      }
      if (!param.valid_parameters()) {
        std::cerr << "Invalid filter options. Aborting..." << std::endl;
        return;
      }
      param.output();

      iter = param.iter;

      sOpt.setTetMessage(tet_connectivity, tet_boundary, handle_indices);
      sOpt.initSolver(problem_def_poisitions, use_fast_svd, param);

      is_initialized = true;
    }

    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(screen_width, screen_height);
    glutCreateWindow("Tetrahegren Simulation");
    glutDisplayFunc(Handle_Display);
    glutReshapeFunc(Handle_Reshape);
    glutKeyboardFunc(Handle_Keypress);
    glutMouseFunc(Handle_Mouse_Click);
    glutMotionFunc(Handle_Mouse_Move);
    glutSpecialFunc(Handle_SpecialKeypress);
    glutIdleFunc(Handle_Idle);
    Handle_Reshape(screen_width, screen_height);

    glutMainLoop();
  }

  static void Handle_Display() {
    glLoadIdentity();
    glClearColor(0.8, 0.8, 0.8, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    gluLookAt(0, 0, zoom, 0, 0, 0, 0, 1, 0);
    glRotated(elevate_angle, 1, 0, 0);
    glRotated(swing_angle, 0, 1, 0);
    glTranslatef(-center[0], -center[1], -center[2]);

    sOpt.Render();
    glutSwapBuffers();
  }

  static void Handle_Idle() {
    if (idle_run)
      Update();
    glutPostRedisplay();
  }

  static void Handle_Reshape(int w, int h) {
    screen_width = w, screen_height = h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(4, (double) w / (double) h, 1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glEnable(GL_DEPTH_TEST);
    {
      GLfloat LightDiffuse[] = { 1.0, 1.0, 1.0, 1 };
      GLfloat LightPosition[] = { 0, 0, -100 };
      glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
      glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
      glEnable(GL_LIGHT0);
    }
    {
      GLfloat LightDiffuse[] = { 1.0, 1.0, 1.0, 1 };
      GLfloat LightPosition[] = { 0, 0, 100 };
      glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
      glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
      glEnable(GL_LIGHT1);
    }
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glutPostRedisplay();
  }

  static void Handle_Keypress(unsigned char key, int mousex, int mousey) {
    switch (key) {
      case 27:
        exit(0);
      case 'a': {
        zoom -= 2;
        if (zoom < 0.3)
          zoom = 0.3;
        break;
      }
      case 'z': {
        zoom += 2;
        break;
      }
      case 't': {
        render_mode = (render_mode + 1) % 4;
        break;
      }
      case 's': {
        int n = 1;
        for (int i = 0; i < n; i++) {
          Update();
        }
        break;
      }
      case '1': {
        idle_run = true;
        break;
      }
      case 'e': {
        idle_run = false;
        break;
      }
    }
    glutPostRedisplay();
  }

  static void Handle_SpecialKeypress(int key, int x, int y) {
    if (key == 100)
      swing_angle += 3;
    else if (key == 102)
      swing_angle -= 3;
    else if (key == 103)
      elevate_angle -= 3;
    else if (key == 101)
      elevate_angle += 3;
    Handle_Reshape(screen_width, screen_height);
    glutPostRedisplay();
  }

  static void Handle_Mouse_Move(int x, int y) {
    if (motion_mode != NO_MOTION) {
      if (motion_mode == ROTATE_MOTION) {
        swing_angle += (double) (x - mouse_x) * 360 / (double) screen_width;
        elevate_angle += (double) (y - mouse_y) * 180 / (double) screen_height;
        if (elevate_angle > 90)
          elevate_angle = 90;
        else if (elevate_angle < -90)
          elevate_angle = -90;
      }
      if (motion_mode == ZOOM_MOTION)
        zoom += 0.05 * (y - mouse_y);
      if (motion_mode == TRANSLATE_MOTION) {
        center[0] -= 0.1 * (mouse_x - x);
        center[2] += 0.1 * (mouse_y - y);
      }
      mouse_x = x;
      mouse_y = y;
      glutPostRedisplay();
    }
  }

  static void Handle_Mouse_Click(int button, int state, int x, int y) {
    if (state == GLUT_UP) {
      motion_mode = NO_MOTION;
    }
    if (state == GLUT_DOWN) {
      // Set up the motion target
      int modif = glutGetModifiers();
      if (modif & GLUT_ACTIVE_SHIFT)
        motion_mode = ZOOM_MOTION;
      else if (modif & GLUT_ACTIVE_CTRL)
        motion_mode = TRANSLATE_MOTION;
      else
        motion_mode = ROTATE_MOTION;
      mouse_x = x;
      mouse_y = y;
    }

    glutPostRedisplay();
  }
};

int OPENGL_DRIVER::screen_width = 1024;
int OPENGL_DRIVER::screen_height = 768;
int OPENGL_DRIVER::mesh_mode = 0;
int OPENGL_DRIVER::render_mode = 0;
double OPENGL_DRIVER::zoom = 30;
double OPENGL_DRIVER::swing_angle = -0;
double OPENGL_DRIVER::elevate_angle = 0;
double OPENGL_DRIVER::center[3] = { 0, 0, 0 };
int OPENGL_DRIVER::motion_mode = NO_MOTION;
int OPENGL_DRIVER::mouse_x = 0;
int OPENGL_DRIVER::mouse_y = 0;
