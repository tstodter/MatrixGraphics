/*********************************************************************
 *  CMPSC 457                                                        *
 *  HW 4                                           *
 *  March 27, 2013                                                   *
 *  Travis Stodter                                                    *
 *  tms5274@psu.edu                                                    *                                                                
 *********************************************************************/  


#include <GL/glut.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "matrices.h"


// for your convenience while debugging
using std::cout;
using std::cerr;
using std::endl;


// glut callbacks
void reshape(int w, int h);
void display();
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char key, int x, int y);


// helpers
void init();
void initObj();
void initCam();
void drawFaces();
void DeviceToWorld(double u, double v, double& x, double& y);
void reset();
bool isBackface(int i);


// transformations
void Rotate(double dx, double dy);
void Translate_xy(double tx, double ty);
void Translate_xz(double tx, double ty);
void Scale(double s);
void rotateCamera(double dx, double dy);


// projection
void SetViewMatrix();
void SetOrthoMatrix();
void SetPerspMatrix();


// default device window size
int win_w = 512;
int win_h = 512;
Matrix4 Mvp = {{ {win_h / 2.0, 0, 		    0, (win_h - 1) / 2.0}, 
			  {0, 		 win_h / 2.0, 0, (win_h - 1) / 2.0}, 
			  {0, 		 0, 		    1, 0                }, 
			  {0, 		 0, 		    0, 1                } 
		    }};


// for tracking mouse events
struct MouseTracker
{
  int modifiers;
  int button;
  double initx, inity;
  double finalx, finaly;
};

MouseTracker mtracker;


// for camera parameters
struct Camera
{
  bool perspective;               /* projection method */
  bool backface;			    /* backface culling toggle */
  double l, r, b, t, n, f;        /* view volume */
  Point3 eye;                     /* eye position */
  Vector3 u, v, w;                /* eye coordinate system */
  Matrix4 Mo;                     /* orthographic projection */
  Matrix4 Mv;                     /* view matrix for arbitrary view*/
  Matrix4 Mp;                     /* perspective matrix */
};

Camera cam;


// for objects
const int MAXNAMESTRING = 20;
const int MAXVERTICES = 1000;
const int MAXEDGES = 500;
const int MAXFACES = 50;

struct Object3D {
  char name[MAXNAMESTRING];       /* The name of object for printing */
  int Nvertices;                  /* number of vertices */
  int Nfaces ;                    /* number of faces */
  Matrix4 frame;                  /* the local to world coord transform */
  Point3 center;                  /* center of mass */
  HPoint3 vertices[MAXVERTICES];  /* coodrdinates of each vertex */
  int faces[MAXFACES][6];         /* If face has N vertices, give N + 1 
			  	     numbers -> first the number of vertices
			  	     in the face, then the index numbers of
				     each vertex as it appears in the 
				     "vertices"  array. */
};



// Note: We will keep the initial coordinates of the vertices
//       as originally given.  In other words, we will not change
//       the given coordinates of the vertices, even after any
//       transformation.  All the transformation will be recorded
//       in frame.  This way, you can reset the object to the
//       original position at any time, even after applying many
//       transformations.

Object3D obj = {
  "house", 10, 7,
  
  // initially identity matrix (no transformations applied yet)
  {{ {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0},
     {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}  }},
  
  // center of the object is at origin
  {0.0, 0.0, 0.0},
  
  // vertices of the object in no particular order
  {  {0.0, 1.0, 2.0, 1.0},    {-1.0, 0.5, 2.0, 1.0},
     {-1.0, -1.0, 2.0, 1.0},  {1.0, -1.0, 2.0, 1.0},
     {1.0, 0.5, 2.0, 1.0},    {0.0, 1.0, -2.0, 1.0},
     {-1.0, 0.5, -2.0, 1.0},  {-1.0, -1.0, -2.0, 1.0},
     {1.0, -1.0, -2.0, 1.0},  {1.0, 0.5, -2.0, 1.0}   },
  
  // faces
  {  {5,   0, 1, 2, 3, 4},
     {5,   9, 8, 7, 6, 5},
     {4,   4, 3, 8, 9},
     {4,   0, 4, 9, 5},
     {4,   1, 0, 5, 6},
     {4,   2, 1, 6, 7},
     {4,   3, 2, 7, 8}    }
};




// OpenGL/glut programs will have the structure shown here
//    although with different args and callbacks.
//
// You should not modify main().
// If you really want to modify it, do it at your own risk.
//
// For complete description of each glut functions used, see
// glut manual page on class website.

int main(int argc, char *argv[])
{
  // initialize glut
  glutInit(&argc, argv);

  // use double buffering with RGB colors
  // double buffer removes most of the flickering
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // set window size and position
  glutInitWindowSize(win_w, win_h);
  glutInitWindowPosition(100, 100);

  // now, create window with title "Viewing"
  glutCreateWindow("Viewing");


  // other stuffs like background color, viewing, etc will be
  // set up in this function.
  init();
  
  // initialize (arrange) the object
  initObj();

  // initialize the camera
  initCam();

  
  // register callbacks for glut
  glutDisplayFunc(display);   // for display
  glutReshapeFunc(reshape);   // for window move/resize
  glutMouseFunc(mouse);       // for mouse buttons
  glutMotionFunc(motion);     // for mouse movement while mouse button pressed
  glutKeyboardFunc(keyboard); // for keyboard


  // start event processing, i.e., accept user inputs
  glutMainLoop();

  return 0;
}


//
// implementation for glut callbacks
//

// called when the window is resized/moved (and some other cases)
void reshape(int w, int h)
{
  // change window size
  win_w = w;
  win_h = h;

  // set the new viewport
  glViewport(0, 0, (GLint)win_w, (GLint)win_h);

  // we will use orthographic projection when drawing the object.
  //
  // NOTE: This has nothing to do with the projections you are
  //       to implement in this assignment.  We only need this
  //       when you draw 2D lines.  In other words, find the 2D
  //       projections of the end points of a given 3D line using
  //       the projection matrices you implemented and then draw
  //       a 2D line between the projected end-points.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 511.0, 0.0, 511.0, -1.0, 1.0);
}


// called when the window needs to be redrawn
void display()
{
  // clear the buffer with bg color set in init()
  // you can think of the buffer as a raster array provided by GL
  glClear(GL_COLOR_BUFFER_BIT);

  // draw the object on the buffer you just cleared
  drawFaces();

  // swap the buffers.
  // we are using 2 buffers provided by GL (see main) -- double buffer.
  // they are called front / back buffers.
  // what you see on the screen is the content of front buffer
  // what you clear/draw above is done on back buffer
  // once drawing is done on the back buffer, 
  //       you need to display the content of the back buffer.
  // swapping buffers means swapping back buffer with front buffer
  //       so that front buffer becomes back buffer and
  //       back buffer becomes front buffer.
  // once back buffer becomes front buffer, the content of it will be
  //       drawn on the screen, so you can see it.
  glutSwapBuffers();
}


// called when a mouse event (button pressed/released) occurs in glut, 
//     mouse buttons are represented as
//           GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, and GLUT_RIGHT_BUTTON
//     status mouse buttons are represented as
//           GLUT_UP and GLUT_DOWN
// 
void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {  // mouse pressed.  record the details
    // which button?
    mtracker.button = button;
    // any modifiers (keys like shift/ctrl/alt) pressed?
    mtracker.modifiers = glutGetModifiers();
    // mouse position in world
    DeviceToWorld(double(x), double(y), mtracker.initx, mtracker.inity);
  }
}


// called when a mouse moves with a button pressed
void motion(int x, int y)
{
  // get the mouse position in world
  DeviceToWorld(double(x), double(y), mtracker.finalx, mtracker.finaly);

  // now, process the user input, i.e., mouse movement
  switch (mtracker.button) {
  case GLUT_LEFT_BUTTON:
    if (mtracker.modifiers & GLUT_ACTIVE_SHIFT) {
      // shift + left button
      Translate_xy(mtracker.finalx - mtracker.initx,
		   mtracker.finaly - mtracker.inity);
    }
    else if (mtracker.modifiers & GLUT_ACTIVE_CTRL) {
      // CTRL + left button
	 rotateCamera(mtracker.finalx - mtracker.initx,
			 		mtracker.finaly - mtracker.inity);
    }
    else {
      // left button
      Translate_xz(mtracker.finalx - mtracker.initx,
		   mtracker.finaly - mtracker.inity);
    }
    break;
  case GLUT_RIGHT_BUTTON:
    Scale(mtracker.finalx - mtracker.initx);
    break;
  case GLUT_MIDDLE_BUTTON:
    Rotate(mtracker.finalx - mtracker.initx,
	   mtracker.finaly - mtracker.inity);
    break;
  }
  
  // redisplay after transformation
  glutPostRedisplay();

  // reset the mouse position
  mtracker.initx = mtracker.finalx;
  mtracker.inity = mtracker.finaly;
}  


// called when a keyboard event (key typed) occurs
void keyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'Q':  // quit the program
  case 'q':
    exit(0);
  case 'P':  // toggle the projection method
  case 'p':  // between orthographic and perspective projections
    cam.perspective = !cam.perspective;
    SetPerspMatrix();
    glutPostRedisplay();
    break;
  case 'R':  //Reset the object
  case 'r':
	reset();	
	glutPostRedisplay();
	break;
  case 'B':
  case 'b':
     cam.backface = !cam.backface;
	glutPostRedisplay();
  }
}


//
// implementation for helpers
//

void init()
{
  // set background color to black
  glClearColor(0.0, 0.0, 0.0, 0.0);
}


// arrange the object to its initial position
// NOTE: you may change the parameters as you like
void initObj()
{
  Vector3 n;

  // rotate around y-axis
  n.x = 0.0;  n.y = 1.0;  n.z = 0.0;
  double angle = M_PI / 6.0;
  Matrix4 m1 = SetRotMatrix(n, angle);

  // rotate around x-axis
  n.x = 1.0;  n.y = 0.0;  n.z = 0.0;
  angle = M_PI / 6.0;
  Matrix4 m2 = SetRotMatrix(n, angle);

  // translate so that the object is inside view volume
  // (see initCam() for the view volume)
  Matrix4 m3 = SetTransMatrix(0.0, 0.0, -5.0);

  // notice the order of the transformations applied
  //  i.e., Ry -> Rx -> T  becomes  TRxRy in matrix multiplication
  obj.frame = Mult4(m3, Mult4(m2, Mult4(m1, obj.frame)));
}


// initialize camera parameters
// NOTE: you may change the parameters as you like
void initCam()
{
  // use orthographic projection as default
  cam.perspective = false;

  // no back face culling by default
  cam.backface = false;

  // camera position
  cam.eye.x = 0.0;
  cam.eye.y = 0.0;
  cam.eye.z = 0.0;

  // view volume
  cam.l = -5.0;  cam.r = 5.0;
  cam.b = -5.0;  cam.t = 5.0;
  cam.n = -1.0;  cam.f = -6.0;

  // camera coordinate system
  cam.u.x = 1.0;  cam.u.y = 0.0;  cam.u.z = 0.0;
  cam.v.x = 0.0;  cam.v.y = 1.0;  cam.v.z = 0.0;
  cam.w.x = 0.0;  cam.w.y = 0.0;  cam.w.z = 1.0;

  // set Mv, Mp, Mo
  SetViewMatrix(); 
  SetPerspMatrix();
  SetOrthoMatrix();
}


// draw object faces
// IMPORTANT: you are only allowed to use glVertex2d.
void drawFaces()
{
  Matrix4 overallMat = Mult4( Mvp, cam.Mo);
  overallMat = Mult4(overallMat, cam.Mp);
  overallMat = Mult4( overallMat, Mult4(cam.Mv, obj.frame) );

  int i = 0;
  //For each face in the obj
  while (i < MAXFACES && i < obj.Nfaces) {
	
	if ( !cam.backface || isBackface(i) ) {

		int numVertices = obj.faces[i][0];
		glBegin(GL_LINE_LOOP);
			// For each vertex in the face
			for (int j = 1; j < numVertices + 1; j++) {
				
				HPoint3 point = obj.vertices[obj.faces[i][j]];
				point = TransHPoint3(overallMat, point);
				point = Homogenize(point);

				glVertex2d( point.x, point.y );
			}
		glEnd();
	}
	i++;
  }
}


// convert device coordinate to world coordinate
void DeviceToWorld(double u, double v, double& x, double& y)
{

	//Simplifies from u / win_w = [x - (-1)] / [1 - (-1)]
	//where 1 and -1 refer to the range of x values in world coordinates

	x = (2 * u / win_w) - 1;	
	y = (2 * v / win_h) - 1;

}

//
// implementation for transformations
//


// rotation using the Rolling Ball transformation
void Rotate(double dx, double dy)
{
  double dr = sqrt( dx*dx + dy*dy );

  Vector3 n;
  n.x = -dy / dr;
  n.y = dx / dr;
  n.z = 0;

  double theta = atan(dr / 1);

  obj.frame = Mult4( obj.frame, SetRotMatrix(n, theta) ) ;
}


// translation in xy-plane
void Translate_xy(double tx, double ty)
{
  Matrix4 translateM = SetTransMatrix(tx, -ty, 0);

  obj.frame = Mult4( obj.frame, translateM );
}


// translation in xz-plane
void Translate_xz(double tx, double ty)
{
  Matrix4 translateM = SetTransMatrix(tx, 0, -ty);

  obj.frame = Mult4( obj.frame, translateM );
}


// uniform scale
void Scale(double sx)
{
  Matrix4 scaleM = SetScaleMatrix(1 + .5 * sx, 
						    1 + .5 * sx,
						    1 + .5 * sx );

  obj.frame = Mult4( obj.frame, scaleM);
}

//
// implementation for projection
//

// Mv
void SetViewMatrix()
{

  Matrix4 tempM = {{ {cam.u.x, cam.u.y, cam.u.z, 0},
				 {cam.v.x, cam.v.y, cam.v.z, 0},
				 {cam.w.x, cam.w.y, cam.w.z, 0},
				 {0	    , 0      , 0      , 1} }};

  Matrix4 transM = SetTransMatrix(-cam.eye.x, -cam.eye.y, -cam.eye.z);

  cam.Mv = Mult4( tempM, transM );
}

// Mo
void SetOrthoMatrix()
{
  double rMinusl = cam.r - cam.l;
  double rPlusl = cam.r + cam.l;
  double tMinusb = cam.t - cam.b;
  double tPlusb = cam.t + cam.b;
  double nMinusf = cam.n - cam.f;
  double nPlusf = cam.n + cam.f;


  Matrix4 m;
  m.elem[0][0] = 2.0 / rMinusl ;  m.elem[0][1] =  0.0;  		 
		m.elem[0][2] = 0.0;  		   m.elem[0][3] = -rPlusl / rMinusl ;
  m.elem[1][0] =  0.0;  		    m.elem[1][1] = 2.0 / tMinusb ;  
		m.elem[1][2] = 0.0;  		   m.elem[1][3] = -tPlusb / tMinusb;
  m.elem[2][0] =  0.0;  		    m.elem[2][1] =  0.0;  		 
		m.elem[2][2] = 2.0 / nMinusf ;  m.elem[2][3] = -nPlusf / nMinusf;
  m.elem[3][0] =  0.0; 		    m.elem[3][1] =  0.0;  		 
		m.elem[3][2] = 0.0;  		   m.elem[3][3] = 1.0;

  cam.Mo = m;	
}


// Mp
void SetPerspMatrix()
{
  if (cam.perspective == true) {
	
	Matrix4 newMp = {{  {cam.n, 0, 	0, 		   	0			  },
			   		{0, 	   cam.n, 0, 		   	0			  },
			 	     {0, 	   0, 	cam.n + cam.f, -cam.f * cam.n},
			  		{0, 	   0, 	1, 		   	0			  } }};
	cam.Mp = newMp;
  }  
  else
  	cam.Mp = IdentMatrix();
}

/* Resets the scene */
void reset() {
  obj.frame = IdentMatrix();
  initObj();
}

/* Returns a boolean indicating whether or not face i is currently back face
   to the camera.												  */
bool isBackface(int i) {
	Matrix4 overallMat = Mult4( Mvp, cam.Mo);
 	overallMat = Mult4(cam.Mo, cam.Mp);
  	overallMat = Mult4( overallMat, Mult4(cam.Mv, obj.frame) );

	HPoint3 p1 = obj.vertices[ obj.faces[i][1] ];
	HPoint3 p2 = obj.vertices[ obj.faces[i][2] ];
	HPoint3 p3 = obj.vertices[ obj.faces[i][3] ];

	// Translate to world coordinates
	p1 = TransHPoint3(overallMat, p1);
	p2 = TransHPoint3(overallMat, p2);
	p3 = TransHPoint3(overallMat, p3);

	p1 = Homogenize(p1);
	p2 = Homogenize(p2);
	p3 = Homogenize(p3);

	Vector3 v1, v2;
	v1.x = p2.x - p1.x;
	v1.y = p2.y - p1.y;
	v1.z = p2.z - p1.z;

	v2.x = p3.x - p1.x;
	v2.y = p3.y - p1.y;
	v2.z = p3.z - p1.z;

	Vector3 xProd = { v1.y*v2.z - v1.z*v2.y,
				   v1.x*v2.z - v1.z*v2.x,
				   v1.x*v2.y - v1.y*v2.x  };
	
	return xProd.z > 0;	
}

/* Intended to rotate the camera around the initial object position,
   but does not work.										 */
void rotateCamera(double dx, double dy) {
  
	//Change position of camera
  double dr = sqrt( dx*dx + dy*dy );

  Vector3 n;
  n.x = -(-dy) / dr;
  n.y = (-dx) / dr;
  n.z = 0;

  double theta = atan(dr / 1);

  Matrix4 rotMatrix = SetRotMatrix(n, theta);
  HPoint3 newEyePos = TransHPoint3(rotMatrix, Pt3toHPt3(cam.eye));
  cam.eye = HPt3toPt3(newEyePos);

  //Get gaze vector
  Vector3 gaze;
  gaze.x = -cam.eye.x;
  gaze.y = -cam.eye.y;
  gaze.z = -5 - cam.eye.z;

  //Get new eye axes
  Vector3 u, v, w;
  
  double sizeOfW = sqrt(gaze.x*gaze.x + gaze.y*gaze.y + gaze.z*gaze.z);
  w.x = -gaze.x / sizeOfW;
  w.y = -gaze.y / sizeOfW;
  w.z = -gaze.z / sizeOfW;

  u = cam.u;

  v.x = w.y*u.z - w.z*u.y;
  v.y = w.x*u.z - w.z*u.x;
  v.z = w.x*u.y - w.y*u.x;

  double sizeOfV = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

  v.x = v.x / sizeOfV;
  v.y = v.y / sizeOfV;
  v.z = v.z / sizeOfV;

  cam.u = u;
  cam.v = v;
  cam.w = w;

  SetViewMatrix();
 
}


