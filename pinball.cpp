
////////////////////////////////////////////////////
// pinball.cpp
// Template code for drawing a pinball simulator.
////////////////////////////////////////////////////

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#ifndef WIN32
#include <GL/glut.h>
#endif

#ifdef WIN32
#include "GL/glut.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <deque>
#include <queue>
#include <vector>

#include "Timer.h"

Timer TM;

int Width = 480 ;
int Height = 480 ;
float Time = 0;

void set_colour(float r, float g, float b) ;

#define pi 3.141592653589793
#define moveAngle (pi / 24)
void moveViewVertical(double eye[], int down);

#define X 0
#define Y 1
#define Z 2

// The eye point and look-at point.
int currentview = 0;
#define NVIEWS 5
double views[NVIEWS][3] = {{0, 10.0, 8.0},  {-5, 10.0, 8.0}, 
                           {0, 10.0, -8.0}, {0, 0, 8.0}, {-5, 10.0, 8.0}};
double eye[3] = {views[0][0], views[0][1], views[0][2]};

double ref[3] = {0.0, 0.0, 0.0};

void (*normfns[])(GLenum) = {glEnable, glDisable};
int nfp;

/////////////////////////////////////////////////////
//    PROC: drawCylinder()
//    DOES: this function 
//			render a solid cylinder  oriented along the Z axis. Both bases are of radius 1. 
//			The bases of the cylinder are placed at Z = 0, and at Z = 1.
//	 
//          
// Don't change.
//////////////////////////////////////////////////////
void drawCylinder(void)
{
  static GLUquadricObj *cyl = NULL ;
  if( cyl == NULL )
    {
      cyl = gluNewQuadric() ;
    }
  if( cyl == NULL )
    {
      fprintf(stderr,"Cannot allocate cylinder.\n") ;
      return ;
    }
  gluQuadricDrawStyle(cyl,GLU_FILL) ;
  gluQuadricNormals(cyl,GLU_SMOOTH) ;
  gluCylinder(cyl,1.0,1.0,1.0,10,10) ;
  
}

//////////////////////////////////////////////////////
//    PROC: drawCone()
//    DOES: this function 
//			render a solid cone oriented along the Z axis with base radius 1. 
//			The base of the cone is placed at Z = 0, and the top at Z = 1. 
//         
// Don't change.
//////////////////////////////////////////////////////
void drawCone(void)
{
  glutSolidCone(1,1,20,20) ;
}


//////////////////////////////////////////////////////
//    PROC: drawCube()
//    DOES: this function draws a cube with dimensions 1,1,1
//          centered around the origin.
// 
// Don't change.
//////////////////////////////////////////////////////

void drawCube(void)
{
  glutSolidCube(1.0) ;
}


//////////////////////////////////////////////////////
//    PROC: drawSphere()
//    DOES: this function draws a sphere with radius 1
//          centered around the origin.
// 
// Don't change.
//////////////////////////////////////////////////////

void drawSphere(void)
{ 
  glutSolidSphere(1.0, 50, 50) ;
}


//////////////////////////////////////////////////////
//    PROC: myKey()
//    DOES: this function gets called for any keypresses
// 
//////////////////////////////////////////////////////

void myKey(unsigned char key, int x, int y)
{
  float time ;
  double oldCoord;
  switch (key) {
  case 'q':
    exit(0); 
    break ;

  case 'v':
    printf("pinball: changing view\n");
    currentview = (currentview + 1) % NVIEWS;
    eye[0] = views[currentview][0];
    eye[1] = views[currentview][1];
    eye[2] = views[currentview][2];
    printf("pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
    break;
	    
  case 'n':
    normfns[nfp++%2](GL_AUTO_NORMAL);
    printf (nfp%2 ? "one\n" : "two\n");
    break;
    
	// EXTRA CREDIT FUNCTIONALITY:
	// similar to vim, h, j, k, and l move the view to the left, right
	// up and down.
        case 'h':
            printf("robot: changing view\n");
	    oldCoord = eye[0];
            eye[0] = eye[0] * cos(moveAngle) - eye[2] * sin (moveAngle);
            eye[2] = oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
            printf("robot: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'j':
            printf("robot: changing view\n");
	    moveViewVertical(eye, 1);
            printf("robot: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'k':
            printf("robot: changing view\n");
	    moveViewVertical(eye, -1);
            printf("robot: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'l':
            printf("robot: changing view\n");
	    oldCoord = eye[0];
            eye[0] = eye[0] * cos(moveAngle) + eye[2] * sin (moveAngle);
            eye[2] = -oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
            printf("robot: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
	// n and m zoom in and zoom out
	// NOTE: this only works with a perspective projection, not the
	// default orthogonal projection!
	case 'u':
	    printf("robot: changing view\n");
	    eye[0] /= 1.1;
	    eye[1] /= 1.1;
	    eye[2] /= 1.1;
	    break;
	case 'i':
	    printf("robot: changing view\n");
	    eye[0] *= 1.1;
	    eye[1] *= 1.1;
	    eye[2] *= 1.1;
	    break;
  }

  glutPostRedisplay() ;
}


/*********************************************************
    PROC: myinit()
    DOES: performs most of the OpenGL intialization
     -- change these with care, if you must.
**********************************************************/

void myinit(void)
{
  GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat position[] = { 0.0, 0.0, 3.0, 1.0 };
    
  GLfloat lmodel_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat local_view[] = { 0.0 };

  /**** set lighting parameters ****/
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE) ;

  /*    glFrontFace (GL_CW); */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
	
  glPixelStorei(GL_PACK_ALIGNMENT,1) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT,1) ;
  glShadeModel(GL_SMOOTH) ;

}

/*********************************************************
    PROC: set_colour();
    DOES: sets all material properties to the given colour
    -- don't change
**********************************************************/

void set_colour(float r, float g, float b)
{
  float ambient = 0.2f;
  float diffuse = 0.7f;
  float specular = 1.4f;
  GLfloat mat[4];

  /**** set ambient lighting parameters ****/
  mat[0] = ambient*r;
  mat[1] = ambient*g;
  mat[2] = ambient*b;
  mat[3] = 1.0;
  glMaterialfv (GL_FRONT, GL_AMBIENT, mat);

  /**** set diffuse lighting parameters ******/
  mat[0] = diffuse*r;
  mat[1] = diffuse*g;
  mat[2] = diffuse*b;
  mat[3] = 1.0;
  glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);

  /**** set specular lighting parameters *****/
  mat[0] = specular*r;
  mat[1] = specular*g;
  mat[2] = specular*b;
  mat[3] = 1.0;
  glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
  glMaterialf (GL_FRONT, GL_SHININESS, 1);
}


/*********************************************************
    PROC: draw_axes();
    DOES: draws unit vectors for the x, y, and z axes,
					colored red, green, and blue, respectively.
**********************************************************/
void draw_axes()
{
  // x-axis
  set_colour(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(1.0,0.0,0.0);
  glEnd();

  // y-axis
  set_colour(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,1.0,0.0);
  glEnd();

  // z-axis
  set_colour(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,0.0,1.0);
  glEnd();
}

/*********************************************************
    PROC: draw_ground( x, y, z );
    DOES: draws a large square plane orthogonal to the
					y axis to represent the ground, centered at
					the designated point.
**********************************************************/
void draw_ground( float x, float y, float z )
{
  set_colour( 0.1, 0.5, 0.25 );
  glPushMatrix();

  glTranslatef( x, y, z );
  glScalef( 12, 0.2, 12 );
  drawCube();

  glPopMatrix();
}

/******************************************************************
 * ALL additional functions should be below here
 *******************************************************************/

vector<Object> obstacles;

////////////////////////////////////////////////////
// Classes 

#define FABS(A, B) (A)<0?-1*(A):(A)

GLfloat pt_norm (pt s, pt e) {
  return sqrt (powf (FABS (s.x - e.x), 2)
	       + powf (FABS (s.y - e.y), 2));
}

GLfloat pt_slope (pt s, pt e) {
  return  (s.y-e.y) / (s.x-e.x);
}

struct pt {
  GLfloat x, y;
};

class Object
{
public:
  virtual bool collides(pinball b);
  virtual GLfloat time_of_collision(pinball);
};

class Wall : public Object
{
public:
  Wall(struct pt st, struct pt end, double height): 
    _st(st), _end(end), _height(height) {}
  
  ~Wall() {}

  void draw() {
    glBegin( GL_POLYGON );
    glVertex3d(  _st.x,       0,  _st.y );
    glVertex3d( _end.x,       0, _end.y );
    glVertex3d(  _st.x, _height,  _st.y );
    glVertex3d( _end.x, _height, _end.y );
    glEnd();
  }
  
private:
  double _height;
  struct pt _st, _end;
};

bool collides(pinball b) {
  GLfloat a, b, v;
  a = pt_slope( b.pos,   _st );
  b = pt_slope( b.pos,  _end );
  v = pt_slope( b.pos, b.dir );
  
  return (a <= v && v <= b) || (b <= v && v <= a);
}

GLfloat Wall::time_of_collision (pinball p) {
  ;
}

class Pin : public Object {
public:
  bool collides(double start[3], double direction[3], double collision[3]) {
    return true;
  }
};

GLfloat Pin::time_of_collision (pinball p) {
  return 3.1415;
}

struct boundary {
  struct pt pos;
  GLfloat dir;
  GLfloat time;
};

class pinball {
  struct pt pos, dir;
  unsigned int lnq;
  GLfloat velocity;
  std::deque<struct boundary> boundary_q;

  void populate_collisions (GLint N);
public:
  void update_location (void);

  pinball(GLfloat x, GLfloat y, GLfloat dir) {
    pos.x = x;
    pos.y = y;
    this->dir = dir;
  }
};

void pinball::update_location () {
  time_t elapsed_time = Time;

  for (;;) {
    if (elapsed_time < boundary_q.front().time)
      { /* translate the pinball */
	GLfloat distance = velocity * elapsed_time;
	pos.x += cos (dir) * distance;
	pos.y += sin (dir) * distance;
	break;
      } 
    else if (elapsed_time > boundary_q.back().time)
      { /* Overshot the collisions; update, repopulate and repeat */
	struct boundary b = boundary_q.back ();
	elapsed_time -= b.time;
	pos = b.pos;
	dir = b.dir;
	
	boundary_q.clear ();

	populate_collisions (lnq);
      }
    else /* search for new trajectory */
      while (!boundary_q.empty() && boundary_q[0].time < elapsed_time) {
	boundary_q.pop_front ();
	this->populate_collisions (1);
	break;
      }
  }
}

void pinball::populate_collisions (GLint N) {
  for (int i=0; i<N; ++i) {
    struct boundary last = boundary_q.back ();
    GLfloat min_time = obstacles[0].point_of_collision (this);
    int k=0;

    for (int j=1; j<obstacles.size(); ++j) {
      GLfloat tmp = obstacles[j].time_of_collision (this);
      if (tmp < min_time) {
	min_time = tmp;
	k=j;
      }
    }

    //boundary_q.push_back (,min_time);
}

// extra credit function... see the myKey function above.
void moveViewVertical(double* eye, int down)
{
	if(eye[0] < 0) down = -down;
	double sinA = sin(down*moveAngle);
	double cosA = cos(down*moveAngle);
	// the angle is 90 degrees minus theta, where 
	// theta = arctan(z/x).
	double angle = pi/2 - atan(eye[2]/eye[0]);
	double wx = cos(angle);
	double wy = 0;
	double wz = -sin(angle);
	// rodriguez's formula
	double rotation[3][3];
	rotation[0][0] = cosA + wx*wx*(1 - cosA);
	rotation[0][1] = -wz*sinA;
	rotation[0][2] = wx*wz*(1-cosA);
	rotation[1][0] = wz*sinA;
	rotation[1][1] = cosA;
	rotation[1][2] = -wx*sinA;
	rotation[2][0] = wx*wz*(1-cosA);
	rotation[2][1] = wx*sinA;
	rotation[2][2] = cosA+wz*wz*(1-cosA);
	// matrix multiple the eye by this
	int row, col;
	double newEye[3] = {0,0,0};
	for(row = 0; row < 3; row++)
	   for(col = 0; col < 3; col++)
	      newEye[row] += rotation[row][col]*eye[col];
	eye[0] = newEye[0];
	eye[1] = newEye[1];
	eye[2] = newEye[2];
}

/******************************************************************
 * ALL additional functions should be above here
 ******************************************************************/


/*********************************************************

    PROC: display()
    DOES: this gets called by the event handler to draw
          the scene 
      
        MAKE YOUR CHANGES AND ADDITIONS HERE

    Add other procedures if you like.

**********************************************************/
void display(void)
{
  glViewport(0, 0, Width, Height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glFrustum(-5,5,-5,5,4,100);
  //glOrtho(-6,6,-6,6,-500,500) ;

  glMatrixMode(GL_MODELVIEW) ;
  glLoadIdentity();

  gluLookAt (eye[X], eye[Y], eye[Z],
             ref[X], ref[Y], ref[Z], 0.0,1.0,0.0);

  /* glClearColor (red, green, blue, alpha)          */
  /* Ignore the meaning of the 'alpha' value for now */
  glClearColor(0.1f,0.1f,0.7f,1.0f);   /* set the background colour */
  /* OK, now clear the screen with the background colour */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();    
  glScalef(10,10,10);
  draw_axes();
  glPopMatrix();

  /************************************************
   * Start your drawing code here
   *************************************************/

  double width = 3; double height = 2;
  Wall w(width, height);
  double cornerA[3] = { -7.5, 1.0, 6.0 };
  double cornerB[3] = { -7.5, 1.0, -6.0 };
  double cornerC[3] = { 4.5, 1.0, -6.0 };
  double cornerD[3] = { 4.5, 1.0, 6.0 };
  w.draw(cornerA, cornerB);
  w.draw(cornerB, cornerC);
  w.draw(cornerC, cornerD);
  w.draw(cornerD, cornerA);

  /************************************************
   * End your drawing code here
   *************************************************/

  draw_ground( 0.0, 0, 0.0 );

  glutSwapBuffers();
}

/**********************************************
    PROC: myReshape()
    DOES: handles the window being resized 
**********************************************************/

void myReshape(int w, int h)
{
  Width = w;
  Height = h;
}

void idleFunc(void)
{
  Time = TM.GetElapsedTime();
  glutPostRedisplay();
}

/*********************************************************
     PROC: main()
     DOES: calls initialization, then hands over control
	   to the event handler, which calls 
	   display() whenever the screen needs to be redrawn
**********************************************************/

int main(int argc, char** argv) 
{
  
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition (0, 0);
  glutInitWindowSize(Width,Height);
  glutInit(&argc, argv);
  glutCreateWindow(argv[0]);

  myinit();

  glutIdleFunc(idleFunc);
  glutReshapeFunc (myReshape);
  glutKeyboardFunc( myKey );

  glutDisplayFunc(display);
  glutMainLoop();
  
  TM.Reset();
  return 0;         // never reached
}
