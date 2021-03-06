 
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

void drawCylinder(void);
void drawCone(void);
void drawCube(void);
void drawSphere(void);

#define pi 3.141592653589793
#define moveAngle (pi / 24)
void moveViewVertical(double eye[], int down);

#define NUMBER_PINS 10
#define PIN_WIDTH 0.2
#define PINBALL_RADIUS 0.4

#define X 0
#define Y 1
#define Z 2

// EXTRA CREDIT: more flythroughs
int currentfly = 0;
#define NFLYS 3

// The eye point and look-at point.
int currentview = 0;
#define NVIEWS 3
double eye[3] = {1,8,0};

double ref[3] = {0.0, 0.0, 0.0};

void (*normfns[])(GLenum) = {glEnable, glDisable};
bool wireframe = false;
int nfp;

#define FABS(A) (A)<0?-1*(A):(A)


struct point2D {
  GLdouble x, y;
};

GLdouble pt_norm (struct point2D s, struct point2D e) {
  return sqrt (powf (FABS (s.x - e.x), 2)
	       + powf (FABS (s.y - e.y), 2));
}

GLdouble pt_slope (struct point2D s, struct point2D e) {
  return  (s.y-e.y) / (s.x-e.x);
}

class Pin {
public:
  GLdouble radius, height;
  struct point2D center;
  
  void draw (void);

  Pin (struct point2D c) {
    radius = PIN_WIDTH;
    height = 1.5;
    center.x = c.x;
    center.y = c.y;
  }
};

class Pinball {
public:
  struct point2D pos;
  struct point2D v;
  GLdouble radius;
  GLdouble velocity;
  GLdouble time2;
  GLdouble dt;
  bool robot_takeover;

  void update (GLdouble time2);
  void draw (void);
  void random ();
  void bounce (Pin p, GLdouble t);
  void increaseVelocity() { velocity++; }
  void decreaseVelocity() { velocity--; }
  void transmogrify (void) {robot_takeover = !robot_takeover;}
  
  void init (GLdouble x, GLdouble y, GLdouble theta) {
    pos.x = x;
    pos.y = y;
    v.x = cos (theta);
    v.y = sin (theta);
    radius = PINBALL_RADIUS;
    velocity = 5;
    robot_takeover = false;
  }
};

void Pinball::random () {
  srand(time(NULL));
  GLdouble w = 6 - PINBALL_RADIUS;
  GLdouble x = (rand() / (double)RAND_MAX) * 2 * w - w;
  GLdouble y = (rand() / (double)RAND_MAX) * 2 * w - w;
  GLdouble theta = (rand() / (double)RAND_MAX) * 360;
  init (x, y, theta);
}


void full_cylinder (GLfloat w, GLfloat l) {
  glPushMatrix ();

  glPushMatrix (); /* bottom of the cylinder */
  glScalef     (w - 1E-5, w - 1E-5, .02);
  drawCone     ();
  glPopMatrix  ();

  glPushMatrix (); 
  glScalef     (w, w, l);
  drawCylinder ();
  glPopMatrix  (); 
  
  glTranslatef (0, 0, l);
  
  glPushMatrix ();  /* top of the cylinder */
  glScalef     (w - 1E-5, w - 1E-5, .02);
  glRotatef    (180, 1, 0, 0); /* upside down */
  drawCone     ();
  glPopMatrix  (); 
  
  glPopMatrix  ();
}

void stiff_arm () {
  GLfloat shoulderw = 1;
  GLfloat shoulderlen = 1;
  GLfloat armw = .25;
  GLfloat armlen = 3.5;
  GLfloat handw = 1.5;
  GLfloat fingerw = .5;

  glPushMatrix ();
  
  glPushMatrix ();
  full_cylinder (shoulderw, shoulderlen);
  glPopMatrix  ();

  glTranslatef (shoulderw, 0, (shoulderlen - (armw/2))/2 );
  
  glPushMatrix ();
  glRotatef    (90, 0, 1, 0);
  glScalef     (armw, armw, armlen);
  drawCylinder ();
  glPopMatrix  ();
  
  glTranslatef (armlen, 0, 0);

  glPushMatrix ();
  glScalef     (fingerw, handw, fingerw);
  drawCube     ();
  glPopMatrix  ();

  glPushMatrix ();
  glTranslatef (fingerw, handw/2 - fingerw/2, 0);
  glScalef     (handw/2, fingerw, fingerw);
  drawCube     ();
  glPopMatrix  ();

  glPushMatrix ();
  glTranslatef (fingerw, -(handw/2 - fingerw/2), 0);
  glScalef     (handw/2, fingerw, fingerw);
  drawCube     ();
  glPopMatrix  ();

  glPopMatrix  ();
}

/* His name is Bob */
void bob (bool odd) {
  GLfloat legw = 0.25;
  GLfloat leglen = 2;
  GLfloat torsolen = 8;
  GLfloat torsow = 2;
  GLfloat headlen = 3;
  GLfloat eye_seperation = 1.3;
  GLfloat eye_height = 2;
  GLfloat antennaw = 0.10;
  GLfloat antennalen = 1.5;
  GLfloat antenna = 0.35;

  glPushMatrix ();
  
  glPushMatrix (); // wheel
  set_colour   (0, 0, .5); // black
  drawSphere   ();
  glPopMatrix  (); // wheel

  glRotatef    (-90, 1, 0, 0); // z vertical

  glTranslatef (0, 0, sqrt (1 - legw*legw));

  glPushMatrix (); // leg
  set_colour   (.5, .5, 0); // yellow
  glScalef     (legw, legw, leglen);
  drawCylinder ();
  glPopMatrix  ();

  glTranslatef (0, 0, leglen);

  set_colour   (0, .5, 0);
  full_cylinder (torsow, torsolen);

  glTranslatef (0, 0, torsolen);

  glPushMatrix (); // save matrix for later
  
  glPushMatrix ();
  set_colour   (.2, .2, .2);
  glTranslatef (0, 0, headlen/2);
  glScalef     (headlen, headlen, headlen);
  drawCube     ();
  glPopMatrix  ();

  glTranslatef (0, 0, headlen);
  
  glPushMatrix ();
  set_colour   (.7, .3, .2);
  glScalef     (antennaw, antennaw, antennalen);
  drawCylinder ();
  glPopMatrix  ();

  glTranslatef (0, 0, antennalen);

  glPushMatrix ();
  set_colour   (.5, 0, 0);
  glScalef     (antenna, antenna, antenna);
  drawSphere   ();
  glPopMatrix  ();

  glPopMatrix  (); // return to torso height

  glPushMatrix (); // eyes

  set_colour   (.5, .5, .5);
  glTranslatef (eye_seperation/2, -headlen/2, eye_height);
  glRotatef    (90, 1, 0, 0);
  full_cylinder (.25, .25);
  glPopMatrix  ();

  glPushMatrix (); // eyes
  glTranslatef (-eye_seperation/2, -headlen/2, eye_height);
  glRotatef    (90, 1, 0, 0);
  full_cylinder (.25, .25);
  glPopMatrix  ();

  glPushMatrix (); // right arm
  set_colour   (.5, 0, .5);
  glTranslatef (torsow, 0, -1);
  glRotatef    (90, 0, 1, 0);
  if (odd)
    glRotatef    (-90, 0, 0, 1);
  stiff_arm    ();
  glPopMatrix  ();

  glPushMatrix (); // left arm
  glTranslatef (-torsow, 0, -1);
  glRotatef    (-90, 0, 1, 0);
  glRotatef    (180, 0, 0, 1);
  stiff_arm    ();
  glPopMatrix  ();
    
  glPopMatrix  ();
}

void Pinball::draw() {
  glPushMatrix ();

  set_colour (1, 0, 0);
  glTranslatef (pos.x, radius, pos.y);

  if (!robot_takeover) {
    glScalef (radius, radius, radius);
    drawSphere ();
  } else {
    glScalef (radius/2., radius/2., radius/2.);
    bob (!!((int)round(Time*5)%2));
  }

  glPopMatrix ();
}

class Wall {
public:
  GLdouble height, width;
  struct point2D inside[4], outside[4];

  void face (struct point2D st, struct point2D end);

  void draw (void);
  
  struct collision *next_collision (Pinball p);

  Wall (double w, double h, struct point2D walls[4]) {
    width = w;
    height = h;
    GLdouble dx, dy;
    dx = dy = width;

    for (int i=0; i<4; i++) {
      inside[i].x  = walls[i].x;
      inside[i].y  = walls[i].y;
      outside[i].x = walls[i].x+ dx;
      outside[i].y = walls[i].y+ dy;
      if (i%2) dx = -dx;
      else dy = -dy;
    }
  }
};


void Wall::face (struct point2D st, struct point2D end) {

  glBegin (GL_POLYGON);
  set_colour (1.0, 0.5, 0.0);
  glVertex3d (  st.x,      0,  st.y );
  glVertex3d ( end.x,      0, end.y );
  glVertex3d ( end.x, height, end.y );
  glVertex3d (  st.x, height,  st.y );
  glNormal3d (0, 0, 1);
  glEnd ();
  if (glIsEnabled (GL_AUTO_NORMAL)) {
    glBegin(GL_LINES);
    if (st.x == end.x ) // since everything is parallel
    {
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2);
      glVertex3d((st.x + end.x)/2 - 1, height/2, (st.y+end.y)/2);
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2);
      glVertex3d((st.x + end.x)/2 + 1, height/2, (st.y+end.y)/2);
      glNormal3d (0, 0, 1);
    }
    else
    {
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2);
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2 - 1);
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2);
      glVertex3d((st.x + end.x)/2, height/2, (st.y+end.y)/2 + 1);
      glNormal3d (0, 0, 1);
    }
    glEnd();
  }
}

void Wall::draw() {
  int i, j;
  for (i=0;i<4;i++) {
    j = (i+1)%4;
    this->face (  inside[j],  inside[i] );
    this->face ( outside[i], outside[j] );

    glBegin (GL_POLYGON);
    set_colour(0.0, 0.5, 1.0);
    glVertex3d ( outside[i].x, height, outside[i].y );
    glVertex3d ( outside[j].x, height, outside[j].y );
    glVertex3d (  inside[j].x, height, inside[j].y  );
    glVertex3d (  inside[i].x, height, inside[i].y  );
    glNormal3d ( 0, 0, 1);
    glEnd ();

  }
}


void Pin::draw() {
  set_colour (0, 0, 1);
  glPushMatrix ();
  glTranslated (center.x, 0, center.y);

  glPushMatrix ();
  glRotatef (-90, 1, 0, 0);
  glScaled (radius, radius, height);
  if (glIsEnabled (GL_AUTO_NORMAL)) {
    glBegin(GL_LINES);
      for(int i = 0; i < 24; i++)
      {
        glVertex3d(radius*cos(i*2*pi/24), radius*sin(i*2*pi/24), height/2);
        glVertex3d((1+radius)*cos(i*2*pi/24)/radius, (1+radius)*sin(i*2*pi/24)/radius, height/2);
      }
    glEnd();
  }
  drawCylinder ();
  glPopMatrix ();

  glTranslated (0, height, 0);
  glScaled (radius, radius, radius);
  drawSphere ();
  glPopMatrix ();
}

GLdouble width = 0.5, height = 2;
struct point2D I   = { 6,  6};
struct point2D II  = { 6, -6};
struct point2D III = {-6, -6};
struct point2D IV  = {-6,  6};
 
struct point2D corners[] = {I, II, III, IV}; 
Wall wall (width, height, corners);
std::vector<Pin> pins;
Pinball pinball;

void Pinball::bounce (Pin p, GLdouble t) {
  struct point2D incident = v;
  struct point2D normal;
  normal.x = pos.x - p.center.x;
  normal.y = pos.y - p.center.y;
  normal.x /= pt_norm(pos, p.center);
  normal.y /= pt_norm(pos, p.center);
  struct point2D newV;
  // dot product of normal with incident
  double dp = incident.x*normal.x + incident.y*normal.y;
  newV.x = -2*dp*normal.x + incident.x;
  newV.y = -2*dp*normal.y + incident.y;
  v.x = newV.x;
  v.y = newV.y;
  double norm = sqrt(v.x*v.x+v.y*v.y);
  v.x /= norm;
  v.y /= norm;
}

void Pinball::update (GLdouble t) {
  GLdouble oldx = pos.x;
  GLdouble oldy = pos.y;
  pos.x += velocity * v.x * (t-time2);
  pos.y += velocity * v.y * (t-time2);

  if (pos.x + radius >= 6)
  {
    pos.x = oldx;
    pos.y = oldy;
    v.x = -1 * v.x;
  }
  if (pos.x - radius <= -6)
  {
    pos.x = oldx;
    pos.y = oldy;
    v.x = -1 * v.x;
  }
  if (pos.y + radius >= 6)
  {
    pos.x = oldx;
    pos.y = oldy;
    v.y = -1 * v.y;
  }
  if (pos.y - radius <= -6)
  {
    pos.x = oldx;
    pos.y = oldy;
    v.y = -1 * v.y;
  }

  for (int i=0;i<NUMBER_PINS; i++)
    if (pt_norm (pins[i].center, pos) <= pins[i].radius + radius) {
      pos.x = oldx;
      pos.y = oldy;
      this->bounce (pins[i], t);
      break;
      return;
    }

  time2 = t;
}


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
  // EXTRA CREDIT: increase and decrease pinball velocity
  case 'p':
    pinball.increaseVelocity();
    break;
  case 'o':
    pinball.decreaseVelocity();
    break;
  // EXTRA CREDIT: more flythrough
  case 'f':
    currentfly = (currentfly + 1 ) % NFLYS;
    break;
  case 'r':
    pinball.random ();
    break;

  // EXTRA CREDIT: robot instead of a pinball
  case 'x':
    pinball.transmogrify ();
    break;

  case 's':
    if (wireframe)
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    else
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    wireframe = !wireframe;
    break;
    
  case 'q':
    exit(0); 
    break ;

  case 'v':
    printf("Pinball: changing view\n");
    currentview = (currentview + 1) % NVIEWS;
    
    if(currentview == 0)
      {
        eye[0] = 1;
        eye[1] = 8;
        eye[2] = 0;
      }
    printf("Pinball: eye is %.2f %.2f %.2f, at view %d\n", eye[0], eye[1], eye[2], currentview);
    break;

  case 'n':
    normfns[nfp++%2](GL_AUTO_NORMAL);
    break;
    
    // EXTRA CREDIT FUNCTIONALITY:
    // similar to vim, h, j, k, and l move the view to the left, right
    // up and down.
  case 'h':
    printf("Pinball: changing view\n");
    oldCoord = eye[0];
    eye[0] = eye[0] * cos(moveAngle) - eye[2] * sin (moveAngle);
    eye[2] = oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
    printf("Pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
    break;
  case 'j':
    printf("Pinball: changing view\n");
    moveViewVertical(eye, 1);
    printf("Pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
    break;
  case 'k':
    printf("Pinball: changing view\n");
    moveViewVertical(eye, -1);
    printf("Pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
    break;
  case 'l':
    printf("Pinball: changing view\n");
    oldCoord = eye[0];
    eye[0] = eye[0] * cos(moveAngle) + eye[2] * sin (moveAngle);
    eye[2] = -oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
    printf("Pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
    break;

    // n and m zoom in and zoom out
    // NOTE: this only works with a perspective projection, not the
    // default orthogonal projection!
  case 'u':
    printf("Pinball: changing view\n");
    eye[0] /= 1.1;
    eye[1] /= 1.1;
    eye[2] /= 1.1;
    break;
  case 'i':
    printf("Pinball: changing view\n");
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
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
	
  glPixelStorei(GL_PACK_ALIGNMENT,1) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT,1) ;
  glShadeModel(GL_SMOOTH) ;

  /* obstacles */
  double pinCoords[][2] = {
    {5.5, 3.5}, {3, 0}, {4, -4}, {1, -3}, {3.75, -1},
    {-1.25, 3}, {-3.5, 5}, {-0.5, 1.8}, {-1, -3.2}, {-4, -2.7}
  };

  struct point2D point;
  for(int i = 0; i < NUMBER_PINS; i++) {
    point.x = pinCoords[i][0];
    point.y = pinCoords[i][1];
    Pin pin (point);
    pins.push_back (pin);
  }

  pinball.random ();
  
  srand (time (NULL));
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

////////////////////////////////////////////////////
// Classes 

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

  if(currentview == 0) {
    glOrtho(-8,8,-8,8,-500,500);
    ref[X] = ref[Y] = ref[Z] = 0;
  }
  else if(currentview == 1)
    {
      // flythrough
      glFrustum(-5,5,-5,5,4,100);
      int a = 16;
      int b = 8;
      int period = 10; // seconds
      int yPeriod = 5; // seconds
      // EXTRA CREDIT: different flythroughs
      switch(currentfly)
      {
      case 0:
        a = 16; b = 8; period = 10; yPeriod = 5;
        break;
      case 1:
        a = 20; b = 4; period = 5; yPeriod = 2.5;
        break;
      case 2:
        a = 10; b = 10; period = 20; yPeriod = 10;
        break;
      }
      eye[X] = a * cos(Time * (2*pi) / period);
      eye[Y] = 7.5 + sin(Time * (2*pi) / yPeriod);
      eye[Z] = b * sin(Time * (2*pi) / period);
    }
  else
  {
    //glFrustum(-5,5,-5,5,4,100);
    gluPerspective (60, 1, .1, 100);
    ref[X] = pinball.pos.x;
    ref[Y] = pinball.radius*2.25;
    ref[Z] = pinball.pos.y;

    eye[X] = pinball.pos.x - pinball.v.x;
    eye[Y] = pinball.radius*2.25;
    eye[Z] = pinball.pos.y- pinball.v.y;
  }
  
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

  wall.draw ();
  for (int i=0; i<NUMBER_PINS; i++)
    pins[i].draw ();

  pinball.draw ();

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
  pinball.update (Time);
  glutPostRedisplay();
}


/*********************************************************
     PROC: main()
     DOES: calls initialization, then hands over control
	   to the event handler, which calls 
	   display() whenever the screen needs to be redrawn
**********************************************************/

void reshape (int width, int height)
{
  if (height == 0) height =1;
  glViewport (0,0,width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

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
  glutReshapeFunc (reshape);
  glutDisplayFunc(display);
  glutMainLoop();
  
  TM.Reset();
  return 0;         // never reached
}
