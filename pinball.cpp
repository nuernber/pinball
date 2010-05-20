
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
#include <queue>
#include "Timer.h"

Timer TM;

int Width = 480 ;
int Height = 480 ;
float Time = 0;

void set_colour(float r, float g, float b) ;

#define pi 3.141592653589793
#define moveAngle (pi / 24)
void moveViewVertical(double eye[], int down);

#define NUMBER_PINS 10
#define PIN_WIDTH 0.3

#define X 0
#define Y 1
#define Z 2

// The eye point and look-at point.
int currentview = 0;
#define NVIEWS 3
double eye[3] = {1,8,0};

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
    
    if(currentview == 0)
    {
        eye[0] = 1;
        eye[1] = 8;
        eye[2] = 0;
    }
    printf("pinball: eye is %.2f %.2f %.2f, at view %d\n", eye[0], eye[1], eye[2], currentview);
    break;
	    
  case 'n':
    normfns[nfp++%2](GL_AUTO_NORMAL);
    printf (nfp%2 ? "one\n" : "two\n");
    break;
    
	// EXTRA CREDIT FUNCTIONALITY:
	// similar to vim, h, j, k, and l move the view to the left, right
	// up and down.
        case 'h':
            printf("pinball: changing view\n");
	    oldCoord = eye[0];
            eye[0] = eye[0] * cos(moveAngle) - eye[2] * sin (moveAngle);
            eye[2] = oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
            printf("pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'j':
            printf("pinball: changing view\n");
	    moveViewVertical(eye, 1);
            printf("pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'k':
            printf("pinball: changing view\n");
	    moveViewVertical(eye, -1);
            printf("pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
        case 'l':
            printf("pinball: changing view\n");
	    oldCoord = eye[0];
            eye[0] = eye[0] * cos(moveAngle) + eye[2] * sin (moveAngle);
            eye[2] = -oldCoord * sin(moveAngle) + eye[2] * cos (moveAngle);
            printf("pinball: eye is %.2f %.2f %.2f\n", eye[0], eye[1], eye[2]);
            break;
	// n and m zoom in and zoom out
	// NOTE: this only works with a perspective projection, not the
	// default orthogonal projection!
	case 'u':
	    printf("pinball: changing view\n");
	    eye[0] /= 1.1;
	    eye[1] /= 1.1;
	    eye[2] /= 1.1;
	    break;
	case 'i':
	    printf("pinball: changing view\n");
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

class Object
{
public:
	// returns the collision coordinates in collision
	virtual bool collides(double start[3], double direction[3], double collision[3]) = 0;
};

class Wall : public Object
{
public:
	Wall(double width, double height): m_width(width), m_height(height) {}
	~Wall() {}
	// draw wall from start coordinates to end coordinates
	void draw(double start[3], double end[3])
	{
		m_start[0] = start[0]; m_start[1] = start[1]; m_start[2] = start[2];
		m_end[0] = end[0]; m_end[1] = end[1]; m_end[2] = end[2];
		double xIncr = start[0] < 0 ? m_width/2 : -m_width/2;
		double zIncr = start[2] < 0 ? m_width/2 : -m_width/2;

		// draw sides
		if(start[0] == end[0])
		{
			// draw the inner wall, i.e. closer to the origin
			set_colour(1.0, 0.5, 0.0);
			glBegin(GL_POLYGON);
				glVertex3d(start[0] + xIncr, start[1] + m_height/2, start[2]);
				glVertex3d(start[0] + xIncr, start[1] - m_height/2, start[2]);
				glVertex3d(end[0] + xIncr, end[1] - m_height/2, end[2]);
				glVertex3d(end[0] + xIncr, end[1] + m_height/2, end[2]);
				glNormal3d(1.0, 0.0, 0.0);
			glEnd();
			// draw the top
			set_colour(0.0, 0.5, 1.0);
			glBegin(GL_POLYGON);
				glVertex3d(end[0] + xIncr, end[1] + m_height/2, end[2]);
				glVertex3d(end[0] - xIncr, end[1] + m_height/2, end[2]);
				glVertex3d(start[0] - xIncr, start[1] + m_height/2, start[2]);
				glVertex3d(start[0] + xIncr, start[1] + m_height/2, start[2]);
				glNormal3d(0.0, 3.0, 0.0);
			glEnd();
			// draw the outer wall, i.e. further to the origin
			set_colour(1.0, 0.5, 0.0);
			glBegin(GL_POLYGON);
				glVertex3d(end[0] - xIncr, end[1] + m_height/2, end[2]);
				glVertex3d(end[0] - xIncr, end[1] - m_height/2, end[2]);
				glVertex3d(start[0] - xIncr, start[1] - m_height/2, start[2]);
				glVertex3d(start[0] - xIncr, start[1] + m_height/2, start[2]);
				glNormal3d(-3.0, 0.0, 0.0);
			glEnd();
			// draw the ends
			glBegin(GL_POLYGON);
				glVertex3d(end[0] - xIncr, end[1] + m_height/2, end[2]);
				glVertex3d(end[0] + xIncr, end[1] + m_height/2, end[2]);
				glVertex3d(end[0] + xIncr, end[1] - m_height/2, end[2]);
				glVertex3d(end[0] - xIncr, end[1] - m_height/2, end[2]);
				//glNormal3d(-3.0, 0.0, 0.0);
			glEnd();
			glBegin(GL_POLYGON);
				glVertex3d(start[0] - xIncr, start[1] + m_height/2, start[2]);
				glVertex3d(start[0] + xIncr, start[1] + m_height/2, start[2]);
				glVertex3d(start[0] + xIncr, start[1] - m_height/2, start[2]);
				glVertex3d(start[0] - xIncr, start[1] - m_height/2, start[2]);
				//glNormal3d(-3.0, 0.0, 0.0);
			glEnd();
			// draw normals
			if(glIsEnabled(GL_AUTO_NORMAL))
			{
				set_colour(0.0, 0.0, 0.0);
				glBegin(GL_LINES);
					glVertex3d(start[0] + xIncr, start[1], (start[2]+end[2])/2);
					glVertex3d(start[0] + xIncr + 1, start[1], (start[2]+end[2])/2);
					glVertex3d(end[0], end[1] + m_height/2, (start[2]+end[2])/2);
					glVertex3d(end[0], end[1] + m_height/2 + 1, (start[2]+end[2])/2);
					glVertex3d(end[0] - xIncr, end[1], (start[2]+end[2])/2);
					glVertex3d(end[0] - xIncr - 1, end[1], (start[2]+end[2])/2);
				glEnd();
			}
		} else {
			// draw the inner wall, i.e. closer to the origin
			set_colour(1.0, 0.5, 0.0);
			glBegin(GL_POLYGON);
				glVertex3d(start[0], start[1] + m_height/2, start[2] + zIncr);
				glVertex3d(start[0], start[1] - m_height/2, start[2] + zIncr);
				glVertex3d(end[0], end[1] - m_height/2, end[2] + zIncr);
				glVertex3d(end[0], end[1] + m_height/2, end[2] + zIncr);
			glEnd();
			// draw the top
			set_colour(0.0, 0.5, 1.0);
			glBegin(GL_POLYGON);
				glVertex3d(end[0], end[1] + m_height/2, end[2] + zIncr);
				glVertex3d(end[0], end[1] + m_height/2, end[2] - zIncr);
				glVertex3d(start[0], start[1] + m_height/2, start[2] - zIncr);
				glVertex3d(start[0], start[1] + m_height/2, start[2] + zIncr);
			glEnd();
			// draw the outer wall, i.e. further to the origin
			set_colour(1.0, 0.5, 0.0);
			glBegin(GL_POLYGON);
				glVertex3d(end[0], end[1] + m_height/2, end[2] - zIncr);
				glVertex3d(end[0], end[1] - m_height/2, end[2] - zIncr);
				glVertex3d(start[0], start[1] - m_height/2, start[2] - zIncr);
				glVertex3d(start[0], start[1] + m_height/2, start[2] - zIncr);
			glEnd();
			// draw the ends
			glBegin(GL_POLYGON);
				glVertex3d(end[0], end[1] + m_height/2, end[2] - zIncr);
				glVertex3d(end[0], end[1] + m_height/2, end[2] + zIncr);
				glVertex3d(end[0], end[1] - m_height/2, end[2] + zIncr);
				glVertex3d(end[0], end[1] - m_height/2, end[2] - zIncr);
				//glNormal3d(-3.0, 0.0, 0.0);
			glEnd();
			glBegin(GL_POLYGON);
				glVertex3d(start[0], start[1] - m_height/2, start[2] - zIncr);
				glVertex3d(start[0], start[1] - m_height/2, start[2] + zIncr);
				glVertex3d(start[0], start[1] + m_height/2, start[2] + zIncr);
				glVertex3d(start[0], start[1] + m_height/2, start[2] - zIncr);
				//glNormal3d(-3.0, 0.0, 0.0);
			glEnd();
			// draw normals
			if(glIsEnabled(GL_AUTO_NORMAL))
			{
				set_colour(0.0, 0.0, 0.0);
				glBegin(GL_LINES);
					glVertex3d((start[0]+end[0])/2, start[1], start[2] + zIncr);
					glVertex3d((start[0]+end[0])/2, start[1], start[2] + zIncr + 1);
					glVertex3d((start[0]+end[0])/2, end[1] + m_height/2, end[2]);
					glVertex3d((start[0]+end[0])/2, end[1] + m_height/2 + 1, end[2]);
					glVertex3d((start[0]+end[0])/2, end[1], end[2] - zIncr);
					glVertex3d((start[0]+end[0])/2, end[1], end[2] - zIncr - 1);
				glEnd();
			}
		}
	}

	// returns the collision coordinates in collision
	bool collides(double start[3], double direction[3], double collision[3])
	{
		return true;
	}

private:
	double m_width;
	double m_height;
	double m_start[3];
	double m_end[3];
};

class Pin : public Object
{
public:
	void draw(double coords[3], double height)
	{
		int nRect = 24;
		double x1, z1, x2, z2;
		for(int i = 0; i < nRect; i++)
		{
			double scale = PIN_WIDTH;
			x1 = sin(2*pi/nRect*i) * scale + coords[0];
			z1 = cos(2*pi/nRect*i) * scale + coords[2];
			x2 = sin(2*pi/nRect*(i+1)) * scale + coords[0];
			z2 = cos(2*pi/nRect*(i+1)) * scale + coords[2];
			set_colour(1.0, 0.8, 0.3);
			glBegin(GL_POLYGON);
				glVertex3d(x1, 0, z1);
				glVertex3d(x2, 0, z2);
				glVertex3d(x2, height, z2);
				glVertex3d(x1, height, z1);
			glEnd();
			// draw normals for sides and top
			if(glIsEnabled(GL_AUTO_NORMAL))
			{
				set_colour(0, 0, 0);
				glBegin(GL_LINES);
					glVertex3d((x1+x2)/2, height/2, (z1+z2)/2);
					double newX = (x1+x2)/2 + ((x1+x2)/2 - coords[0]);
					double newZ = (z1+z2)/2 + ((z1+z2)/2 - coords[2]);
					glVertex3d(newX, height/2, newZ);
				glEnd();
			}
		}
		// draw sphere for top
		set_colour(1.0, 0.8, 0.3);
		glPushMatrix();
			glTranslated(coords[0], coords[1] + height, coords[2]);
			glScalef(PIN_WIDTH, PIN_WIDTH, PIN_WIDTH);
			drawSphere();
		glPopMatrix();
		// draw normal on top of sphere
		if(glIsEnabled(GL_AUTO_NORMAL))
		{
			glBegin(GL_LINES);
				glVertex3d(coords[0], height + PIN_WIDTH, coords[2]);
				glVertex3d(coords[0], height + PIN_WIDTH + 1, coords[2]);
			glEnd();
		}

		m_coords[0] = coords[0];
		m_coords[1] = coords[1];
		m_coords[2] = coords[2];
	}

	void getCoords(double coords[3])
	{
		coords[0] = m_coords[0];
		coords[1] = m_coords[1];
		coords[2] = m_coords[2];
	}

	// returns the collision coordinates in collision
	bool collides(double start[3], double direction[3], double collision[3])
	{
		return true;
	}
private:
	double m_coords[3];
};


struct pt {
  GLfloat x, y;
};

struct boundary {
  struct pt pos, direction;
  GLfloat time;
};

class pinball {
public:
  struct pt pos, direction;

  void location (GLfloat elapsed_time);

};


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

  if(currentview == 0)
      glOrtho(-6,6,-6,6,-500,500);
  else if(currentview == 1)
  {
      // flythrough
      glFrustum(-5,5,-5,5,4,100);
      int a = 16;
      int b = 8;
      const int period = 10; // seconds
      const int yPeriod = 5; // seconds
      eye[X] = a * cos(Time * (2*pi) / period);
      eye[Y] = 7.5 + sin(Time * (2*pi) / yPeriod);
      eye[Z] = b * sin(Time * (2*pi) / period);
  }
  else
      glFrustum(-5,5,-5,5,4,100);

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
  double width = 0.5; double height = 2;
  Wall w1(width, height);
  Wall w2(width, height);
  Wall w3(width, height);
  Wall w4(width, height);
  double cornerA[3] = { -6.0, 1.0, 5.5 };
  double cornerB[3] = { -6.0, 1.0, -5.5 };
  double cornerC[3] = { 6.0, 1.0, -5.5 };
  double cornerD[3] = { 6.0, 1.0, 5.5 };
  double cornerA2[3] = { -6.25, 1.0, 5.75 };
  double cornerB2[3] = { -6.25, 1.0, -5.75 };
  double cornerC2[3] = { 6.25, 1.0, -5.75 };
  double cornerD2[3] = { 6.25, 1.0, 5.75 };
  w1.draw(cornerA, cornerB);
  w2.draw(cornerB2, cornerC2);
  w3.draw(cornerC, cornerD);
  w4.draw(cornerD2, cornerA2);

  Pin p[NUMBER_PINS];
  double pinCoords[NUMBER_PINS][2] = {
    {5.5, 3.5}, {3, 0}, {4, -4}, {1, -3}, {3.75, -1},
    {-1.25, 3}, {-3.5, 5}, {-0.5, 1.8}, {-1, -3.2}, {-4, -2.7}
  };
  for(int i = 0; i < NUMBER_PINS; i++)
  {
  	double x = pinCoords[i][0];
	double z = pinCoords[i][1];
  	double coords[3] = {x, 0.0, z};
	double pHeight = 3.0;
	p[i].draw(coords, pHeight);
  }
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
