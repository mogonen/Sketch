#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <time.h>
#include <string.h>
#include <GL/glut.h>

#include "Vector.h"
#include "Camera.h"
#include "Stroke.h"
#include "Curve2Mesh.h"
using namespace std;

Vec3 (*ControlPoint::toView)(Vec3 p);
int WIDTH = 1200;
int HEIGHT = 900;
int persp_win;

bool SHOW_CONTROL_POINTS = true;
bool SHOW_MESH = true;
Camera* camera;
bool showGrid = true;
// draws a simple grid
StrokeManager * strokeman;
CurveBoss * boss;
LayerManager * layerman;

void makeGrid() {
  glColor3f(0.93, 0.93, 0.93);
  glLineWidth(1.0);

  for (float i=-2; i<2; i+=0.1) {
    for (float j=-2; j<2; j+=0.1) {
      glBegin(GL_LINES);
      glVertex3f(i, 0, j);
      glVertex3f(i, 0, j+1);
      glEnd();
      glBegin(GL_LINES);
      glVertex3f(i, 0, j);
      glVertex3f(i+1, 0, j);
      glEnd();

      if (j == 1.9){
	glBegin(GL_LINES);
	glVertex3f(i, 0, j+1);
	glVertex3f(i+1, 0, j+1);
	glEnd();
      }
      if (i == 1.9){
	glBegin(GL_LINES);
	glVertex3f(i+1, 0, j);
	glVertex3f(i+1, 0, j+1);
	glEnd();
      }
    }
  }

  glLineWidth(2.0);
  glBegin(GL_LINES);
  glVertex3f(-2, 0, 0);
  glVertex3f(2, 0, 0);
  glEnd();
  glBegin(GL_LINES);
  glVertex3f(0, 0, -2);
  glVertex3f(0, 0, 2);
  glEnd();
  glLineWidth(1.0);
}

Vec3 toView(Vec3 p){
	return camera->toView(p);
}

void init() {

  // set up camera
  camera = new Camera(Eye::get()->P, Eye::get()->N, Eye::get()->U, PZ, 1000, FOV);
  //glClearColor(0.98, 0.98, 0.98, 0.00);
  glClearColor(0.781, 0.781, 0, 0.00);
  glShadeModel(GL_SMOOTH);
  glDepthRange(0.0, 1.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);

  strokeman = StrokeManager::getManager();
  boss = CurveBoss::getBoss();
  ControlPoint::toView = toView;
}

Vec3 genEllipse(Vec3 p0, Vec3 * p, int num, double t){

	double ta = t*2*PI;
	Vec3 n0(p[0] - p0);
	n0 = n0.normalize();
	Vec3 ny = n0 % (p[1]-p0).normalize();
	Vec3 n1 = (ny%n0).normalize();
	Vec3 n = n0*cos(ta) + n1*sin(ta);;

	double base_ang = 0;
	int sec = 0;
	for(int i=0; i<num; i++){
		Vec3 ni(p[i] - p0);
		ni = ni.normalize();
		double d = n0 * ni;
		double ba = acos(d);
		ba = ( n1 * ni >=0)? ba: (2*PI - ba);
		if (ta > ba){
			base_ang = ba;
			sec = i;
		}
	}
	double aa = (ta - base_ang);
	double ang = acos( (p[sec]-p0).normalize() * (p[(sec+1)%num]-p0).normalize() );
	double s = aa  / ang;
	int sec1 = (sec+1)%num;

	Vec3 t0 = ( ny % (p[sec] - p0).normalize()).normalize()*0.5;
	Vec3 t1 = ( ny % (p[sec1]- p0).normalize()).normalize()*0.5;

	float h1 =  2*s*s*s - 3*s*s + 1; // calculate basis function 1
	float h2 = -2*s*s*s + 3*s*s;   // calculate basis function 2
	float h3 =  s*s*s - 2*s*s + s; // calculate basis function 3
	float h4 =  s*s*s -  s*s;   
	//cout<<ta<<":"<<aa<<" "<<sec<<endl;
	//double rad =  (p[sec]-p0).norm()*(1-tt) + (p[(sec+1)%num]-p0).norm()*tt;

	return p[sec]*h1 + p[sec1]*h2 + t0*h3 + t1*h4;
}

void drawGLScene(){

	if (SHOW_CONTROL_POINTS)
		ControlPoint::drawAll();

	strokeman->drawGL();
	if (boss)
		boss->drawGL();
/*	int N = 2;
	Vec3 * p = new Vec3[N];
	p[0] = Vec3(-0.5, 0, -PZ);
	p[1] = Vec3(0, 0.25, -PZ);
	/*p[2] = Vec3(0.3, 0.1, -PZ);
	/*p[3] = Vec3(0.2, -0.5, -PZ);
	p[4] = Vec3(-0.3, -0.4, -PZ);
	
	Vec3 p0 = Vec3(0, 0, -PZ);
	glBegin(GL_LINE_LOOP);
	for(double t=0; t<0.95; t+=0.01){
		Vec3 pp = genEllipse(p0, p, N, t);
		glVertex3f(pp.x, pp.y, pp.z);
	}
	glEnd();

	glBegin(GL_POINTS);
	for(int i=0; i<N; i++)
		glVertex3f(p[i].x, p[i].y, p[i].z);
	glEnd();*/
}

void PerspDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // draw the camera created in perspective
  camera->PerspectiveDisplay(WIDTH, HEIGHT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (showGrid) 
    makeGrid();

  drawGLScene();
  glutSwapBuffers();
}

void mouseEventHandler(int button, int state, int x, int y) {
  int mode = glutGetModifiers();
  // let the camera handle some specific mouse events (similar to maya)
  Vec3 mp = camera->HandleMouseEvent(button, state, x, y);
  if ( mode != GLUT_ACTIVE_ALT && mode != GLUT_ACTIVE_CTRL){
	  strokeman->HandleMouseEvent(button, state, mp);
	  /*if (state == GLUT_UP)
		  boss->step();*/
  }
  glutPostRedisplay();
}

void motionEventHandler(int x, int y) {
  // let the camera handle some mouse motions if the camera is to be moved
  Vec3 mp = camera->HandleMouseMotion(x, y);
  strokeman->HandleMouseMotion(mp);
  if (camera->isInactive())
	ControlPoint::handleMouseAll(0,mp);
  glutSetWindow(persp_win);
  glutPostRedisplay();
}

void keyboardEventHandler(unsigned char key, int x, int y) {
  switch (key) {
  case '1':
    // reset the camera to its initial position
		camera->Reset();
		Eye::setCam(camera);
    break;
  case '2':
	  camera->setAll(TOP_P, AIM, TOP_U);
		Eye::setCam(camera);
  break;
  case '3':
	  camera->setAll(SOL_P, AIM, SOL_U);
		Eye::setCam(camera);
  break;
  case 'f':
    //camera->SetCenterOfFocus(Vec3(0, 0, 0));
    break;
  case 'g':
    showGrid = !showGrid;
    break; 
  case 'c':
	  SHOW_CONTROL_POINTS=!SHOW_CONTROL_POINTS;
  break;

  case 13:
	  boss->reset();
  break;

  case ' ':
	  boss->step();
    break;
  case '+':
	  //c2s->buildCorrespondance();
	  break;
  case 'p':
	  /*if(c2s)
		  c2s->rebuildRingsByPerspective();*/
	  //project
	  camera->SetClippingPlanes(PZ,1000);
	  Eye::setCam(camera);
	  boss->reproject(camera);
  break;

   case 'a':
		::DRAW_AUX = !::DRAW_AUX;
	   cout<<"aux:"<<DRAW_AUX<<endl;
	break;

   case 's':
	   boss->saveMesh("e:/temp/out.obj");
   break;

   case 'i':
	   strokeman->insertStrokes(0);
	break;

   case 'x':
	   strokeman->deleteLast();
   break;
   case 'n':
	   boss->reset();
	   strokeman->reset();
	   ControlPoint::clear();
	   camera->Reset();
	   Eye::setCam(camera);
   break;
   case 'm':
	  if (boss)
		  boss->buildMesh();
	  ::DRAW_AUX = true;
   break;
   case ',':
	   SHOW_MESH=!SHOW_MESH;
   break;

   case 'q':
	   camera->SetClippingPlanes(0.1, 1000);
   break;
     case 'w' :
	  camera->SetClippingPlanes(PZ,1000);
  break;

  }

  glutPostRedisplay();
}

static void initGLScene(){

	glShadeModel(GL_SMOOTH);
	glDepthRange(0.0, 1.0);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);

	//glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	GLfloat l0pos[]= { 0.0f, 0.0f, 5.0f, 1.0f };   
	glLightfv(GL_LIGHT0, GL_POSITION, l0pos);
	//glEnable(GL_LIGHT1);

	/*
	glEnable(GL_LIGHT2);
	glEnable(GL_LIGHT3);
	glEnable(GL_LIGHT4);
	glEnable(GL_LIGHT5);
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_BLEND);
	//glEnable(GL_CULL_FACE);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	*/

	//glEnable(GL_TEXTURE_2D);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEX7TURE_ENV_MODE, GL_MODULATE);
}


void initScene(){
}

int main(int argc, char *argv[]) {

  initScene();
  // set up opengl window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowSize(WIDTH, HEIGHT);
  glutInitWindowPosition(50, 50);
  persp_win = glutCreateWindow("Scene");

  // initialize the camera and such
  init();
  initGLScene();
  //updateGLScene();

  // set up opengl callback functions
  glutDisplayFunc(PerspDisplay);
  glutMouseFunc(mouseEventHandler);
  glutMotionFunc(motionEventHandler);
  glutKeyboardFunc(keyboardEventHandler);

  glutSetWindow(persp_win);

  glutMainLoop();
  return(0);
}