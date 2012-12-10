#ifndef __LAYER_H__
#define __LAYER_H__

#include <list>
//#include <GL/glut.h>
#include "Matrix.h"'
#include <QtOpenGL/QGLShaderProgram>
#define LayerPtr ShapeLayer*

class ShapeLayer{
	Vec3 _p0;
public:
	virtual bool isFlat(){return true;};
	virtual void drawGL(){};
};

class LayerManager{
	list<LayerPtr> layers;
public:

	void insert(LayerPtr lp){layers.push_back(lp);};
	void remote(LayerPtr lp){layers.remove(lp);};
	void sendFront(LayerPtr lp){
		layers.remove(lp);
		layers.push_front(lp);
	};

	void sendBack(LayerPtr lp){
		layers.remove(lp);
		layers.push_front(lp);
	};

	void draw(){
		/*int flat =0;
		for(list<LayerPtr>::iterator it = layers.begin(); it!=layers.end(); it++)
			if ((*it)->isFlat())
				flat++;*/
		double z_off = 0;
		glPushMatrix();
		for(list<LayerPtr>::iterator it = layers.begin(); it!=layers.end(); it++){
			glTranslatef(0, 0, z_off);
			(*it)->drawGL();
			z_off+= 0.01;
		}
		glPopMatrix();
	};
};

#endif