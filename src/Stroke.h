#ifndef _H_Stroke
#define _H_Stroke

#include <list>
#include <stdio.h>
#include <stdlib.h>
//#include <GL/glut.h>
#include "Vector.h"
#include "Curve.h"
//#include "Eye.h"

using namespace std;
class Dot;
class Stroke;
class Curve;

#define PZ 1.0
//#define DotPtr Dot*
#define StrokePtr Stroke*

class Stroke:public Curve{

	list<Vec3> _pts;
	CurvePtr _curve;
	CurvePtr _hcurve;

public:

	Stroke();
	Stroke(CurvePtr cp){ _curve = cp; };

	void drawGL();
	void add(const Vec3& p);
	void finalize();

	double length();
	int size();
	const Vec3& getP(double t);
	Vec3* toArr(int len = 0);

	CurvePtr curve(){return _curve;};
	void makeNURBS();
};

class StrokeManager{

	list<StrokePtr> strokes;
	StrokePtr active;
	Vec3 m_prev;
	StrokeManager();
	static StrokeManager * _manager;

public: 

	void add(StrokePtr s);
	void add(const Vec3& _p);
	StrokePtr newStroke();
	list<StrokePtr>* getStrokes();

	void drawGL();
	void HandleMouseEvent(int button, int state, const Vec3& p);
	void HandleMouseMotion(const Vec3& p);

	void insertStrokes(int i);
	void makeNURBS();
	void deleteLast();
	void reset();
	void processStrokes();
	void setCurves(list<CurvePtr>);

	static StrokeManager * getManager();
};

#endif