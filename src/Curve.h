#ifndef __CURVE_H__
#define __CURVE_H__

#include <vector>
#include <stdio.h>
#include "Vector.h"
#include "Matrix.h"

using namespace std;

class Curve;

#define CurvePtr Curve*
#define CPPtr CurvePoint* 

class Curve{

protected:
	bool _isclosed;

public:

	Curve(){_isclosed = true;};

	virtual const Vec3& getP(double t)=0;
	virtual const Vec3& getT(double t);
	virtual double length();
	virtual int size()=0;
	virtual Vec3 * toArr(int len=0);
	virtual void drawGL();

	void close(){_isclosed = true;cout<<"closed"<<endl;};
	bool isClosed(){return _isclosed;};
};

class CurvePoint{

protected:
	CurvePtr _c;
	CPPtr _next;
	CPPtr _prev;
	Vec3 _p0, _t0, _n0;

	void updateNT();

public:

	CPPtr next(){return _next;};
	CPPtr prev(){return _prev;};
	CPPtr last();
	CPPtr go(int);
	Vec3 go(double);
	int find(CPPtr);

	void setPrevNext(CPPtr p, CPPtr n);
	Vec3 P(){return _p0;};
	Vec3 T(){return _t0;};
	Vec3 N(){return _n0;};

	void drawAll();
};

class ControlPoint{

	ControlPoint(CurvePtr c, Vec3 * p, int type = 0);

protected:
	CurvePtr _curve;
	bool _selected;
	bool _draggable;
	double _rad;
	Vec3* _p;
	static vector<ControlPoint*> _controls;
	int _type;


public:
	
	ControlPoint(){};

	int handleMouseEvent(int state, Vec3 p);
	void handleMouseMotion(Vec3 p);
	virtual void drawGL();
	virtual void clicked(){};

	static int handleMouseAll(int, Vec3);
	static void create(CurvePtr c, Vec3 * p, int typ=0);
	static void add(ControlPoint * cp);
	static void create(CurvePtr c, Vec3 p, int typ=0);
	static Vec3 (*toView)(Vec3 p);

	static void drawAll();
	static void clear();
	static void remove(CurvePtr c);
};

//a curve of equadistance array points
class ArrCurve:public Curve{

protected:
	Vec3 * _pts;
	int _size;
	double _len;
	double _step;
	bool _draw_aux;

public:
	ArrCurve(Vec3*, int size, bool resampleit = false);
	ArrCurve(ArrCurve&);
	ArrCurve(double *, int, Vec3 p, Vec3 nx, Vec3 ny, double w, double h);

	const Vec3& getP(double t);
	const Vec3& getP(int);
	const Vec3& getT(int);

	const Vec3& last(){return _pts[_size-1];};

	int size(){return _size;};

	double length();
	Vec3* toArr(int l=0){return _pts;};
	void drawGL();
	void drawGLAux();
	void drawAux(bool da = true){_draw_aux = da;};
	void setP(int i, Vec3 p){_pts[i].set(p);};

	void set(Vec3 * p, int s){delete [] _pts; _pts = p; _size = s;};
	void append(Vec3 * p, int s);

	//make an equadistance resampling
	void movingAverage(int num, bool norm=false);
	void resample(int size=0);
	void resample(double step);
	void reverse();

	static Vec3* resample(Vec3 * p, int size, double step, int& newsize);
	static Vec3* resample(Vec3 * p, int size, int outsize = 0);
	static double length(Vec3 * p, int size);
};

class ArrCurveUtils{

	ArrCurve* _ac;

public:
	double * maxhs;
	ArrCurveUtils(ArrCurve* ac){_ac = ac; };
	int * findMinMax(int step, int& mmi);
	double findMaxH(int i, int step, int& maxi);

	ArrCurve* getHCurve(int step, Vec3 p, int&);
};

class NURBS:public Curve{

	int _degree;
	int _num_knots;
	int _num_controls;
	int _size;

	double * _knots;   //the knots
	Vec3 * _controls;  //the controls points
	Vec3 * _curve;	   //the generated B-spline points
	Vec3 getIntermed(Vec3, Vec3, double, double, double);
	int LookUpKnots(double);
	void updateKnots(int deg);

public:

	NURBS(Vec3 *, int, int, double k);
	static NURBS* create(CurvePtr, double);
	int size(){return _size;};
	const Vec3& getP(double t);
};

class Bezier:public Curve{

	Vec3 * _pts;
	int _pnum;
	int _size;
	int _deg;

	Matrix _controls;	//control points
	Vec3* _cps;		//degree reduced control points
	int _num_controls;
	void Bezier::genControl(int deg, double * knots);
	
public:

	Bezier(Vec3 * cps, int num, int deg);
	const Vec3& getP(double t);
	int size(){return _size;};
	static Bezier* create(CurvePtr, int);
};

double berstein(int n, int i, double t);
int factorial(int i);


class Knot:public ArrCurve{ 

	double _rad;

public:
	Knot(ArrCurve* ac):ArrCurve(ac->toArr(), ac->size()){};
	void solve(double rad);
};

#endif