#ifndef _H_C2M
#define _H_C2M

#include <map>
#include "Stroke.h"
#include "Mesh.h"
#include "Camera.h"
#include "Eye.h"
#include "Layer.h"

class CorrPoint;
class Ring;
#define CorrPtr CorrPoint*
#define RingPtr Ring*
#define STEP 0.005

static const double ZCOMPRESS = 0.25;
extern bool DRAW_AUX;
class CorrPoint:public CurvePoint{

	//corr
	CorrPtr _corr;
	bool _isshooter;
	unsigned int _flag;

	void updateNT();
	double _storeD;

	//spinal fields
	CurvePtr _spine;
	CorrPtr  _sptip;
	int		 _spdir;

public:

	CorrPoint(const Vec3& p, CurvePtr c);
	void drawAll();
	
	CorrPtr corr();
	CorrPtr tipCorr();
	CorrPtr nextCorr();
	CorrPtr prevCorr();
	CorrPtr corrNextCorr();
	CorrPtr corrPrevCorr();

	void setCorr(CPPtr cp);
	void discardCorr();

	bool isTip();
	bool isSplit();
	bool isIntersecting(CorrPtr cp);
	int  getFreeCount(int dir = 1);
	bool isShooter(){return _isshooter;};

	CorrPtr getShoot(Vec3 p, Vec3 n);
	void setSpine(CurvePtr sp, CorrPtr end);
	CurvePtr spine(){return _spine;};
	int spinalDir(){return _spdir;};
	CorrPtr tip(){return _sptip;};

	void storeD(double d){_storeD = d;};
	double retrieveD(){return _storeD;};

	void setFlag(unsigned int b){_flag |= (1 << b);};
	bool getFlag(unsigned int b){return _flag & (1 << b);};

	static CorrPtr sampleCurve(CurvePtr c, double step);
	static CorrPtr sampleCurve(Vec3 * p, int size, double step);
};

class C2S{

	list<CorrPtr> cps;
	StrokeManager* sman;

	CorrPtr spot2Shoot(CorrPtr cp);
	CorrPtr findCorr(CorrPtr);
	CorrPtr findCorr(CorrPtr, CorrPtr,  double&);
	ArrCurve* getSpine(CorrPtr cp, CorrPtr& endcp);

	void bumFilter(CorrPtr);
	void addCurve(CurvePtr c);
	void growCap(CorrPtr);

public:

	C2S();	
	void resampleCurves();
	void drawGL();

	void buildCorrespondances();
	void buildCorrespondance();
	void filter(int num);
	
	void bumFilter();
	void filterSpine();
	list<CurvePtr> buildSpine();
	static const Vec3& getOrgbyRad(CorrPtr, double);
	static double getZbyRadius(CorrPtr, double);
	
	list<CorrPtr> getCorrs(){return cps;};
	void growCaps();
};

class C2M{

protected:
	list<ArrCurve*> _curves;
	list<RingPtr>  _rings;
	list<CurvePtr> _spines;

	StrokeManager* _sman;

	MeshPtr _mesh;
	void buildMesh(RingPtr ring, double, int size=0);
	void buildMesh(RingPtr ring, int size);
	void buildRings(CurvePtr);

	Vec3 _eye;

public:

	ArrCurve* ac;

	C2M();
	void updateCurves();
	void setSpine(list<CurvePtr> &sp){_spines = sp;};
	void addSpine(CurvePtr sp){_spines.push_back(sp);};
	virtual void buildRings();
	void buildMesh(int);
	void buldMesh(double, int&);
	virtual void buildMesh(RingPtr ring, double len, int &size);

	list<CurvePtr> getProjectionCurves(Camera * cam);

	double getHit(Vec3 p, Vec3 n);
	void drawGL();
	void saveMesh(char * fname);
};

struct Cross;
//Knot to mesh
class K2M:public C2M{

	list<Cross> _crosses;
	Vec3 * solveKnots(int);

	int * updown;
	int getUpDown(int i);
	Vec3 * orgp;
	int size; 
	double rad; 

	int * _inds;
	int _curvenum;
	int skip;

public:
	K2M():C2M(){updown = 0;};
	void buildRings();
	void solveKnots(double r);
	void solveKnots2(double r);
	void goForIt();
};

struct Cross{
	ArrCurve* _c0;
	ArrCurve* _c1;
	int i0, i1;
	bool up0;
};

class F2M:public C2M{

	int _vnum;
	CorrPtr _cp;
	VertexPtr * createVerts(Vec3, Vec3);
	VertexPtr * createMidVerts(Vec3, Vec3);
	void buildMesh(CorrPtr cp, VertexPtr * v);
	void setVertexNormals(VertexPtr v0, VertexPtr v1);
public:
	void buildMesh();
	void setCP(CorrPtr cp){_cp = cp;};
	void blendMesh();
};

class CurveBoss:ShapeLayer{

	C2S * _c2s;
	C2M * _c2m;
	StrokeManager * _sm;

	static CurveBoss * _boss;
	int _step;
	Camera * _cam;

public:

	CurveBoss();
	void buildSpine();
	void buildMesh();

	void reset();
	void drawGL();

	void saveMesh(char *fname);
	void reproject(Camera * cam);

	static CurveBoss* getBoss();
	void step();
	void setCam(Camera * cam){ _cam = cam;};

};

class Ring{

	RingPtr _next;
	RingPtr _prev;

	CorrPtr _cp;
	Vec3 _p0, _p1, _p2, _n0;
	double _rad;

	void draw();
	static const int RSIZE = 24;

	Vec3 _org;
	double _orgrad;
	bool _isclosed;

public:

	Ring(CorrPtr cp);
	Ring(Vec3, Vec3, Vec3);
	Ring(Vec3, Vec3, double);
	RingPtr next(){return _next;};
	RingPtr prev(){return _prev;};

	Vec3 * getVertexPos(int size, int rot=0);
	Vec3 * getVertexPos(int size, Vec3 * prex, int rot=0);

	void setRadius(double rad);
	void reorient(Vec3 n);
	void setPrevNext(RingPtr _prev, RingPtr _next);

	double getNiceRadius();

	int getSegmentNum(double len);
	double getSegmentLength(int num);

	void drawGL();
	void drawAll();

	const Vec3& P(){return _p0;};
	const Vec3& N(){return _n0;};
	const Vec3& D(){return (_p2 - _p1);};
	const Vec3& R(){return (_p1 - _p0);};

	Vec3 * getSilhouettePoints(Vec3 p);

	void close(){_isclosed = true;};
	bool isClosed(){return _isclosed;};
};

class RingNode:public Ring{

	list<RingPtr> _rings;

public:
	void addRing(RingPtr r);
};

class CrossPoint:public ControlPoint{

	int * _state0;
	int * _state1;

public:
	void drawGL(){
		//glEnable(GL_POINT_SMOOTH);
		//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glPointSize(1.0);
		if (*_state0 > 0)
			glColor3f(0.0, 0.1, 0.5);
		else 
			glColor3f(1.0, 0.5, 0.1);

		glVertex3f(_p->x, _p->y, _p->z);
		//glVertex3f(pp.x, pp.y, pp.z);
		//glDisable(GL_POINT_SMOOTH);
	};

	void clicked(){*_state0 = -*_state0; *_state1 = -*_state1;cout<<"clicked!"<<(*_state0 )<<(*_state1)<<endl;};
	CrossPoint(int *state0, int *state1, Vec3 p ):ControlPoint(){_rad=0.02; _draggable=false; _state0 = state0; _state1 = state1; _p = new Vec3(p); ControlPoint::add(this);};
};

#endif