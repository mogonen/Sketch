#ifndef _MESH_H__
#define __MESH_H__

#include <list>
#include <vector>
#include <set>
//#include <GL/glut.h>
#include "Vector.h"
#include <QtOpenGL/QGLShaderProgram>

using namespace std;

#define EdgePtr Edge*
#define VertexPtr Vertex*
#define MeshPtr Mesh*
#define FacePtr Face*

class Mesh;
class Edge;
class Vertex;
class Face;

class Mesh{

	set<VertexPtr> _verts;
	list<FacePtr> _faces;
	//list<VertexPtr> edges;
public:

	//EdgePtr addEdge(VertexPtr v0, VertexPtr v1);

	VertexPtr addVertex(Vec3 p);
	void addFace(FacePtr);
	void addQuad(VertexPtr, VertexPtr, VertexPtr, VertexPtr);
	void addTriangle(VertexPtr, VertexPtr, VertexPtr);

	void drawGL();
	void drawGLTri();

	int exportOBJ(char * fname);
	Vec3 eye;

	void blendVertexNormals();
};

class Edge{
	VertexPtr v0; 
	VertexPtr v1; 
};

class Vertex{

	Vec3 _p;
	Vec3 _n0;

public:
	void setP(Vec3 p){_p = p;};
	void setN(Vec3 n){_n0 = n;};

	Vec3 getP(){return _p;};
	Vec3 getN(){return _n0;};

	float * getCoords(){return (float*)&_p[0];};
	vector<FacePtr> faces;
	set<VertexPtr> verts;
	Vec3 normal();

	void drawGL();
};



class Face{
	Vec3 _n0;
public:
	vector<VertexPtr> verts;
	MeshPtr mesh;
	Face(){};
	void addNextVertex(VertexPtr v){verts.push_back(v); v->faces.push_back(this);};
	//void addQuad(VertexPtr v0, VertexPtr v1, VertexPtr v2, VertexPtr v3);
	Vec3 normal();
};


#endif