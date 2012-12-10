#include <iostream> 
#include <fstream>
#include <map>
#include "Mesh.h"
#include "Eye.h"

//Eye* Eye::get();

void Mesh::addFace(FacePtr f){
	f->mesh = this;
	_faces.push_back(f);
}

void Mesh::addQuad(VertexPtr v0, VertexPtr v1, VertexPtr v2, VertexPtr v3){
	FacePtr f = new Face();
	f->mesh = this;
	_faces.push_back(f);
	f->addNextVertex(v0);
	f->addNextVertex(v1);
	f->addNextVertex(v2);
	f->addNextVertex(v3);

	_verts.insert(v0);
	_verts.insert(v1);
	_verts.insert(v2);
	_verts.insert(v3);

	v0->verts.insert(v1);
	v0->verts.insert(v3);

	v1->verts.insert(v0);
	v1->verts.insert(v2);

	v2->verts.insert(v1);
	v2->verts.insert(v3);

	v3->verts.insert(v2);
	v3->verts.insert(v0);
}

void Mesh::addTriangle(VertexPtr v0, VertexPtr v1, VertexPtr v2){
	FacePtr f = new Face();
	f->mesh = this;
	_faces.push_back(f);
	f->addNextVertex(v0);
	f->addNextVertex(v1);
	f->addNextVertex(v2);

	_verts.insert(v0);
	_verts.insert(v1);
	_verts.insert(v2);

	v0->verts.insert(v1);
	v0->verts.insert(v2);

	v1->verts.insert(v0);
	v1->verts.insert(v2);

	v2->verts.insert(v0);
	v2->verts.insert(v1);
}

VertexPtr Mesh::addVertex(Vec3 p){
	VertexPtr v  = new Vertex();
	v->setP(p);
	_verts.insert(v);
	return v;
}

void Mesh::drawGL(){
	//drawGLTri();
	/*for(list<VertexPtr>::iterator itv = _verts.begin(); itv!=_verts.end(); itv++)
		(*itv)->drawGL();
		*/
	for(list<FacePtr>::iterator itf = _faces.begin(); itf!=_faces.end(); itf++){
		FacePtr f = (*itf);
		glColor3f(0, 0, 0);
		
		glBegin(GL_LINE_LOOP);
		for(vector<VertexPtr>::iterator itv = f->verts.begin(); itv!=f->verts.end(); itv++){
			Vec3 p = (*itv)->getP();
			glVertex3f(p.x, p.y, p.z);
		}
		glEnd();
		//*/
		//glEnable(GL_LIGHTING);
			glColor3f(0.5, 0.2, 0.2);
			glBegin(GL_POLYGON);
			for(vector<VertexPtr>::iterator itv = f->verts.begin(); itv!=f->verts.end(); itv++){
				Vec3 p = (*itv)->getP();
				Vec3 n = (*itv)->getN();
				//glNormal3f(n.x, n.y, n.z);
				glColor3f((n.x+1)/2.0 , (n.y+1)/2.0  , (n.z+1)/2.0 );
				//glColor3f(1, 0, 0  );
				glVertex3f(p.x, p.y, p.z);
			}
		glEnd();
		//glDisable(GL_LIGHTING);
	}
}

void Mesh::drawGLTri(){
	
	for(list<FacePtr>::iterator itf = _faces.begin(); itf!=_faces.end(); itf++){

		FacePtr f = (*itf);
				//glEnable(GL_LIGHTING);
		
		/*for(vector<VertexPtr>::iterator itv = f->verts.begin(); itv!=f->verts.end(); itv++){
			Vec3 p = (*itv)->getP();
			glVertex3f(p.x, p.y, p.z);
		}*/

		int j=0;
		int size = f->verts.size();
		//glColor3f(1.0f,0.0f,0.0f);  
		while(j < size - 1){

			glColor3f(0, 0, 0);
			glBegin(GL_LINE_LOOP);
			for(int i=j; i < j+3; i++){
				Vec3 p = f->verts[i%size]->getP();
				glVertex3f(p.x, p.y, p.z);
			}
			glEnd();
			

			Vec3 p0;
			for(int i=j; i < j+3; i++){
				p0 = p0+ f->verts[i%size]->getP();
			}

			Vec3 n0 = ((f->verts[(j+1)%size]->getP() -f->verts[(j)%size]->getP() ) % (f->verts[(j+2)%size]->getP() -f->verts[(j+1)%size]->getP() )).normalize();


			p0 = p0/3.0;
			//Vec3 ray = (p0 - Eye::get()->P).normalize();
			//double d = (ray*n0);
			//double d = Vec3(0,0,-1)*f->normal();
			/*
			glColor3f(0, 0 ,0);
			glBegin(GL_LINES);
				Vec3 np1 = p0 + f->normal()*0.1;
				glVertex3f(p0.x, p0.y, p0.z);
				//glColor3f(0, 0 ,1)
				glVertex3f(np1.x, np1.y, np1.z);
			glEnd();
			*/
			glBegin(GL_TRIANGLES);
			for(int i=j; i < j+3; i++){

				Vec3 p = f->verts[i%size]->getP();
				Vec3 ray = (p - Eye::get()->P).normalize();
				double d = (ray*f->verts[i%size]->normal());
				glColor3f( (d)>0?(d):0 , d<0?abs(d):0 ,0);
				
				glVertex3f(p.x, p.y, p.z);
			}
			glEnd();
			j+=2;
		}

		//glDisable(GL_LIGHTING);
	}
}

Vec3 Face::normal(){
	if (!_n0.norm()){
		_n0 = ( verts[1]->getP() - verts[0]->getP() ).normalize() % ( verts[2]->getP() - verts[1]->getP() ).normalize();
		_n0 = _n0.normalize();
	}
	return _n0;
}

Vec3 Vertex::normal(){
	if (!_n0.norm()){
		for(vector<FacePtr>::iterator itf = faces.begin(); itf!=faces.end(); itf++)
			_n0 = _n0 + (*itf)->normal();
		_n0 = _n0.normalize();
	}
	return _n0;
}

void Vertex::drawGL(){
	glColor3f(1.0, 0, 0);
	glBegin(GL_LINES);
		glVertex3f(_p.x, _p.y, _p.z);
		Vec3 p1 =_p + _n0*0.05;
		glVertex3f(p1.x, p1.y, p1.z);
	glEnd();
	/*_p.print();
	p1.print();
	cout<<"..."<<endl;*/
}

int Mesh::exportOBJ(char * fname){

  ofstream out(fname); 
  if(!out) { 
    cout << "Cannot open file.\n"; 
    return -1; 
  } 

  map<VertexPtr, int> vtable;

  out<<"#"<<fname<<endl;
  out<<"#vertices"<<endl;
  int ind = 1;
  for(set<VertexPtr>::iterator itv = _verts.begin(); itv != _verts.end(); itv++){
	  VertexPtr v = (*itv);
	  out<<"v "<< v->getP().x <<" "<< v->getP().y <<" "<< v->getP().z << endl;
	  vtable[v]=ind++;
  }

  out<<"#faces"<<endl;
  for(list<FacePtr>::iterator itf = _faces.begin(); itf!=_faces.end(); itf++){
	  out<<"f ";
	  for(vector<VertexPtr>::iterator itv = (*itf)->verts.begin(); itv != (*itf)->verts.end(); itv++)
		  out<<vtable[(*itv)]<<" ";
	  out<<endl;
  }

  out.close();
  cout <<" meash saved to "<<fname; 
  return 0;
}


void Mesh::blendVertexNormals(){	
	Vec3 * blended = new Vec3[_verts.size()];
	int i = 0;
	for(set<VertexPtr>::iterator itv = _verts.begin(); itv != _verts.end(); itv++){
		VertexPtr v = *itv;
		Vec3 nv(0,0,0);
		for(set<VertexPtr>::iterator itvv = v->verts.begin(); itvv != v->verts.end(); itvv++){
			nv = nv + (*itvv)->getN();
		}
		nv = nv / v->verts.size();
		blended[i++] = v->getN()*0.33 +  nv*0.67;
	}

	i = 0;
	for(set<VertexPtr>::iterator itv = _verts.begin(); itv != _verts.end(); itv++)
		(*itv)->setN(blended[i++]);
}