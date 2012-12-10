#include <stdlib.h>
//#include <GL/glut.h>
#include "Vector.h"
#include "Curve.h"
#include <QtOpenGL/QGLShaderProgram>

vector<ControlPoint*> ControlPoint::_controls;

Vec3 * Curve::toArr(int len){
	if(!len)
		len = size();
	Vec3 * arr = new Vec3[len];
	for(int i=0; i<len; i++)
		arr[i] = getP(i*1.0/(len-1.0));
	return arr;
}

const Vec3& Curve::getT(double t){
	double e = 0.025;
	double t1 = (t+e)>1.0?1.0:(t+e);
	double t0 = (t-e)<0?0:(t-e);
	return (getP(t1) - getP(t0));
}

double Curve::length(){
	double len=0;
	int SAMPLES = size();
	Vec3 p_pre;
	for(int i=0; i<SAMPLES; i++){
		Vec3 p = getP(i/(SAMPLES-1.0));
		if(i>0)
			len+=(p-p_pre).norm();
		p_pre.set(p);
	}
	return len;
}

void Curve::drawGL(){
	glBegin(GL_LINE_STRIP);
	for(int i=0; i<100; i++){
		Vec3 p = getP(i/99.0);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}

void CurvePoint::setPrevNext(CPPtr p, CPPtr n){
	_prev  = p;
	_next = n;
	if (_prev)
		_prev->_next = this;
	if (_next)
		_next->_prev = this;
}

CPPtr CurvePoint::go(int d){
	int inc = d>0?1:(d<0?-1:0);
	CPPtr cp = this;
	for(int i=0; i!=d; i+=inc){
		if (!cp)
			return cp;
		if (inc>0)
			cp = cp->next();
		else
			cp = cp->prev();
	}
	return cp;
}

CPPtr CurvePoint::last(){
	CPPtr last = this->next();
	for(;last->next();last = last->next())
		if (last->next() == this)
			break;
	return last;
}

int CurvePoint::find(CPPtr cp){

	CPPtr it_n = this;
	CPPtr it_p = this;
	int i_n = 0;
	int i_p = 0;

	while(it_n || it_p){

		if (it_n == cp)
			return i_n;

		if (it_p == cp)
			return i_p;

		if (it_n){
			it_n = it_n->_next;
			i_n++;
		}
		if (it_p){
			it_p = it_p->prev();
			i_p--;
		}
	}
	return 0;
}

Vec3 CurvePoint::go(double dist){
	if (dist*dist < 0.00001)
		return _p0;
	int dir =  (dist > 0 )?1:-1;
	CPPtr it0 = this;
	CPPtr it1 = this->go(dir);
	while(it0 && it1 && dist){
		Vec3 n = (it1->_p0 - it0->_p0);
		double d = n.norm();
		if ( d < dist)
			dist-=d;
		else{
			return it0->_p0 + n * (dist/d);
		}
		it0  = it1;
		it1 = it1->go(dir);
	}
	return it0->_p0;
}

void CurvePoint::updateNT(){
	_t0 = ((_next?_next->P():_p0) - (_prev?_prev->P():_p0)).normalize();
	//_n0 = (Eye::get()->N%_t0).normalize();
}

void CurvePoint::drawAll(){
glColor3f(0.0f,1.0f,0.0f);
	glPointSize(4.0f);
	glBegin(GL_POINTS);
	for(CPPtr cp = this; cp; cp = cp->next())	
		glVertex3f(cp->_p0.x, cp->_p0.y, cp->_p0.z);
	glEnd();
	glBegin(GL_LINE_STRIP);
	for(CPPtr cp = this; cp; cp = cp->next())	
		glVertex3f(cp->_p0.x, cp->_p0.y, cp->_p0.z);
	glEnd();

/*	glBegin(GL_LINES);
	for(CPPtr cp = this; cp; cp = cp->next()){			
		Vec3 p1 = cp->_p0 + cp->_n0*0.02;
		glVertex3f(cp->_p0.x, cp->_p0.y, cp->_p0.z);
		glVertex3f(p1.x, p1.y, p1.z);
	}
	glEnd();*/
	
}

NURBS::NURBS(Vec3 * cps, int size, int deg, double k){
	_size = size; 
	//_num_controls = 8;
	double len = ArrCurve::length(cps, size);
	double step = len*k + 0.005;
	step = 0.05;//(step>len)?len-0.0000001:step; 
	_controls = ArrCurve::resample(cps, size, step, _num_controls);

	for(int i=0; i<_num_controls; i++)
		ControlPoint::create(this, &_controls[i]);

	_degree    = deg;	
	_num_knots = _num_controls + _degree - 1;
	_knots = new double[_num_knots];
		
	for(int i=0; i<_degree; i++){
		_knots[i] = 0;
		_knots[_num_knots-i-1] = 1;
	}

	for(int i = _degree - 1; i<_num_controls-1; i++)
		_knots[i] = (i - _degree + 1.0)/(_num_controls - _degree);
	cout<<"new nurbs!"<<endl;
}

NURBS* NURBS::create(CurvePtr c, double k){
	return new NURBS(c->toArr(), c->size(), 2, k);
}

const Vec3& NURBS::getP(double t){
	Vec3 ** pyramid;
	pyramid = new Vec3*[_degree+1];
	int i,j=0, pos;
	pos = LookUpKnots(t);
	pyramid[0] = new Vec3[_degree+1];
	for(i=0;i<_degree+1;i++)
		pyramid[0][i] = _controls[i + pos - _degree + 1];
	for(i=1; i<_degree+1; i++){
		pyramid[i] = new Vec3[_degree+1-i];
		for(j=0;j<_degree+1-i;j++)
			pyramid[i][j] = getIntermed(pyramid[i-1][j], pyramid[i-1][j+1], t, _knots[pos - _degree+i+j], _knots[pos+1+j]);
	}
	return pyramid[i-1][j-1];
}

int NURBS::LookUpKnots(double t){
	int i = _degree-1;
	while(_knots[i]<t)			
		i++;
	if(t<_knots[0]+0.00001&&t>_knots[0]-0.00001)
		return i;
	return i-1;
}

void NURBS::updateKnots(int deg){
	if(_degree == deg)
		return;
		//generate new pesudo knots
		int i;
		_num_knots = _num_controls+_degree-1;
		double * new_knots = new double[_num_knots];			
		for(i=0;i<deg;i++){
			new_knots[i] = 0;
			new_knots[_num_knots-i-1] = 1;
		}
			
		for(i=deg-1;i<_num_controls+deg-1;i++)
			new_knots[i] = _knots[i-deg+_degree];
			
		_knots = new_knots;
		_degree = deg;	
}

Vec3 NURBS::getIntermed(Vec3 p1, Vec3 p2, double t, double t1, double t2) {
	if(t2==t1)
		return p1;
	return (p1*(t2-t) + p2*(t-t1)) / (t2-t1);		
}

ArrCurve::ArrCurve(Vec3 * p, int s, bool resampleit){

	if (resampleit)
		_pts = ArrCurve::resample(p, s,0);
	else 
		_pts = p;
	cout<<"size:"<<s<<endl;
	_size = s;
	_len = 0;
	length();
	_draw_aux = false;
	_isclosed = false;
}

ArrCurve::ArrCurve(ArrCurve& ac){
	_size = ac.size();
	_pts = new Vec3[_size];
	for(int i = 0; i<_size; i++ )
		_pts[i] = ac._pts[i];
	_len = 0;
	length();
	_draw_aux = false;
	_isclosed = ac.isClosed();
}

ArrCurve::ArrCurve(double* hs, int size, Vec3 pos, Vec3 nx, Vec3 ny, double width, double height){

	_draw_aux = false;
	double min =  9999999; 
	double max = -9999999;
	for(int i=0; i<size; i++){
		if (hs[i] < min)
			min = hs[i];
		if (hs[i] > max)
			max = hs[i];
	}

	_pts = new Vec3[size];
	for(int i=0; i<size; i++)
		_pts[i] = pos + ny*( (hs[i]/(max-min))*height) + nx*(((i*1.0)/size)*width);

	_size  = size;
	cout<<"buuuum!"<<size<<endl;

}

const Vec3& ArrCurve::getP(double t){
	t = (t>1)?1.0:((t<0)?0:t);
	double tt =  (t*(_size-1));
	int i = (int)tt;
	double ttt = tt - i;
	return _pts[i]*(1-ttt) + _pts[i+1]*ttt;
}

const Vec3& ArrCurve::getP(int i){

	if (_isclosed)
		i = (i+_size)%_size;
	else
		i = (i>=_size)?(_size-1):((i<0)?0:i);

	return _pts[i];
}

const Vec3& ArrCurve::getT(int t){


	int t1 = (t+1>=_size)?(_size-1):t+1;
	int t0 = (t-1<0)?0:(t-1);

	if (_isclosed){
		t1 = (t+1+_size)%_size;
		t0 = (t-1+_size)%_size;
	}

	if (t ==0 || t == _size)
		cout<<t0<<"  "<<t<<" "<<t1<<" "<<_isclosed<<endl;

	return (_pts[t1] - _pts[t0]).normalize();
}

double ArrCurve::length(){
	double len = 0;
	for(int i=1; i<_size; i++)
		len+=(_pts[i]-_pts[i-1]).norm();
	_len = len;
	return len;
}

void ArrCurve::drawGL(){
	if (_draw_aux)
		drawGLAux();
	else{
		glBegin(GL_LINE_STRIP);
		for(int i=0; i<_size; i++)
			glVertex3f(_pts[i].x, _pts[i].y, _pts[i].z);
		if (_isclosed)
			glVertex3f(_pts[0].x, _pts[0].y, _pts[0].z);
		glEnd();
	}
}

void ArrCurve::drawGLAux(){

	glColor3f(1.0, 0, 0);
	glPointSize(3.0);
	glBegin(GL_POINTS);
	for(int i=0; i<_size; i++)
		glVertex3f(_pts[i].x, _pts[i].y, _pts[i].z);
	glEnd();

	glColor3f(0, 0, 0);
	glBegin(GL_LINE_STRIP);
	for(int i=0; i<_size; i++)
		glVertex3f(_pts[i].x, _pts[i].y, _pts[i].z);
	glEnd();

	glColor3f(0, 0, 1.0);
	glBegin(GL_LINES);
	for(int i=0; i<_size; i++){
		Vec3 t = getT(i);
		Vec3 p1 = _pts[i] + t*0.02;
		//_pts[i].print();p1.print();cout<<endl;
		glVertex3f(_pts[i].x, _pts[i].y, _pts[i].z);
		glVertex3f(p1.x, p1.y, p1.z);
	}
	glEnd();

}

void ArrCurve::resample(int size){
	if (size<=0)
		size = _size;
	_len = length();
	_step = _len/(size-1.0);

	Vec3 * newp = ArrCurve::resample(_pts, _size);
	delete [] _pts;
	_pts = newp;
	_size = size;
	length();
}

void ArrCurve::resample(double step){
	int newsize = 0;
	Vec3 * newp = ArrCurve::resample(_pts, _size, step, newsize);
	delete [] _pts;
	_pts = newp;
	_size = newsize;
	length();
}

Vec3* ArrCurve::resample(Vec3 * p, int size, int newsize){
	return ArrCurve::resample(p, size, 0, (!newsize)?size:newsize );
}

Vec3* ArrCurve::resample(Vec3 * p, int size, double step, int& newsize){

	double len = ArrCurve::length(p, size);
	if (!step)
		step = len/(newsize-1.0);
	else 
		newsize = (int)(len/step)+1;

	cout<<"sampling:"<<size<<" step:"<<step<<" ns:"<<newsize<<" l:"<<len<<endl;

	Vec3 * newp = new Vec3[newsize];
	int j = 1;
	double d_tot = 0;
	Vec3 pre(p[0]);
	newp[0].set(pre);
	for(int i = 1; i < newsize; i++){
		double d = 0;
		while(d_tot<step && j<size){
			d = (p[j]-pre).norm();
			d_tot+=d;
			if (d_tot<step){
				pre.set(p[j]);
				j++;
			}
		}
		double t = (d_tot - step) / d;
		newp[i]  = p[j]*(1-t) + pre*t;
		pre.set(newp[i]);
		d_tot = 0;
	}
	newp[newsize-1]  = p[size-1];
	return newp;
}

double ArrCurve::length(Vec3 *p, int size){
	double len = 0;
	for(int i=1; i<size; i++){
		len+=(p[i]-p[i-1]).norm();
		//cout<<i<<" -l:"<<len<<endl;
	}return len;
}

void ArrCurve::movingAverage(int num, bool norm){

	int end = (_isclosed)?_size+1:_size-1;
	for(int j=0; j<num; j++){
		Vec3 p0(_pts[0]); 
		Vec3 p1(_pts[1]);
		
		for(int i = 1; i<end; i++){
			int ind = i%_size;
			Vec3 p2(_pts[(i+1)%_size]);
			if (norm)
				_pts[ind] = p1.normalize()*((p0 + p1 + p2 ).norm()/3.0);
			else
				_pts[ind] = (p0 + p1 + p2 ) / 3.0; 
			p0 = p1;
			p1 = p2;
		}
	}
}

void ArrCurve::reverse(){
	for(int i =0; i<_size/2; i++){
		Vec3 tmp = _pts[i];
		_pts[i] = _pts[_size-i-1];
		_pts[_size-i-1] = tmp; 
	}
}

void ArrCurve::append(Vec3 * p, int s){

	Vec3 * newp = new Vec3[_size + s];
	for(int i=0; i<_size; i++)
		newp[i] = _pts[i];

	for(int i=0; i<s; i++)
		newp[i+_size] = p[i];
	delete [] _pts;
	_pts = newp;
	_size+= s;
}

Bezier* Bezier::create(CurvePtr c, int deg){
	return new Bezier(c->toArr(), c->size(), deg);
}


Bezier::Bezier(Vec3 * cps, int num, int deg){
	_deg = deg;//3;
	_size = 100;
	_num_controls = num;
	_controls = Matrix(_num_controls,3);
	//copy control points
	for(int i=0;i<_num_controls; i++)
		_controls[i] = cps[i];

	//generate pesudo knots, uniformly-distributed
	double * knots = new double[_num_controls];
	for(int i=0; i<_num_controls; i++)
		knots[i] = i/(_num_controls-1.0);
	genControl(_deg, knots);
}

void Bezier::genControl(int deg, double * knots){

	Matrix C = Matrix(_num_controls, deg+1);
	for(int i=0; i<_num_controls; i++)
		for(int j=0; j<deg+1; j++)
				C[i][j] = berstein(deg, j, knots[i]);
			
	Matrix U, V;
	Vector D;
	C.svd(U, D, V);
	Matrix invC = ((V*diag(D).inv())*U.transpose());
	//C.print(); cout<<"---"<<endl; invC.print();
	
	Matrix _deg_controls = invC *_controls;
	//cout<<"="<<endl;
	
	_cps = new Vec3[_deg_controls.nrows()];
	for(int i=0; i<_deg_controls.nrows();i++){
		_cps[i] = Vec3(_deg_controls[i]);
		ControlPoint::create(this, &_cps[i]);
	}
}

//n-total degree, i-current points, t-time
double berstein(int n, int i, double t){
	return (factorial(n)/factorial(n-i)/factorial(i))*pow(1-t,n-i)*pow(t,i);
}
	
int factorial(int i){
	int x = 1;
	while(i>0)
		x*=i--;
	return x;
}	

const Vec3& Bezier::getP(double t){
	//cout<<"t:"<<t<<endl;
	int i;
	Vec3 p(0,0,0);
	for(i=0; i < _deg+1; i++){
		double b = berstein(_deg, i, t);
		p = p + b*_cps[i];
	}
	//p.print();cout<<", t:"<<t<<endl;
	return Vec3(p.x, p.y, p.z);
}

ControlPoint:: ControlPoint(CurvePtr c, Vec3 * p, int tip){ 
	_curve = c; _p = p; 
	_selected = false;
	_rad = 0.015;
	if (tip)
		_draggable = false;
	else
		_draggable = true;
}

int ControlPoint::handleMouseEvent(int state, Vec3 p){
	_selected = false;
	Vec3 pp = toView(*_p);
	double rad = (p - pp).norm();
	if (rad>_rad || !state)
		return 0;
	_selected = true;
	return 1;
}

void ControlPoint::handleMouseMotion(Vec3 p){
	if (!_selected || !_draggable)
		return;
	Vec3 pre = toView(*_p);
	_p->set( (*_p) + (p-pre) );
}

int ControlPoint::handleMouseAll(int state, Vec3 p){
	int i = 0;
	for (vector<ControlPoint*>::iterator it = ControlPoint::_controls.begin(); it != ControlPoint::_controls.end(); it++)
		if (state==0)
			(*it)->handleMouseMotion(p);
		else{
			int clicked = (*it)->handleMouseEvent(state, p);
			if (clicked && state == 2)
				(*it)->clicked();
			i = i || clicked;
		}
	return i;
}

void ControlPoint::create(CurvePtr c, Vec3 *p, int type){
	ControlPoint::_controls.push_back(new ControlPoint(c, p, type));
}

void ControlPoint::create(CurvePtr c, Vec3 p, int type){
	ControlPoint::_controls.push_back(new ControlPoint(c, new Vec3(p), type));
}

void ControlPoint::add(ControlPoint * cp){
	ControlPoint::_controls.push_back(cp);
}

void ControlPoint::drawGL(){
	if (_draggable){
			glColor3f(1.0, 0.7, 0);
			glPointSize(8.0);
		} else {
			glColor3f(1.0, 0.1, 0);
			//glPointSize(6.0);
		}
		glVertex3f(_p->x, _p->y, _p->z);
}

void ControlPoint::drawAll(){
	
	glColor3f(1.0, 0.7, 0);
	glPointSize(8.0);
	glBegin(GL_POINTS);
		for (vector<ControlPoint*>::iterator it = ControlPoint::_controls.begin(); it != ControlPoint::_controls.end(); it++)
			(*it)->drawGL();
	glEnd();
}

void ControlPoint::clear(){
	_controls.clear();
}

void ControlPoint::remove(CurvePtr c){
	vector<ControlPoint*>::iterator it = ControlPoint::_controls.begin(); 		
	while(it != ControlPoint::_controls.end()){
		if ((*it)->_curve == c){
			ControlPoint * cp = (*it);
			it = _controls.erase(it);
			delete cp;
		} else
			it++;
	}
}

double ArrCurveUtils::findMaxH(int i0, int i1, int& maxi){

	Vec3 p0 = _ac->getP(i0);
	Vec3 n = (_ac->getP(i1) - p0);
	double len = n.norm();
	n = n.normalize();

	double maxh = 0;
	for(int i = i0+1; i < i1; i++){
		Vec3 n1 = (_ac->getP(i) - p0);
		double hip = n1.norm();
		double a = n1*n;
		double h = sqrt(hip*hip - a*a);
		if (maxh < h){
			maxh = h;
			maxi = i;
		}
	}
	//cout<<endl<<"mnaxh:"<<maxh<<endl;
	return maxh;
}

int * ArrCurveUtils::findMinMax(int step, int& mmi){

	int i1 = _ac->size()-step;
	maxhs = new double[i1];
	int * maxis = new int[i1];

	for(int i = 0; i < i1; i++)
		maxhs[i] = findMaxH(i, i+step, maxis[i]);

	double  * maxhs_f = new double[i1];
	for(int k=0;k<3; k++){
		for(int i=1; i < i1-1; i++)
			maxhs_f[i] = (maxhs[i-1]+maxhs[i+1]+maxhs[i+1])/3.0;
		maxhs = maxhs_f;
	}

	/*int * maxmins = new int[20];
	mmi = 0;
	for(int i = 1; i < i1-1; i++)
		if ( maxhs[i] > maxhs[i-1] && maxhs[i] > maxhs[i-1]  )
			maxmins[mmi++] = maxis[i];
	
	int * newmaxmins = new int[mmi];
	for(int i = 0; i < mmi; i++)
		newmaxmins[i] = maxmins[i];
		*/

	mmi = 1;
	double hmax=0;
	int * newmaxmins = new int[1];
	for(int i=0; i<i1-1; i++)
		if (hmax < maxhs[i]){
			hmax = maxhs[i];
			newmaxmins[0] = maxis[i];
		}

	return newmaxmins;
}

ArrCurve* ArrCurveUtils::getHCurve(int step, Vec3 pos, int &locmax){
	int i1 = _ac->size()-step;
	maxhs = new double[i1];
	int * maxis = new int[i1];

	for(int i = 0; i < i1; i++)
		maxhs[i] = findMaxH(i, i+step, maxis[i]);

	double  * maxhs_f = new double[i1];
	for(int k=0; k<10; k++ ){
		for(int i=1; i < i1; i++)
			maxhs_f[i] = (maxhs[i-1] + maxhs[i] + maxhs[i+1])/3.0;
	
		maxhs_f[0] = maxhs[0];
		maxhs_f[i1-1] = maxhs[i1-1];
		maxhs = maxhs_f;
	}
	locmax = 0;
	for(int i=1; i<i1-1; i++)
		if ( (maxhs[i] > maxhs[i-1]) && (maxhs[i] > maxhs[i+1]) )
			locmax++;

	cout<<"local maxes:"<<locmax<<endl;

	return new ArrCurve(maxhs, i1, pos, Vec3(1,0,0), Vec3(0,1,0), 0.5, 0.2);
}

void Knot::solve(double r){

	Vec3 up(0.0, 0.0, 1.0);
	double eps = 0.0001;
	_rad = r;
	double diam = _rad*2;
	int istop = 1;
	int precoll = 0;
	int lowest = 5;
	resample(diam/2 + eps);

	int * colls = new int[_size];
	for(int i=0; i < _size; i++)
		colls[i] = 0;
	Vec3 ax(0,0,1);
	Vec3 * newp = new Vec3[_size];
	for(int i = 0; i < _size; i++)
		newp[i] = _pts[i];
	
	bool crossed = false;

	for(int i = lowest; i < _size; i++){

		Vec3 dir = (_pts[i]-_pts[i-2]);
		double len  = dir.norm();
		Vec3 ncross = dir%ax;
		dir = dir.normalize();

		for(int j=lowest; j < _size; j++){

			if (abs(i-j) < lowest)
				continue;
			if (!crossed){
				/*Vec3 dir2 = (_pts[j] - _pts[j-2]);
				double len2 = dir2.norm();
				dir2 = dir2.normalize();
				double d = getIntersectionDist(Vec2(dir),Vec2(_pts[i-2]), Vec2(dir2), Vec2(_pts[j-2]));
				double d2 = getIntersectionDist(Vec2(dir2),Vec2(_pts[j-2]), Vec2(dir), Vec2(_pts[i-2]));

				if ( (d>=0 && d<=len) && (d2>=0 && d2<=len2) )
					crossed = true;
					*/
				Vec3 v0 = _pts[j-2] - _pts[i-1];
				Vec3 v1 = _pts[j] - _pts[i-1];
				double d0 = v0*dir;  
				double d1 = v1 * dir;

				if ( ((d0 < len) && (d1 < len) && (d0 >0 ) && (d1 > 0)) && ((v0*ncross)*(v1*ncross)<0) )
					crossed = true;
			}

			Vec3 v =  (newp[i] - newp[j]);
			double d = v.norm();
			if (d + eps < diam ){
				newp[i]  = newp[i] + up*(sqrt(diam*diam - d*d)*istop);
				colls[i] = j;
				colls[j] = i;
			}

		}

		bool collision_end = precoll && (abs(colls[i] - precoll) > lowest);

		if (collision_end){
			cout<<"corss:"<<crossed<<endl;
			//if (crossed)
				istop = istop*-1;
			crossed = false;
		}

		precoll = colls[i];
	}

	for(int i=0; i<_size; i++)
		cout<<colls[i]<<" ";

	delete [] _pts;
	_pts = newp;
	resample(diam/2);

	cout<<endl;
	this->movingAverage(4);
}