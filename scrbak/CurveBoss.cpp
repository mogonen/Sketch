#include "Curve2Mesh.h"

CurveBoss* CurveBoss:: _boss;
CurveBoss::CurveBoss(){
	_step = 0;
	_sm = StrokeManager::getManager();
	_c2s = new C2S();
	_c2m = new F2M();
}

CurveBoss* CurveBoss::getBoss(){
	if (!_boss)
		_boss =  new CurveBoss();
	return _boss;
}

void CurveBoss::drawGL(){
	if (_c2s)
		_c2s->drawGL();
	if (_c2m)
		_c2m->drawGL();
}

void CurveBoss::reset(){
	delete _c2s;
	delete _c2m;
	_c2s = new C2S();
	_c2m = new F2M();
	_step = 0;
}

void CurveBoss::step(){
	/*_c2s->resampleCurves();
	((K2M*)_c2m)->solveKnots(0.025);
	((K2M*)_c2m)->goForIt();
	_c2m->buildRings();
	_c2m->buildMesh(12);
	*/
	//cout<<"next step"<<endl;

	if (_step == 0 ){
		_c2s->resampleCurves();
		_c2s->buildCorrespondances();
	}

	if (_step == 1){
		_c2s->filter(3);
		_c2s->growCaps();
	}
	if (_step == 2 ){
		//_c2m->setSpine(_c2s->buildSpine());
		//_c2m->buildRings();
		((F2M*)_c2m)->setCP(_c2s->getCorrs().front());
		((F2M*)_c2m)->buildMesh();
		cout<<"2"<<endl;
	}

	if (_step == 3)
		((F2M*)_c2m)->blendMesh();

	_step++;
}

void CurveBoss::buildMesh(){
	_c2m->buildMesh(32);
}

void CurveBoss::buildSpine(){
	_c2m->setSpine(_c2s->buildSpine());
}

void CurveBoss::saveMesh(char * fname){
	_c2m->saveMesh(fname);
}

void CurveBoss::reproject(Camera * cam){
	_sm->setCurves(_c2m->getProjectionCurves(cam));
	//_c2m->buildRings();
}