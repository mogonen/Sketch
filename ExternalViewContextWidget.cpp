/***************************************************************************\

  NAME:         EXTERNALVIEWCONTEXT.CPP

  DESCRIPTION:  This sample shows how to instantiate and animate some simple
                mesh objects

  NOTE:         This file is part of the HueSpace visualization API.
                Copyright (C) 2001-2011 Hue AS. All rights reserved.

\***************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////

#include "ExternalViewContextWidget.h"
#include <math.h>
#include "RenderView.h"
#include "RootObject.h"
#include "Workspace.h"
#include "SceneManager.h"
#include "Scene.h"
#include "ProjectManager.h"
#include "Project.h"
#include "VMVolumeProgram.h"
#include "ProjectLibrary.h"
#include "StringObjectMask.h"
#include "ObjectMaskManager.h"
#include "Viewer.h"
#include "ViewerManager.h"
#include "Camera.h"
#include "CameraManagerSingle.h"
#include "ViewContextContainer.h"
#include "ViewContextContainerManagerSingle.h"
#include "WorldBackgroundViewContext.h"
#include "ViewContextManager.h"
#include "ShapeInstanceViewContext.h"
#include "ExternalViewContext.h"
#include "MarkerShape.h"
#include "MarkerShapeManagerSingle.h"
#include "OrbitCameraController.h"
#include "OrbitCameraControllerManagerSingle.h"
#include "FlyCameraController.h"
#include "FlyCameraControllerManagerSingle.h"
#include "Edit3DController.h"
#include "Edit3DControllerManagerSingle.h"
#include "PaintPolylineController.h"
#include "PaintPolylineControllerManagerSingle.h"
#include "EditTransferFunctionController.h"
#include "EditTransferFunctionControllerManagerSingle.h"
#include "ViewSplitController.h"
#include "ViewSplitControllerManagerSingle.h"
#include "RequestVCExternalEventContext.h"
#include "DeliverVCExternalEventContext.h"
#include "RenderOpaqueExternalEventContext.h"
#include "RenderTransparentExternalEventContext.h"
#include "RenderDepthPeelExternalEventContext.h"
#include "DeInitExternalEventContext.h"
#include "Trace3DExternalEventContext.h"
#include "IntersectionExternalManager.h"
#include "RenderLayer3D.h"
#include "HueQtDefaultViewerWidget.h"
#include "VDSGenNoise.h"
#include "VDSManager.h"
#include "MovableObjectManager.h"

#include "RootObject.h"
#include "Workspace.h"
#include "SceneManager.h"
#include "Scene.h"
#include "ProjectManager.h"
#include "Project.h"
#include "VMVolumeProgram.h"
#include "ProjectLibrary.h"
#include "StringObjectMask.h"
#include "ObjectMaskManager.h"
#include "Viewer.h"
#include "ViewerManager.h"
#include "Camera.h"
#include "MovableObjectManager.h"
#include "Horizon.h"
#include "VDSGenNoise.h"
#include "VDSManager.h"
#include "Edit3DController.h"
#include "ControllerManager.h"
#include "OrbitCameraController.h"
#include "TransferFunction.h"
#include "TransferFunctionManager.h"
#include "RenderLayerScreengrab.h"
#include "RenderLayerManager.h"
#include "RenderLayerSupersample.h"
#include "RenderLayerManagerSingle.h"
#include "RenderLayer3D.h"
#include "ViewContextContainer.h"
#include "ViewContextManager.h"
#include "ShapeInstanceViewContext.h"
#include "WorldBackgroundViewContext.h"
#include "HorizonViewContext.h"
#include "TorusShape.h"
#include "MovableObjectManager.h"
#include "ShapeInstance.h"
#include "ViewContextManager.h"
#include "ViewContextContainer.h"
#include "ShapeInstanceViewContext.h"
#include "HueQtDefaultViewerWidget.h"
#include "Workspace.h"
#include "SceneManager.h"
#include "Scene.h"
#include "ProjectManager.h"
#include "Project.h"
#include "ViewerManager.h"
#include "Viewer.h"
#include "ShapeManager.h"
#include "VolumeBox.h"
#include "VolumeSlice.h"
#include "RenderStateTree.h"
#include "VolumeRSManager.h"
#include "RenderStateTreeManager.h"
#include "RSTViewContext.h"
#include "SetFrontFaceStateNode.h"
#include "RSTNodeManagerSingle.h"
#include "DirectShape.h"
#include "ShrinkWrapShape.h"
#include "Trace2DEventContext.h"
#include "PaintPolylineController.h"

#include "Property.h"
#include "PropertyManager.h"
#include "PropertyData.h"




#include <QtOpenGL/QGLShaderProgram>
#include <gl/glu.h>
#include <QtGui/QMouseEvent>
// forward-declaration
void ExternalViewContext_RenderOpaqueExternalHandler(Hue::ProxyLib::RenderOpaqueExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID);
void ExternalViewContext_RenderTransparentExternalHandler(Hue::ProxyLib::RenderTransparentExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID);
void ExternalViewContext_RenderDepthPeelExternalHandler(Hue::ProxyLib::RenderDepthPeelExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID);
void ExternalViewContext_Trace3DExternalHandler(Hue::ProxyLib::Trace3DExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID);
void ExternalViewContext_DeInitExternalHandler(Hue::ProxyLib::DeInitExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID);
static void InitSphere(float radius);

using namespace Hue::ProxyLib;

Hue::ProxyLib::Viewer * ExternalViewContextWidget::pViewer;
StrokeManager * ExternalViewContextWidget::strokeman;

// METHODS //////////////////////////////////////////////////////////////////
Vec3 (*ControlPoint::toView)(Vec3 p);

/////////////////////////////////////////////////////////////////////////////
// ExternalViewContextWidget constructor
Vec3 toView(Vec3 p){
	Hue::Util::DoubleMatrix4x4
		cModelViewMatrix = ExternalViewContextWidget::pViewer->RenderLayer3DOutsideMagnifier()->LastModelViewMatrix(),
	  cProjMatrix =  ExternalViewContextWidget::pViewer->RenderLayer3DOutsideMagnifier()->LastProjectionMatrix();

  Hue::Util::IntVector2
	  cView = ExternalViewContextWidget::pViewer->RenderView()->Size();
  GLdouble vx, vy, vz;

  GLint viewport[4] =  {0, 0, cView.X, cView.Y };

  gluProject(p.x, p.y, p.z, (double *)&cModelViewMatrix, (double *)&cProjMatrix, viewport, &vx, &vy, &vz);

  return Vec3(vx, vy, vz);
}

Vec3 toWorld(int x, int y){

	Hue::Util::DoubleMatrix4x4
		cModelViewMatrix = ExternalViewContextWidget::pViewer->RenderLayer3DOutsideMagnifier()->LastModelViewMatrix(),
	  cProjMatrix =  ExternalViewContextWidget::pViewer->RenderLayer3DOutsideMagnifier()->LastProjectionMatrix();

 GLdouble ProjMatrix[16];
  glGetDoublev(GL_PROJECTION_MATRIX, ProjMatrix);
  GLdouble mvMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

  Hue::Util::IntVector2
 	  cView = ExternalViewContextWidget::pViewer->RenderView()->Size();

	double winX = (double)x;
	double winY = (double)cView.Y - (double)y;
	double winZ = 0;
	double wx, wy, wz;
	GLint viewport[4] = {0, 0, cView.X, cView.Y };

	gluUnProject( winX, winY, winZ, (double*)&cModelViewMatrix, (double*)&cProjMatrix , viewport, &wx, &wy, &wz);
 	return Vec3(wx, wy, wz);
}

ExternalViewContextWidget::ExternalViewContextWidget(QWidget* pParent, int argc, char ** argv) : HueQtDefaultViewerWidget(pParent)
{

  sketching_on = false;
  ControlPoint::toView = toView;
  ExternalViewContextWidget::strokeman = StrokeManager::getManager();
  // Find root scene object:
  Hue::ProxyLib::Scene *
    pScene = Hue::ProxyLib::Workspace::Instance()->Scenes()[0];

  Hue::ProxyLib::Project *
    pProject = pScene->Projects().Create("SampleProject");

  // Create viewer:
    pViewer = pProject->Viewers().Create(); // This is a viewer object which represents a viewport in a running Huespace-enabled application. It is useful for holding on to render layers, controllers and other objects which might conveniently be associated with a single viewport.The viewer is also responsible for receiving analog XY input and distributing it among the associated controllers. When SetPivotPointOnDown is true, the viewer directly handles AnalogXYInputAction.Down to set its OrbitCameraController's pivot point location to the analog XY input position.

  pViewer->SetHandleTestViewer(true); // At TRUE this object is prepared to handle events from Hue's internal test application(HueViewer). Keep at FALSE in any other framework.

 
//  gluProject

  // Create camera:
  Hue::ProxyLib::Camera *
    pCamera = pViewer->OwnedCamera(); // A camera is HueViewer's way of defining a vantage point. It can be moved and rotated freely in 3D, and also defines the field of view (zoom). In order for a camera to used, it must be assigned to a Render layer 3D

  pCamera->SetName("Camera"); // This property contains the name of the object. This is provided only for the benefit of the user, as all object relationships are internally resolved using the unique IDs of each object. (Default: "Camera (2)")
  pCamera->SetPosition(DoubleVector3(-0.843, -4.856, 5.136)); // This is the 3D position of the object, relative to the parent object in the object tree. This means that if you move or rotate the parent object (or another object higher up in the hierarchy) this object will move/rotate accordingly. If you want the object to be fixed in place independently of all other objects, place it as a direct child of the Objects container.
  pCamera->SetOrientation(DoubleMatrix3x3( 0.9701, -0.2426,  0.0,
                                           0.1562,  0.6247,  0.7650,
                                          -0.1856, -0.7421,  0.6439)); // This is the 3D orientation of the object, relative to the parent object in the object tree. This means that if you move or rotate the parent object (or another object higher up in the hierarchy) this object will move/rotate accordingly. If you want the object to be fixed in place independently of all other objects, place it as a direct child of the Objects container.
  pCamera->SetLocalMatrix(DoubleMatrix4x4( 0.9701, -0.2426,  0.0,  0.0,
                                           0.1562,  0.6247,  0.7650,  0.0,
                                          -0.1856, -0.7421,  0.6439,  0.0,
                                          -0.8434, -4.8561,  5.1364,  1.0)); // This property defines the local matrix of this object, i.e. its location relative to the parent object
  pCamera->SetFOV(90.0); // The vertical field of view of the camera (in degrees) when using 3D projection (Default: 60.0)

  
  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer = pViewer->ViewContextContainerInsideMagnifier(); // This is a view context container object which can hold on to one or more view contexts. By e.g. referring to this from a ViewContextInclude object, you can easily reuse view contexts across multiple Render layers.

  pViewContextContainer->SetName("Inside magnifier (left)"); // This parameter contains the name of this View context container object
  SetActiveViewer(pViewer);


pViewer->Edit3DController()->SetActive(true);
  pViewer->Edit3DController()->SetEdit3DTargetType(Edit3DTargetType::ThreeDObject);

  //Create noise (small procedural volume)
  Hue::ProxyLib::VDSGenNoise *
    vdsNoise = pProject->VDSs().CreateVDSGenNoise();

  Hue::ProxyLib::VDS *
    //restoreVDS = pProject->RestoreVDSFromFileName("c://project//PQT//data//Stack_Training");
    restoreVDS = pProject->RestoreVDSFromFileName("c://converted//entenschnabel_poststack");
//*/
  Hue::ProxyLib::VolumeBox *
	  volumeBox = pProject->MovableObjects().CreateVolumeBox();

  volumeBox->SetValueVDS(restoreVDS);
  volumeBox->SetScaleToCubicVoxelsMode(Hue::ProxyLib::ScaleToCubicVoxelsMode::Once);
  //volumeBox->SetValueReadoutMethod(Hue::ProxyLib::ReadoutMethod::Discrete);
  volumeBox->SetVBCacheBias(0);
  
  Hue::ProxyLib::VolumeSlice *
    volumeSlice = volumeBox->Attached3DObj().CreateVolumeSlice();
  Hue::ProxyLib::VolumeSlice *
    volumeSlice2 = volumeBox->Attached3DObj().CreateVolumeSlice();

 volumeSlice2->SetPlane(AxisAlignedPlane::YZ);

 volumeSlice->SetCacheBias(1);
 volumeSlice2->SetCacheBias(1);

 volumeBox->SetOrientation(DoubleMatrix3x3(DoubleVector3(0.0,0.0,-1.0),
                                           DoubleVector3(0.0,1.0,0.0),
			                   DoubleVector3(1.0,0.0,0.0)));
 
 
 Hue::ProxyLib::TransferFunction *
    pTransferFunction = pProject->TransferFunctions().Create(); // A transfer function is an object which contains information about how to map from raw, original data values to final colors and transparency. They can be applied to volume boxes and mesh surfaces using VM programs.

 Hue::ProxyLib::VolumeRS *
    pVolumeRS = pProject->VolumeRenderStates().Create(); // This is a volume render state which defines how volume rendering should be performed.

  StringObjectMask
    *pStringObjectMask1D1SReadout = RootObject::Instance()->ObjectMasks().CreateStringObjectMask(SubObjectEnumString::NamedObject_Name, "1D1S_Readout");

  Hue::ProxyLib::ProxyObjectList<VMVolumeProgram*>
    vmVolumePrograms1D1SReadout = VMVolumeProgram::StaticFindAll(ProjectLibrary::Instance(), true, pStringObjectMask1D1SReadout);

  pVolumeRS->SetVMVolumeProgram(vmVolumePrograms1D1SReadout[0]); // This property defines a reference to the volume VM program to be used by this particular render state
  pVolumeRS->SetVolumeBox0(volumeBox); // This property defines a reference to a volume dataset used by a VM program
  pVolumeRS->SetTransferFunction0(pTransferFunction); // This property defines a reference to a 1D transfer function used by a VM program

   // Create render state tree:
  Hue::ProxyLib::RenderStateTree *
    pRenderStateTree = pProject->RenderStateTrees().Create(); // This object is used to determine which render states should be used in what contexts during rendering

 Hue::ProxyLib::SetFrontfaceStateNode *
    pSetFrontfaceStateNode = pRenderStateTree->ChildNode().CreateSetFrontfaceStateNode(); // This node type defines which surface render state to use for the front faces when rendering a 3D shape object. All objects are by default not visible in the Render state tree, so this node must be used to instantiate an object.

  pSetFrontfaceStateNode->SetSurfaceStateSet(pVolumeRS); // When the active node type is Set front-/backface state, this parameter defines which surface render state to set
  pSetFrontfaceStateNode->SetShapeInstance(volumeSlice); // This property refers to the Shape instance object to which the render state should be applied

Hue::ProxyLib::SetFrontfaceStateNode *
    pSetFrontfaceStateNode2 = pSetFrontfaceStateNode->ChildNode().CreateSetFrontfaceStateNode(); // This node type defines which surface render state to use for the front faces when rendering a 3D shape object. All objects are by default not visible in the Render state tree, so this node must be used to instantiate an object.

  pSetFrontfaceStateNode2->SetSurfaceStateSet(pVolumeRS); // When the active node type is Set front-/backface state, this parameter defines which surface render state to set
  pSetFrontfaceStateNode2->SetShapeInstance(volumeSlice2); // This property refers to the Shape instance object to which the render state should be applied



  // Create view context (external):
  Hue::ProxyLib::ExternalViewContext *
	  pExternalViewContext = pViewer->ViewContextContainer3DOverlay()->ViewContexts().CreateExternalViewContext(); // This is an view context which provides a hook (through custom events) for performing direct OpenGL rendering from outside the component


  // Subscribe to render 3D event on external view context:
  pExternalViewContext->SubscribeEvent(&ExternalViewContext_RenderOpaqueExternalHandler);
  pExternalViewContext->SubscribeEvent(&ExternalViewContext_RenderTransparentExternalHandler);
  pExternalViewContext->SubscribeEvent(&ExternalViewContext_RenderDepthPeelExternalHandler);

  // Also handle trace 3D event:
  pExternalViewContext->SubscribeEvent(&ExternalViewContext_Trace3DExternalHandler);

  pExternalViewContext->SubscribeEvent(ExternalViewContext_DeInitExternalHandler);

  pExternalViewContext->SetEnableDepthPeeling(true);
  _pExternalViewContext = pExternalViewContext;

  InitSphere(1.0f);
  pExternalViewContext->SetGlobalBoundingBox(BoundingBox(DoubleVector3(-1.0, -1.0, -1.0), DoubleVector3(1.0, 1.0, 1.0)));

 

}

#define SPHERE_TESS_X 64
#define SPHERE_TESS_Y 32

static GLfloat SphereVertices[SPHERE_TESS_X * (SPHERE_TESS_Y + 1)][3];
static GLuint SphereIndices[SPHERE_TESS_X * SPHERE_TESS_Y][3 * 2];

static void InitSphere(float r)
{
  // prepare vertex-buffer
  for (int y = 0; y < SPHERE_TESS_Y + 1; ++y)
  {
    float phi = float(y) * (M_PI / SPHERE_TESS_Y);
    for (int x = 0; x < SPHERE_TESS_X; ++x)
    {
      float theta = float(x) * ((2 * M_PI) / SPHERE_TESS_X);
      SphereVertices[y * SPHERE_TESS_X + x][0] = r * cos(theta) * sin(phi);
      SphereVertices[y * SPHERE_TESS_X + x][1] = r * sin(theta) * sin(phi);
      SphereVertices[y * SPHERE_TESS_X + x][2] = r * cos(phi);
    }
  }

  // prepare index-buffer
  for (int y = 0; y < SPHERE_TESS_Y; ++y)
  {
    for (int x = 0; x < SPHERE_TESS_X; ++x)
    {
      // first triangle
      SphereIndices[y * SPHERE_TESS_X + x][0] = y * SPHERE_TESS_X + x;
      SphereIndices[y * SPHERE_TESS_X + x][1] = y * SPHERE_TESS_X + ((x + 1) % SPHERE_TESS_X);
      SphereIndices[y * SPHERE_TESS_X + x][2] = (y + 1) * SPHERE_TESS_X + ((x + 1) % SPHERE_TESS_X);

      // second triangle
      SphereIndices[y * SPHERE_TESS_X + x][3] = y * SPHERE_TESS_X + x;
      SphereIndices[y * SPHERE_TESS_X + x][4] = (y + 1) * SPHERE_TESS_X + ((x + 1) % SPHERE_TESS_X);
      SphereIndices[y * SPHERE_TESS_X + x][5] = (y + 1) * SPHERE_TESS_X + x;
    }
  }
}

static void DrawSphere()
{
  glVertexPointer(3, GL_FLOAT, 0, SphereVertices);
  glNormalPointer(GL_FLOAT, 0, SphereVertices);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glDrawElements(GL_TRIANGLES, SPHERE_TESS_X * SPHERE_TESS_Y * 3 * 2, GL_UNSIGNED_INT, SphereIndices);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
}

struct PerLightHashData
{
  QGLShaderProgram *prog;
  QGLShaderProgram *prog_depthpeel;
};

struct PerRenderContextData
{
  std::map<Hue::ProxyLib::uint64, PerLightHashData> perLightHashData;

  void DeleteShaders()
  {
    // delete PerRenderContextData
    std::map<Hue::ProxyLib::uint64, PerLightHashData>::iterator it;
    for (it = perLightHashData.begin(); it != perLightHashData.end(); ++it)
    {
      delete it->second.prog;
      delete it->second.prog_depthpeel;
    }
    perLightHashData.clear();
  }
};

std::map<Hue::ProxyLib::int32, PerRenderContextData> shaders;

QGLShaderProgram *CreateShaderProgram(const RenderExternalParameters *params, bool isDepthPeeling)
{
  QString vertexShaderSource =
"\
varying vec3 Position;\n\
varying vec3 Normal;\n\
void main()\n\
{\n\
  gl_Position = ftransform();\n\
  Position = (gl_ModelViewMatrix * gl_Vertex).xyz;\n\
  Normal = gl_NormalMatrix * gl_Normal;\n\
}\n\
";
  QGLShader vertexShader(QGLShader::Vertex);
  if (!vertexShader.compileSourceCode(vertexShaderSource))
  {
    fprintf(stderr, "failed to compile vertex shader:\n%s\n", vertexShader.log().toAscii());
    return NULL;
  }

  QString fragmentShaderSource = params->LightingShaderFunction.c_str();
  if (isDepthPeeling)
  {
    fragmentShaderSource.append(params->DepthPeelShaderDeclarations.c_str());
  }
  fragmentShaderSource.append(
"\
varying vec3 Position;\n\
varying vec3 Normal;\n\
void main()\n\
{\n\
");
  if (isDepthPeeling)
  {
    fragmentShaderSource.append("  DepthPeel();\n");
  }
  fragmentShaderSource.append(
"\
  vec3 N = normalize(Normal);\n\
  if (gl_FrontFacing) N = -N;\n\
  gl_FragColor = LightFragment(Position, N, vec4(1.0, 0.0, 0.0, 0.5));\n\
}\n\
");
  QGLShader fragmentShader(QGLShader::Fragment);
  if (!fragmentShader.compileSourceCode(fragmentShaderSource))
  {
    fprintf(stderr, "failed to compile fragment shader:\n%s\n", fragmentShader.log().toAscii());
    return NULL;
  }

  QGLShaderProgram *prog = new QGLShaderProgram(QGLContext::currentContext());
  prog->addShader(&vertexShader);
  prog->addShader(&fragmentShader);
  if (!prog->link())
  {
    fprintf(stderr, "failed to link shader program:\n%s\n", prog->log().toAscii());
    delete prog;
    return NULL;
  }

  if (isDepthPeeling)
  {
    prog->setUniformValue("DepthPeelingTexture", 7);
  }

  return prog;
}


QGLShaderProgram *GetShaderProgram(const RenderExternalParameters *params, bool isDepthPeeling)
{
  // make sure we have a PerRenderContextData struct for the
  std::map<Hue::ProxyLib::int32, PerRenderContextData>::iterator it = shaders.find(params->RenderContextIdentifier);
  if (it == shaders.end())
  {
    PerRenderContextData temp;
    it = shaders.insert(std::pair<Hue::ProxyLib::int32, PerRenderContextData>(params->RenderContextIdentifier, temp)).first;
  }

  std::map<Hue::ProxyLib::uint64, PerLightHashData>::iterator lightHashIterator;
  if ((lightHashIterator = it->second.perLightHashData.find(params->LightingShaderHash)) == it->second.perLightHashData.end())
  {
    PerLightHashData temp = { NULL, NULL };
    lightHashIterator = it->second.perLightHashData.insert(std::pair<Hue::ProxyLib::uint64, PerLightHashData>(params->LightingShaderHash, temp)).first;
  }

  if (!isDepthPeeling)
  {
    if (lightHashIterator->second.prog == NULL)
    {
      lightHashIterator->second.prog = CreateShaderProgram(params, false);
    }
    return lightHashIterator->second.prog;
  }
  else
  {
    if (lightHashIterator->second.prog_depthpeel == NULL)
    {
      lightHashIterator->second.prog_depthpeel = CreateShaderProgram(params, true);
    }
    return lightHashIterator->second.prog_depthpeel;
  }
}

void ExternalViewContext_DeInitExternalHandler(Hue::ProxyLib::DeInitExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID)
{
  // delete PerRenderContextData
  std::map<Hue::ProxyLib::int32, PerRenderContextData>::iterator it = shaders.find(context->RenderContextIdentifier());
  if (it != shaders.end())
  {
    it->second.DeleteShaders();
    shaders.erase(it);
  }
}

void RenderTransparent(const RenderExternalParameters *params, bool isDepthPeeling)
{
  QGLShaderProgram *prog = GetShaderProgram(params, isDepthPeeling);

//  glEnable(GL_CULL_FACE);
  glDisable(GL_CULL_FACE);
  if (isDepthPeeling)
  {
    prog->setUniformValue("ViewportScale", params->ViewportScale.X, params->ViewportScale.Y);
  }

  prog->bind();
  DrawSphere();
  prog->release();
}

void ExternalViewContext_RenderTransparentExternalHandler(Hue::ProxyLib::RenderTransparentExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID)
{
  RenderTransparent(&context->RenderParameters(), false);
}

void ExternalViewContext_RenderDepthPeelExternalHandler(Hue::ProxyLib::RenderDepthPeelExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID)
{
  RenderTransparent(&context->RenderParameters(), true);
}

void
ExternalViewContext_RenderOpaqueExternalHandler(Hue::ProxyLib::RenderOpaqueExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID)
{
  glDisable(GL_CULL_FACE);
  glColor3f(1.0f, 1.0f, 1.0f);

  glDisable(GL_LIGHTING);

  // draw a single gourard-shaded triangle
  glBegin(GL_TRIANGLES);
  glColor3f(1.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);

  glColor3f(0.0f, 1.0f, 0.0f);
  glVertex3f(1.0f, 0.0f, 0.0f);

  glColor3f(0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 1.0f, 0.0f);
  glEnd();
	
  if (ExternalViewContextWidget::strokeman)
	 ExternalViewContextWidget::strokeman->drawGL();

}

static bool IntersectLineWithTriangle(
  double *prWriteU, double *prWriteV, double *prWriteT,
  const DoubleVector3 &cStart, const DoubleVector3 &cDirection,
  const DoubleVector3 &cV0, const DoubleVector3 &cV1, const DoubleVector3 &cV2)
{
  DoubleVector3
    cEdge1,
    cEdge2;

  // Find vectors for two edges sharing vert0
  cEdge1 = cV1 - cV0;
  cEdge2 = cV2 - cV0;

  DoubleVector3
    cVec;

  // Begin calculating determinant - also used to calculate U parameter
  cVec.CrossProduct(cDirection, cEdge2);

  double
    rDet = cEdge1.DotProduct(cVec);

  // If determinant is near zero, ray lies in plane of triangle
  if (rDet > -0.0000001f && rDet < 0.0000001f)
  {
    return false;
  }

  double
    rInvDet = 1.0 / rDet;

  // Calculate distance from vert0 to ray origin
  DoubleVector3
    cTVec = cStart - cV0;

  double
    rU = cTVec.DotProduct(cVec) * rInvDet;

  // Calculate U parameter and test bounds
  if (rU < 0.0f || rU > 1.f)
  {
    return false;
  }

  // Prepare to test V parameter
  DoubleVector3
    cQVec;

  cQVec.CrossProduct(cTVec, cEdge1);

  double
    rV = cDirection.DotProduct(cQVec) * rInvDet;

  if (rV < 0.0f || (rU + rV) > 1.f)
  {
    return false;
  }

  double
    rT = cEdge2.DotProduct(cQVec) * rInvDet;

  if (rT < 0.f)
  {
    return false;
  }

  if (prWriteU) *prWriteU = rU;
  if (prWriteV) *prWriteV = rV;
  if (prWriteT) *prWriteT = rT;

  return true;
}

void
ExternalViewContext_Trace3DExternalHandler(Hue::ProxyLib::Trace3DExternalEventContext *context, Hue::ProxyLib::int64 subscriptionID)
{
  // prepare vertex-buffer
  for (int i = 0; i < SPHERE_TESS_X * SPHERE_TESS_Y; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      DoubleVector3
        v1 = DoubleVector3(
               SphereVertices[SphereIndices[i][j * 3 + 0]][0],
               SphereVertices[SphereIndices[i][j * 3 + 0]][1],
               SphereVertices[SphereIndices[i][j * 3 + 0]][2]),
        v2 = DoubleVector3(
               SphereVertices[SphereIndices[i][j * 3 + 1]][0],
               SphereVertices[SphereIndices[i][j * 3 + 1]][1],
               SphereVertices[SphereIndices[i][j * 3 + 1]][2]),
        v3 = DoubleVector3(
               SphereVertices[SphereIndices[i][j * 3 + 2]][0],
               SphereVertices[SphereIndices[i][j * 3 + 2]][1],
               SphereVertices[SphereIndices[i][j * 3 + 2]][2]);

      double u, v, t;
      if (IntersectLineWithTriangle(
        &u, &v, &t,
        context->StartPos(), context->Direction(),
        v1, v2, v3))
      {
        IntersectionExternal *intersection = context->IntersectionList().Create();
        intersection->SetGlobalPosition(context->StartPos() + context->Direction() * t);
        intersection->SetIntersectionType(Hue::ProxyLib::IntersectionType::Triangle);
        intersection->SetGlobalNormal(FloatVector3(0.0f, 0.0f, 1.0f));
        intersection->SetHitRGBAColor(FloatVector4(1.0f, 1.0f, 1.0f, 1.0f));
      }
    }
  }

  DoubleVector3
    v1 = DoubleVector3(0.0, 0.0, 0.0),
    v2 = DoubleVector3(1.0, 0.0, 0.0),
    v3 = DoubleVector3(0.0, 1.0, 0.0);
  double u, v, t;
  if (IntersectLineWithTriangle(
    &u, &v, &t,
    context->StartPos(), context->Direction(),
    v1, v2, v3))
  {
    IntersectionExternal *intersection = context->IntersectionList().Create();
    intersection->SetGlobalPosition(context->StartPos() + context->Direction() * t);
    intersection->SetIntersectionType(Hue::ProxyLib::IntersectionType::Triangle);
    intersection->SetGlobalNormal(FloatVector3(0.0f, 0.0f, 1.0f));
    intersection->SetHitRGBAColor(FloatVector4(1.0f, 1.0f, 1.0f, 1.0f));
  }
}

/////////////////////////////////////////////////////////////////////////////
// ExternalViewContextWidget destructor

ExternalViewContextWidget::~ExternalViewContextWidget()
{
  _pExternalViewContext->Delete();
}


void ExternalViewContextWidget::mouseMoveEvent(QMouseEvent *pMouseEvent)
{   if (sketching_on)
		strokeman->HandleMouseMotion(toWorld(pMouseEvent->x(),pMouseEvent->y()) );
	HueQtDefaultViewerWidget::mouseMoveEvent(pMouseEvent);
}

void ExternalViewContextWidget::mousePressEvent(QMouseEvent * pMouseEvent){

	Hue::Util::DoubleMatrix4x4
		cModelViewMatrix = pViewer->RenderLayer3DInsideMagnifier()->LastModelViewMatrix(),
	  cProjMatrix = pViewer->RenderLayer3DInsideMagnifier()->LastProjectionMatrix();

	//ispressed = true;
 
	if (sketching_on){
		strokeman->HandleMouseEvent(1, 1, toWorld(pMouseEvent->x(),pMouseEvent->y()) );
		//Trace2DEventContext * pcEventContext = this->Perform2DTrace(QPoint(pMouseEvent->x(),pMouseEvent->y()));

	}else
		HueQtDefaultViewerWidget::mousePressEvent(pMouseEvent);
}

void ExternalViewContextWidget::mouseReleaseEvent(QMouseEvent * pMouseEvent){
	//ispressed = false;
	//newstroke = true;
	//issketching = false;
	strokeman->HandleMouseEvent(1, 2, toWorld(pMouseEvent->x(),pMouseEvent->y()) );
	HueQtDefaultViewerWidget::mouseReleaseEvent(pMouseEvent);
}


void ExternalViewContextWidget::keyReleaseEvent(QKeyEvent * pQKeyEvent){

	if (pQKeyEvent->text()==" "){
		sketching_on = !sketching_on;		
	}

	/*if (pQKeyEvent->text()=="c")
		pConvexHullShape->SetInputShape((DirectShape *)m_pShapeInstance->Shape());*/

	HueQtDefaultViewerWidget::keyReleaseEvent(pQKeyEvent);
}

/*

// Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer2 = pViewer->ViewContextContainerOutsideMagnifier();

  pViewContextContainer2->SetName("Outside magnifier (left)");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer3 = pViewer->ViewContextContainer();

  pViewContextContainer3->SetName("Inside and outside magnifier (left)");

  // Create view context (world background):
  Hue::ProxyLib::WorldBackgroundViewContext *
    pWorldBackgroundViewContext = pViewContextContainer3->ViewContexts().CreateWorldBackgroundViewContext(); // This object instantiates a world background object

  // Create view context (shape instance):
  Hue::ProxyLib::ShapeInstanceViewContext *
    pShapeInstanceViewContext = pViewContextContainer3->ViewContexts().CreateShapeInstanceViewContext(); // This object instantiates a Shape instance object

  pShapeInstanceViewContext->SetEnableCollision(false); // This property defines whether object collision detection (i.e. trace line functionality) is enabled or not for this object

  // Create view context (external):
  Hue::ProxyLib::ExternalViewContext *
	  pExternalViewContext = pViewer->ViewContextContainer3DOverlay()->ViewContexts().CreateExternalViewContext(); // This is an view context which provides a hook (through custom events) for performing direct OpenGL rendering from outside the component

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer4 = pViewer->ViewContextContainerInsideMagnifierRight();

  pViewContextContainer4->SetName("Inside magnifier (right)");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer5 = pViewer->ViewContextContainerOutsideMagnifierRight();

  pViewContextContainer5->SetName("Outside magnifier (right)");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer6 = pViewer->ViewContextContainerRight();

  pViewContextContainer6->SetName("Inside and outside magnifier (right)");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer7 = pViewer->ViewContextContainerInsideAndOutsideMagnifier();

  pViewContextContainer7->SetName("Inside and outside magnifier");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer8 = pViewer->ViewContextContainerLightweightObjects();

  pViewContextContainer8->SetName("Inside and outside magnifier (lightweight)");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer9 = pViewer->ViewContextContainerSelectedObjects();

  pViewContextContainer9->SetName("Selected objects");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer10 = pViewer->ViewContextContainer3DOverlay();

  pViewContextContainer10->SetName("3D overlay");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer11 = pViewer->ViewContextContainerOrthogonalOverlay();

  pViewContextContainer11->SetName("Orthogonal overlay inside and outside magnifier");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer12 = pViewer->ViewContextContainerOrthogonalOverlayInsideMagnifier();

  pViewContextContainer12->SetName("Orthogonal overlay inside magnifier");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer13 = pViewer->ViewContextContainerOrthogonalOverlayOutsideMagnifier();

  pViewContextContainer13->SetName("Orthogonal overlay outside magnifier");

  // Create view context (container):
  Hue::ProxyLib::ViewContextContainer *
    pViewContextContainer14 = pViewer->ViewContextContainerControllers();

  pViewContextContainer14->SetName("Controllers");

  // Create marker shape:
  Hue::ProxyLib::MarkerShape *
    pMarkerShape = pViewer->AxisShape(); // This is a marker shape, i.e. a mesh object in the shape of arrow/crosshair/axis in one of a few available styles

  pMarkerShape->SetName("Axis"); //  (Default: "Marker shape")
  pMarkerShape->SetMemoryCachePolicy(Hue::ProxyLib::CachePolicy::NoCaching); // This parameter controls the policy used for caching produced data in memory (Default: KeepLastInMemoryWhileValid)
  pMarkerShape->SetStyle(Hue::ProxyLib::MarkerShapeStyle::Axis); // This property defines what kind of arrow to create (Default: RoundArrowFromCenter)

  // Create marker shape:
  Hue::ProxyLib::MarkerShape *
    pMarkerShape2 = pViewer->SunShape();

  pMarkerShape2->SetName("Sun"); //  (Default: "Marker shape")
  pMarkerShape2->SetMemoryCachePolicy(Hue::ProxyLib::CachePolicy::NoCaching); //  (Default: KeepLastInMemoryWhileValid)
  pMarkerShape2->SetStyle(Hue::ProxyLib::MarkerShapeStyle::Sun); //  (Default: RoundArrowFromCenter)

  // Create marker shape:
  Hue::ProxyLib::MarkerShape *
    pMarkerShape3 = pViewer->NorthArrowShape();

  pMarkerShape3->SetName("North arrow"); //  (Default: "Marker shape")
  pMarkerShape3->SetMemoryCachePolicy(Hue::ProxyLib::CachePolicy::NoCaching); //  (Default: KeepLastInMemoryWhileValid)
  pMarkerShape3->SetStyle(Hue::ProxyLib::MarkerShapeStyle::FlatArrowFromCenter); //  (Default: RoundArrowFromCenter)

  // Create orbit camera controller:
  Hue::ProxyLib::OrbitCameraController *
    pOrbitCameraController = pViewer->OrbitCameraController(); // This object is a Orbit camera controller, i.e. an object which provides Orbit behavior to a camera

  pOrbitCameraController->SetActive(true); // This property determines whether the control is active or not
  pOrbitCameraController->SetInertiaFallOff(9.0); // This property defines the speed at which rotation inertia will fall off to stand still, given in half-life time (seconds). (Default: 0.10000)
  pOrbitCameraController->SetPivotPointMode(Hue::ProxyLib::PivotPointMode::Centered); // This property maintains pivot point status. (Default: Floating)
  pOrbitCameraController->SetAutoZoomFactor(0.899f); // Upon rotating in to pivot point, this factor defines zoom in. 1: Keep current distance. Close to 0: Zoom fully in. (Default: 1.0f)
  pOrbitCameraController->SetCameraTarget(pCamera); // This property holds a reference to the camera which is controlled
  pOrbitCameraController->SetViewContextContainerController(pViewContextContainer3); // If the controller should have some 3D representation, this property must be set to indicate a View context container where it can place its required view contexts

      // With that in place, we can create the references to this object:
  pViewer->SetActiveCameraController(pOrbitCameraController); // The active camera controller; may only be set to one of the camera controllers owned by this viewer

  // Create fly camera controller:
  Hue::ProxyLib::FlyCameraController *
    pFlyCameraController = pViewer->FlyCameraController(); // This object is a fly camera controller, i.e. an object which provides fly behavior to a camera

  pFlyCameraController->SetInertiaFallOff(9.0); //  (Default: 0.10000)
  pFlyCameraController->SetCameraTarget(pCamera);
  pFlyCameraController->SetViewContextContainerController(pViewContextContainer14);

  // Create 3d object controller:
  Hue::ProxyLib::Edit3DController *
    pEdit3DController = pViewer->Edit3DController(); // This object is a generic 3D object controller, i.e. an object which provides move/rotate/scale behavior to any 3D object

  pEdit3DController->SetCameraTarget(pCamera);
  pEdit3DController->SetViewContextContainerController(pViewContextContainer14);
  pEdit3DController->SetViewContextContainerSelection(pViewContextContainer9); // If the Edit3D controller should have some 3D representation for selected objects, this property must be set to indicate a View context container where it can place its required view contexts

  // Create paint polyline controller:
  Hue::ProxyLib::PaintPolylineController *
    pPaintPolylineController = pViewer->PaintPolylineController(); // This object is a controller which provides behavior to paint polylines

  pPaintPolylineController->SetCameraTarget(pCamera);
  pPaintPolylineController->SetViewContextContainerController(pViewContextContainer11);

  // Create edit transfer function controller:
  Hue::ProxyLib::EditTransferFunctionController *
    pEditTransferFunctionController = pViewer->EditTransferFunctionController(); // This object is a controller which edits transfer functions

  pEditTransferFunctionController->SetCameraTarget(pCamera);
  pEditTransferFunctionController->SetViewContextContainerController(pViewContextContainer11);
  pEditTransferFunctionController->SetViewContextContainerSelection(pViewContextContainer9); // If the Edit Transfer Function controller should have some 3D representation for selected objects, this property must be set to indicate a View context container where it can place its required view contexts
  pEditTransferFunctionController->SetViewContextContainerPrimary(pViewContextContainer3); // The Edit Transfer Function places the created scene objects here.

  // Create view split controller:
  Hue::ProxyLib::ViewSplitController *
    pViewSplitController = pViewer->ViewSplitController(); // This object is a controller which manages interactive split view changes

  pViewSplitController->SetCameraTarget(pCamera);
  pViewSplitController->SetViewContextContainerController(pViewContextContainer11);
  pViewSplitController->SetViewContextContainerOuterSplitLine(pViewContextContainer13); // The ViewSplit controller places the outer split line graphics here.
  pViewSplitController->SetViewContextContainerInnerSplitLine(pViewContextContainer12); // The ViewSplit controller places the inner split line graphics here.


*/