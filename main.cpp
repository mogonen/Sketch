/***************************************************************************\

  main.cpp

\***************************************************************************/

// INCLUDES /////////////////////////////////////////////////////////////////

#include <QtGui/qapplication.h>
#include <QtGui/qlayout.h>
#include <QtGui/qframe.h>
#include <QtGui/qmainwindow.h>
#include "HueLicenseAPIHelper.h"

#include "ExternalViewContextWidget.h"

// FUNCTIONS ////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// main

int
main(int argc, char *argv[])
{
#if defined(_DEBUG) && defined(WIN32)
  _putenv("HUE_DEBUG=1");
#endif

  int
    windowWidth  = 800,
    windowHeight = 600;

  // Parse command line options to see if we should differ from the default window size:
  for (int arg = 1; arg < argc; arg++)
  {
    if (argv &&
        argv[arg] &&
        argv[arg][0] == '-')
    {
      if (argv[arg][1] == 'x')
      {
        windowWidth = atoi(&argv[arg][2]);
      }
      if (argv[arg][1] == 'y')
      {
        windowHeight = atoi(&argv[arg][2]);
      }
    }
  }

  
  QApplication
    qApplication(argc, argv);

  QMainWindow
    qMainWindow;

  QFrame
    *qFrame = new QFrame(&qMainWindow);

  QGridLayout
    *qGridLayout = new QGridLayout(qFrame);

  qGridLayout->setMargin(0);
  qGridLayout->addWidget(new ExternalViewContextWidget(qFrame, argc, argv));

  qFrame->setFocus();

  qMainWindow.setCentralWidget(qFrame);
  qMainWindow.resize(windowWidth, windowHeight);
  qMainWindow.show();

  qApplication.connect(&qApplication, SIGNAL(lastWindowClosed()), &qApplication, SLOT(quit()));

  int
    iReturnCode = qApplication.exec();

  
  return iReturnCode;
}
