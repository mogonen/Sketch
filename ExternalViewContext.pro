# Qmake project file for ExternalViewContext

# Assume a path if HUESPACE3_HOME is not defined
HUESPACE_SDK = $$(HUESPACE3_HOME)
isEmpty(HUESPACE_SDK) {
	HUESPACE_SDK = ../../../..
}

# Define LINUX if not on Windows
!win32 {
	DEFINES += LINUX
}

QT += opengl
TEMPLATE = app
DEPENDPATH += .
INCLUDEPATH += $$HUESPACE_SDK/include/HueSpace3/
INCLUDEPATH += $$HUESPACE_SDK/include/HueSpace3/Qt
LIBS += -L$$HUESPACE_SDK/lib/ -lhueproxy -lhueqtviewerwidget -luuid

HEADERS += ExternalViewContextWidget.h
SOURCES += ExternalViewContextWidget.cpp
SOURCES += main.cpp
TARGET = ExternalViewContext
