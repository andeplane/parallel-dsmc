include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt
CONFIG   += c++11
TEMPLATE = lib

QMAKE_CXXFLAGS += -std=c++11

TARGET = myapp

SOURCES += myclass.cpp \
    system.cpp \
    state.cpp \
    vec3.cpp \
    random.cpp \
    communication.cpp \
    grid.cpp \
    cell.cpp \
    particlemover.cpp \
    settings.cpp \
    geometry.cpp \
    filemanager.cpp
HEADERS += myclass.h \
    system.h \
    defines.h \
    state.h \
    vec3.h \
    random.h \
    communication.h \
    grid.h \
    cell.h \
    particlemover.h \
    settings.h \
    geometry.h \
    filemanager.h
