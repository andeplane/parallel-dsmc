include(../defaults.pri)

CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

TEMPLATE = lib

TARGET = myapp

SOURCES += myclass.cpp \
    system.cpp \
    state.cpp \
    vec3.cpp \
    random.cpp \
    communication.cpp \
    grid.cpp \
    cell.cpp \
    particlemover.cpp
HEADERS += myclass.h \
    system.h \
    defines.h \
    state.h \
    vec3.h \
    random.h \
    communication.h \
    grid.h \
    cell.h \
    particlemover.h
