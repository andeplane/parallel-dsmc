TEMPLATE = app
CONFIG += console
CONFIG -= qt

release {
    DEFINES +=
}

DEFINES +=

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    statisticssampler.cpp \
    random.cpp \
    grid.cpp \
    cutil.cpp \
    system.inc.cpp \
    unitconverter.cpp \
    settings.cpp \
    dsmc_io.cpp \
    dsmctimer.cpp \
    moleculemover.cpp \
    colliderbase.cpp \
    colliderspecular.cpp \
    colliderthermal.cpp \
    collidercercignanilampis.cpp

HEADERS += \
    system.h \
    cell.h \
    statisticssampler.h \
    random.h \
    grid.h \
    cutil.h \
    cinifile.h \
    unitconverter.h \
    settings.h \
    dsmc_io.h \
    dsmctimer.h \
    moleculemover.h \
    colliderbase.h \
    colliderspecular.h \
    colliderthermal.h \
    collidercercignanilampis.h


OTHER_FILES += \
    ../dsmc.ini

mac {
    CONFIG -= app_bundle
    LIBS   +=
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = icpc

    QMAKE_LFLAGS += -xCORE-AVX-I
    QMAKE_CXXFLAGS_RELEASE += -xCORE-AVX-I
}

unix:!mac {
    LIBS   +=
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = icpc
}

QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

QMAKE_LFLAGS += -ipo -no-prec-div -falign-functions=16
QMAKE_CXXFLAGS_RELEASE += -ipo -no-prec-div -falign-functions=16

QMAKE_LFLAGS -= -lm
QMAKE_LFLAGS -= -O2
QMAKE_LFLAGS += -O3
QMAKE_CFLAGS_RELEASE -= -fPIE
QMAKE_CXXFLAGS_RELEASE -= -fPIE
