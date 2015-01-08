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
    unitconverter.cpp \
    settings.cpp \
    dsmc_io.cpp \
    dsmctimer.cpp \
    moleculemover.cpp \
    colliderbase.cpp \
    colliderspecular.cpp \
    colliderthermal.cpp \
    collidercercignanilampis.cpp \
    collidermaxwell.cpp \
    cvector.cpp \
    topology.cpp \
    solver.cpp \
    progressbar.cpp \
    perlin.cpp \
    complexgeometry.cpp \
    cmath.cpp \
    geometry.cpp \
    diamondsquare.cpp \
    statisticalproperty.cpp

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
    collidercercignanilampis.h \
    collidermaxwell.h \
    cvector.h \
    topology.h \
    solver.h \
    progressbar.h \
    perlin.h \
    complexgeometry.h \
    cmatrix.h \
    cmath.h \
    diamondsquare.h \
    statisticalproperty.h \
    statisticalvalue.h \
    defines.h


OTHER_FILES += \
    ../dsmc.ini

mac {
    CONFIG -= app_bundle
    INCLUDEPATH += /usr/local/include/
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = mpic++

    QMAKE_LFLAGS +=
    QMAKE_CXXFLAGS_RELEASE += -fast
    QMAKE_CXXFLAGS_RELEASE -= -fPIE
}

unix:!mac {
    LIBS   +=
    INCLUDEPATH +=
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS -O3
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = mpic++
}

QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$QMAKE_CXXFLAGS
