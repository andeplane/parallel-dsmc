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
    visualizer.cpp \
    testshader.cpp \
    solver.cpp \
    progressbar.cpp \
    perlin.cpp \
    oglshader.cpp \
    moviedata.cpp \
    movie.cpp \
    mesh.cpp \
    marchingcubes.cpp \
    ctexture.cpp \
    cshaders.cpp \
    copengl.cpp \
    complexgeometry.cpp \
    combined.cpp \
    cmath.cpp \
    citmap.cpp \
    cisosurface.cpp \
    cbitmap.cpp \
    camera.cpp \
    geometry.cpp \
    diamondsquare.cpp \
    statisticalproperty.cpp \
    lodepng.cpp

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
    fpsmanager.hpp \
    visualizer.h \
    testshader.h \
    solver.h \
    progressbar.h \
    perlin.h \
    oglshader.h \
    moviedata.h \
    mesh.h \
    marchingcubes.h \
    ctexture.h \
    cshaders.h \
    copengl.h \
    complexgeometry.h \
    cmatrix.h \
    cmath.h \
    cisosurface.h \
    cbitmap.h \
    camera.h \
    diamondsquare.h \
    statisticalproperty.h \
    statisticalvalue.h \
    defines.h \
    lodepng.h


OTHER_FILES += \
    ../dsmc.ini

mac {
    CONFIG -= app_bundle
    INCLUDEPATH += /usr/local/include/
    QMAKE_CXXFLAGS +=
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
    QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS
    QMAKE_CXX = icpc

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
