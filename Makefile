
# target might be WINDOWS, OS_X or LINUX
TARGET = OS_X

# Home directory of the acw program

# point to the inclue directory
INCLUDE = program/dsmc
SOURCEDIR = program/dsmc
INCLUDES = -I$(SOURCEDIR) -I/usr/local/Cellar/open-mpi/1.6.4/include
# library dir
LIBDIR = -L/usr/local/Cellar/open-mpi/1.6.4/lib 

# compiler specific flags
CFLAGS =  -D$(TARGET) -O3 $(INCLUDES) 

FFLAGS = $(LIBDIR) -lmpi_cxx -lmpi -lm

PROJECT = main


_obj = main.o cvector.o cell.o cmath.o colliderbase.o collidercercignanilampis.o collidermaxwell.o colliderspecular.o colliderthermal.o cutil.o dsmc_io.o dsmctimer.o grid.o moleculemover.o random.o settings.o statisticssampler.o system.o topology.o unitconverter.o

obj_main     = $(patsubst %,$(SOURCEDIR)/%, $(_obj))

CC 	= icpc

default: $(PROJECT)

$(PROJECT):  $(obj_main) 
	$(CC) $(INCLUDES) -o $(PROJECT) $(obj_main) $(FFLAGS)  
	
%.o: %.cpp
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   

clean:	
	rm program/dsmc/*.o

all: default
