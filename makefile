objects := $(exe %.o)

CXX = g++

PROGS = optDriving

LDFLAGS = -lomp

CPPFLAGS = -Xpreprocessor -fopenmp #openMP flags (If parallelized with openMP)
MPICXX?= mpicxx #MPI Flags (If parallelized with MPI)
#BOOSTPATH = ../lib/boost_1_72_0/
BOOSTPATH = ../usr/local/boost_1_72_0/
INCLUDE = -I $(BOOSTPATH) -I ./include/
LIB = -L $(BOOSTPATH)
CPPFILES = ./src/main.cpp ./src/InitFunctions.cpp ./src/StationarySolver.cpp ./src/TimeDependentSolver.cpp ./src/CoordinateFunctions.cpp ./src/HelperFunctions.cpp ./src/InterpolationFunctions.cpp ./src/ProblemFunctions.cpp ./src/ProbabilityFunctions.cpp
OBJ = $(SRC:.cpp = .o)
EXAMPLE =

.PHONY: all clean

all : $(PROGS)

optDriving: $(CPPFILES)
	$(CXX) -std=c++17 -O3 -o $@ $^ $(INCLUDE)

run: optDriving
	./optDriving $(EXAMPLE)

cleanMain:
	rm -f *.o $(PROGS) *.txt *.dat

cleanOutput:
	rm ./output/*.txt
	rm ./output/*.dat
