
exed := ./

exe=main.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# C++
#################

cc := mpic++
copt := -std=c++17  -O3
#coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80,mem:managed
coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80,mem:separate
coptdebug := -g -traceback

# CPU environment for debug
#cc := icpx
#copt := -std=c++17  -O3 -lrt
#coptgpu := 
#coptdebug := -g -traceback

.SUFFIXES: .cpp .o

.cpp.o: 
	${cc} ${copt} ${coptgpu} ${coptdebug} -c $< 

OBJS = main.o hydro.o boundary.o mpi_routines.o

mpi_routines.o : mpi_routines.cpp
hydro.o: hydro.cpp hydro.hpp mpi_routines.o
boundary.o: boundary.cpp boundary.hpp mpi_routines.o
main.o: main.cpp hydro.o boundary.o mpi_routines.o

${exed}/${exe}: ${OBJS}
	${cc} ${copt} ${coptgpu} ${coptdebug} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} ${ana} *.o *.mod *~
