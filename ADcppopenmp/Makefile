
exed := ./

exe=main.x

all: ${exed}/${exe}

.PHONY: all clean

#################
# Fortran
#################
cc := nvc++ 
copt := -std=c++17  -O3
coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80,mem:managed
#coptgpu := -mp=gpu -Minfo=mp,accel -gpu=cc80,mem:separate
coptdebug := -g -traceback

#fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

.SUFFIXES: .cpp .o

.cpp.o: 
	${cc} ${copt} ${coptgpu} ${coptdebug} -c $< 

OBJS = main.o hydro.o boundary.o

hydro.o: hydro.cpp hydro.hpp
boundary.o: boundary.cpp boundary.hpp
main.o: main.cpp hydro.o boundary.o

${exed}/${exe}: ${OBJS}
	${cc} ${copt} ${coptgpu} ${coptdebug} ${OBJS} -o $@

#################
# clean up
#################
clean:
	rm -f ${exed}/${exe} ${ana} *.o *.mod *~
