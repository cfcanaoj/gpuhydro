
exe=kh.x

all: ${exe}

.PHONY: all clean allclean

##################################
# setup for CfCA XC50 begin
##################################
fc= ifort -mcmodel=medium
foptopenmp = -qopenmp -qopt-report=5
#foptopenmp = -parallel
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

# setup for CfCA XC50 end
##################################

##################################
# setup for Fugaku
##################################
#fc= mpifrtpx -Free
#foptopenmp = -Kopenmp
#fopt = -Kfast
# setup for Fugaku end
##################################

${exe}: main.f90
	${fc} ${fopt} ${foptopenmp} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

