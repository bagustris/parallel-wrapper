SHELL=/bin/sh

F90=f95
CC=gcc
OPT=-Wall -O0
#OPT=-Wall -O 
#OPT=-Wall -O2 
#OPT=-Wall -O2 -funroll-loops 
#OPT=-Wall -O3
#OPT=-Wall -O3 -ffast-math
#OPT=-Wall -O3 -ffast-math -mpc64 -march=native
#OPT=-Wall -O3 -fexpensive-optimizations -ffast-math -mpc64 -march=native
LIB=-lm

# required input and data files.
INPUTS= argon_108.inp argon_2916.inp argon_78732.inp \
	argon_108.rest argon_2916.rest argon_78732.rest 

ALL=ljmd-c.x ljmd-f.x

all: $(ALL) input

.SUFFIXES:
.PHONY: all clean bench input

input: $(INPUTS)

$(INPUTS): 
	ln -s ../00_input/$@ $@

bench: all $(INPUTS)
	@echo ========================================== >  ljmd.time 
	@echo ljmd-c.x 108 atoms                         >> ljmd.time
	( time ./ljmd-c.x < argon_108.inp )             2>> ljmd.time
	@echo ========================================== >> ljmd.time 
	@echo ljmd-f.x 108 atoms                         >> ljmd.time
	( time ./ljmd-f.x < argon_108.inp )             2>> ljmd.time
	@echo ========================================== >> ljmd.time 

clean:
	-rm -f $(ALL) *.mod *~ fort.* argon_108.dat  argon_108.xyz	\
	lj_108.restart lj_108.xyz log.lj_108 argon_108-f.gmon		\
	argon_108.prof argon_2916.dat argon_2916.xyz log.lammps		\
	lj_2916.restart lj_2916.xyz log.lj_2916 argon_108-c.gmon	\
	$(INPUTS)

ljmd-f.x: ljmd-f.f90
	$(F90) $(OPT) -o $@ $<

ljmd-c.x: ljmd-c.c
	$(CC) $(OPT) -o $@ $< $(LIB)

