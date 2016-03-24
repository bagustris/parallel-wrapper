SHELL=/bin/sh
F90=gfortran
CC=gcc
OPENMP=-fopenmp
NUM=4
#OPT=-Wall -O0 $(OPENMP)
OPT=-Wall -O3 -fexpensive-optimizations -ffast-math -mpc64 -march=native $(OPENMP)
LIB=-lm

# required input and data files.
INPUTS= argon_108.inp argon_2916.inp argon_78732.inp \
	argon_108.rest argon_2916.rest argon_78732.rest

# baseline sources
BASE=ljmd-c.c ljmd-f.f90

ALL=ljmd-c.x ljmd-c0.x ljmd-c1.x ljmd-c2.x ljmd-c3.x ljmd-c4.x ljmd-cc.x \
	ljmd-f.x ljmd-f0.x ljmd-f1.x ljmd-f2.x ljmd-f3.x ljmd-f4.x ljmd-fc.x

all: $(ALL) input $(BASE)

.SUFFIXES:
.SUFFIXES: .c .f90 .x
.PHONY: all clean bench input

input: $(INPUTS)

$(INPUTS): 
	ln -s ../00_input/$@ $@

$(BASE): 
	ln -s ../01_baseline/$@ $@

bench: all $(INPUTS)
	@echo ========================================== >  ljmd.time 
	@for s in c c0 c1 c2 c3 c4 f f0 f1 f2 f3 f4 cc fc; do \
	  echo ========================================== >> ljmd.time ; \
	  echo ljmd-$$s.x 108 atoms                       >> ljmd.time ; \
	  ( env OMP_NUM_THREADS=4 time ./ljmd-$$s.x < argon_108.inp )  2>> ljmd.time ; \
	done
	@echo ========================================== >> ljmd.time 
	@for s in c c0 c1 c2 c3 c4 f f0 f1 f2 f3 f4 cc fc; do \
	  echo ========================================== >> ljmd.time ; \
	  echo ljmd-$$s.x 2916 atoms                      >> ljmd.time ; \
	  ( env OMP_NUM_THREADS=4 time ./ljmd-$$s.x < argon_2916.inp ) 2>> ljmd.time ; \
	done
	@echo ========================================== >> ljmd.time 

clean:
	-rm -f $(ALL) *.mod *~ fort.* argon_108.dat  argon_108.xyz	\
	lj_108.restart lj_108.xyz log.lj_108 argon_108.gmon			\
	argon_108.prof argon_2916.dat argon_2916.xyz log.lammps		\
	lj_2916.restart lj_2916.xyz log.lj_2916 $(INPUTS) $(BASE)

.c.x:
	$(CC) $(OPT) -o $@ $< $(LIB)

.f90.x:
	$(F90) $(OPT) -o $@ $<

# dependencies
ljmd-c.x:  ljmd-c.c
ljmd-c0.x: ljmd-c0.c
ljmd-c1.x: ljmd-c1.c
ljmd-c2.x: ljmd-c2.c
ljmd-c3.x: ljmd-c3.c
ljmd-c4.x: ljmd-c4.c

ljmd-f.x:  ljmd-f.f90
ljmd-f0.x: ljmd-f0.f90
ljmd-f1.x: ljmd-f1.f90
ljmd-f2.x: ljmd-f2.f90
ljmd-f3.x: ljmd-f3.f90
ljmd-f4.x: ljmd-f4.f90

ljmd-cc.x: ljmd-cc.c
ljmd-fc.x: ljmd-fc.f90

