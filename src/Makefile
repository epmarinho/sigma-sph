#!/usr/bin/make
# @<Makefile>::

# WARNING: leave a free line-feed between command blocks and a tail line-feed at the Makefile end-line

project-name	='Anisotropic Covariance-based SPH-TreeCode'
author		='Eraldo Pereira Marinho'
creation-date	='Fri Mar 16 15:44:24 BRT 2012'
local		='Rio Claro, SP'
country		='Brazil'
institute	='Sao Paulo State University, IGCE/UNESP, Department of Statistics, Applied Mathematics and Computing'
address		='Avenida 24A, 1515, Sala 21, Bloco I, DEMAC'
zipcode		='13.506-900'

# memory-model=-mcmodel=large
cpu-type=-march=native

# Uncomment next line if you wish to integrate thermal epecific energies according to the SPH adiabatic heating equation:
# directives=-D_INTEGRATE_THERMAL_ENERGIES_ -D_FIND_SELFCONSISTENT_KNN_CLUSTER_
# directives=-D_INTEGRATE_THERMAL_ENERGIES_
# directives=-D_FIND_SELFCONSISTENT_KNN_CLUSTER_
# directives=

#Optimization method (-g to generate debuger graph) -- do not uncomment:
# optmizing-flags=-g
optmizing-flags=-Ofast -floop-parallelize-all -floop-nest-optimize -fvect-cost-model=cheap -ftree-loop-if-convert
# optmizing-flags=-Ofast -floop-parallelize-all -floop-nest-optimize
# optmizing-flags=-Ofast -floop-parallelize-all
# optmizing-flags=-Ofast

include-dir=./include/

bin-dir=../bin/

openmp-flag=-fopenmp

CFLAGS=$(openmp-flag) $(optmizing-flags) $(directives) $(memory-model) $(cpu-type) -I$(include-dir) -Wall

LDFLAGS=-lm -lgomp

code-prefix=sph-treecode_

header-files=\
		$(include-dir)$(code-prefix)defaults.h\
		$(include-dir)$(code-prefix)eigenvalues.h

C-codes=\
		$(code-prefix)main.c\
		$(code-prefix)binsort.c\
		$(code-prefix)io.c\
		$(code-prefix)sigmaoctree.c\
		$(code-prefix)knn.c\
		$(code-prefix)kernel.c\
		$(code-prefix)eigenvalues.c\
		$(code-prefix)quantities.c\

executable-image=$(bin-dir)sigma-treecode$(model)

relocatable-list=\
		$(code-prefix)main.o\
		$(code-prefix)binsort.o\
		$(code-prefix)io.o\
		$(code-prefix)sigmaoctree.o\
		$(code-prefix)knn.o\
		$(code-prefix)kernel.o\
		$(code-prefix)eigenvalues.o\
		$(code-prefix)quantities.o\

all:
	+make _isotropic
	+make _anisotropic
	+make _isotropic_H_ATOM
	+make _anisotropic_H_ATOM

debug:
	+make clean
	+make _debug
_debug: model='_debug'
_debug: optmizing-flags=-g
_debug: $(relocatable-list)
	$(CC) -o $(executable-image) $^ $(LDFLAGS)


clean:
	$(RM) *.o

mostlyclean: clean
	$(RM) *~ *.tmp
	+(cd $(include-dir); $(RM) *~ *.tmp)

.PHONY: clean mostlyclean default indent

# fast: optmizing-flags=-Ofast -floop-parallelize-all -floop-nest-optimize -fvect-cost-model=unlimited -ftree-vectorize -ftree-loop-if-convert $(cpu-type)
fast: $(relocatable-list)
	$(CC) -o $(executable-image) $^ $(LDFLAGS)

mksphere:
	+(cd mkobjects;make clean;make $@)

mkdisk:
	+(cd mkobjects;make clean;make $@)

INDENTFLAGS=\
	-bad\
	-bap\
	-c0\
	-cd0\
	-cd2\
	-d0\
	-di16\
	-i4\
	-l132\
	-nbbb\
	-nbfda\
	-ncdb\
	-nfca\
	-npsl\
	-nsc\
	-nut\
	-orig\
	-ppi4\
	-sob\
	-ts4\

INDENT=/usr/bin/indent $(INDENTFLAGS)
indent: $(C-codes) $(header-files)
	$(INDENT) $(C-codes) $(header-files)


max_iterations_allowed=4
# both _isotropic and _anisotropic use as default the H_2 molecular gas model
_isotropic:
	+(make clean;make fast 'model=$@')
_anisotropic:
	+(make clean;\
	make fast 'directives=-D_FIND_SELFCONSISTENT_KNN_CLUSTER_ -DMAX_ITERATIONS_ALLOWED=$(max_iterations_allowed)' \
	'model=$@')
# H_ATOM directive switches the simulations to reproduce atomic hydrogen model
_isotropic_H_ATOM:
	+(make clean;make fast 'directives=-DH_ATOM' 'model=$@')
_anisotropic_H_ATOM:
	+(make clean;\
	make fast 'directives=-D_FIND_SELFCONSISTENT_KNN_CLUSTER_ -DMAX_ITERATIONS_ALLOWED=$(max_iterations_allowed) -DH_ATOM' \
	'model=$@')
_calc_integrals_only:
	+(make clean;make fast 'directives=-D__CALC_INTEGRALS_ONLY__' 'model=$@')
