# Makefile for Anisotropic, Covariance-based Smoothed Particle Hydrodynamics Project
# Author: Dr. Eraldo Pereira Marinho

author		 ='Eraldo Pereira Marinho'
project-name	 ='Anisotropic Covariance-based SPH-TreeCode'
to-do		 ='Multiphase self-gravitating gaseous system'
major-version	 ='6'
creation-date	 ='Fri Mar 16 15:44:24 BRT 2012'
last-modification='Mon Nov 11 13:38:00 BRT 2019'
local		 ='Rio Claro, SP'
country		 ='Brazil'
institute	 ='Sao Paulo State University, IGCE/UNESP, Department of Statistics, Applied Mathematics and Computing'
address		 ='Avenida 24A, 1515, Sala 21, Bloco I, DEMAC'
zipcode		 ='13.506-900'

MCMODEL=-mcmodel=medium

DIRECTIVES=-D_ADD_EFFECTIVE_NEIGHBORS_ -D_FULL_TEST_

OPTMIZATIONFLAGS=-g

CFLAGS=$(OPTMIZATIONFLAGS) $(DIRECTIVES) $(MCMODEL)

LDFLAGS=-lm

subdirs=src/

default:
	for f in $(subdirs); do (cd $$f; ${MAKE}); done

mostlyclean:
	for f in $(subdirs); do (cd $$f; ${MAKE} $@); done

clean:
	for f in $(subdirs); do (cd $$f; ${MAKE} $@); done
fast:
	+(cd src; make fast)

#indent-flags=-orig -nut -l131 -ncdb -nfc1 -nfca -bad -bap -bfda -nbfde -brs -ppi8 -npsl
#indent=/usr/bin/indent $(indent-flags)
indent:
	for f in $(subdirs); do (cd $$f; ${MAKE} $@); done

.PHONY: default clean mostlyclean indent $(subdirs)
