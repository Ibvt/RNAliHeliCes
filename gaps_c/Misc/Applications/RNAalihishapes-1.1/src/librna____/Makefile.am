
# CFLAGS=-Wall -g -std=c99
# LDLIBS=-lm
# 
# all: test
# 
# 
# test: test.o rnalib.o


lib_LTLIBRARIES = librna.la
librna_la_SOURCES = rnalib.c rnalib.h vienna/config.h \
            vienna/energy_par.c  vienna/fold_vars.h  vienna/intl21dH.h  vienna/intl22.h    vienna/params.h      vienna/rna_andronescu2007.par  vienna/utils.c \
            vienna/dna_mathews2004.par  vienna/energy_par.h  vienna/intl11dH.h   vienna/intl21.h    vienna/pair_mat.h  vienna/read_epars.c  vienna/rna_turner1999.par      vienna/utils.h \
            vienna/data_structures.h  vienna/energy_const.h       vienna/fold_vars.c   vienna/intl11.h     vienna/intl22dH.h  vienna/params.c    vienna/read_epars.h  vienna/rna_turner2004.par

librna_la_CFLAGS = --pedantic -Wall -g -std=c99
librna_la_LDFLAGS = -module -avoid-version
#librna_la_LDADD = -lm
EXTRA_DIST = readme vienna/COPYING


