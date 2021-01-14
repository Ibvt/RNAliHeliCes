#all: compile
#
#compile:
#	gcc -c -o librna/rnalib.pio -fpic   -std=c99 -O3 -DNDEBUG librna/rnalib.c
#	gcc -c -o librna/vienna/energy_par.pio -fpic   -Wall -Wno-char-subscripts -Wmissing-prototypes -Wtype-limits -g -std=c99  librna/vienna/energy_par.c
#	gcc -c -o librna/vienna/fold_vars.pio -fpic   -Wall -Wno-char-subscripts -Wmissing-prototypes -Wtype-limits -g -std=c99  librna/vienna/fold_vars.c
#	gcc -c -o librna/vienna/read_epars.pio -fpic   -Wall -Wno-char-subscripts -Wmissing-prototypes -Wtype-limits -g -std=c99  librna/vienna/read_epars.c
#	gcc -c -o librna/vienna/params.pio -fpic   -Wall -Wno-char-subscripts -Wmissing-prototypes -Wtype-limits -g -std=c99  librna/vienna/params.c
#	gcc -c -o librna/vienna/utils.pio -fpic   -Wall -Wno-char-subscripts -Wmissing-prototypes -Wtype-limits -g -std=c99  librna/vienna/utils.c
#	gcc -shared librna/rnalib.pio librna/vienna/energy_par.pio librna/vienna/fold_vars.pio librna/vienna/read_epars.pio librna/vienna/params.pio librna/vienna/utils.pio -lm -o librna/librna.so
#	gcc -c -o librna/rnalib.fpio -fpic   -std=c99 -O3 -DNDEBUG -ffast-math librna/rnalib.c
#	gcc -shared librna/rnalib.fpio librna/vienna/energy_par.pio librna/vienna/fold_vars.pio librna/vienna/read_epars.pio librna/vienna/params.pio librna/vienna/utils.pio -lm -o librna/librnafast.so
#
#clean:
#	rm -f librna/librna.so librna/librnafast.so librna/rnalib.fpio librna/librna.a librna/rnalib.o librna/vienna/energy_par.o librna/vienna/fold_vars.o librna/vienna/read_epars.o librna/vienna/read_epars.o librna/vienna/params.o librna/vienna/utils.o librna/rnalib.pio librna/vienna/energy_par.pio librna/vienna/fold_vars.pio librna/vienna/read_epars.pio librna/vienna/read_epars.pio librna/vienna/params.pio librna/vienna/utils.pio


lib_LTLIBRARIES = librna.la
librna_la_SOURCES = rnalib.c rnalib.h vienna/config.h vienna/dna_mathews1999.par  vienna/energy_par.c  vienna/fold_vars.h \
                    vienna/params.h      vienna/rna_andronescu2007.par  vienna/utils.c \
            vienna/dna_mathews2004.par  vienna/energy_par.h  vienna/intloops.h   vienna/pair_mat.h  \
            vienna/read_epars.c  vienna/rna_turner1999.par      vienna/utils.h vienna/data_structures.h  vienna/energy_const.h \
            vienna/fold_vars.c   vienna/params.c    vienna/read_epars.h  vienna/rna_turner2004.par
#            vienna/findpath.h vienna/findpath.c vienna/fold.h vienna/fold.c circfold.inc
#-O2 -DTURNER1999
librna_la_CFLAGS = --pedantic -Wall -std=c99
librna_la_LDFLAGS = -module -avoid-version
EXTRA_DIST = readme vienna/COPYING


