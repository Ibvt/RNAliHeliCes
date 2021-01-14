
##bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
##		RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
##		RNAalifold RNAplfold RNAup RNAaliduplex
#bin_PROGRAMS = RNAalifold RNAalieval Pathfinder
#
#LDADD = ../lib/libRNA.a -lm
#
#INCLUDES = -I$(srcdir)/../H





if WITH_LIBSVM
  RNALfold_LINK = ${CXXLD} \
                  $(AM_CXXFLAGS) $(CXXFLAGS) \
                  $(AM_LDFLAGS) $(LDFLAGS) -o $@
else
  RNALfold_LINK = ${CCLD} \
                  $(AM_CFLAGS) $(CFLAGS) \
                  $(AM_LDFLAGS) $(LDFLAGS) -o $@
endif


#@SUBDIRS = libs librna
#RNA2Dfold 
bin_PROGRAMS = RNAalifold RNAalieval Pathfinder 
#@RNAalihipath

LDADD = $(OPENMP_CFLAGS) ../lib/libRNA.a -lm

AM_CFLAGS = $(OPENMP_CFLAGS)

AM_CPPFLAGS = -I$(srcdir)/../H -I$(srcdir)/../lib 
#@ -I$(srcdir)/rtlib -I$(srcdir)/librna

#              RNA2Dfold_cmdl.c RNA2Dfold_cmdl.h
GENGETOPT_CMDL = RNAalifold_cmdl.c RNAalifold_cmdl.h

#RNA2Dfold.ggo
GENGETOPT_FILES =  RNAalifold.ggo 

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)


RNAalifold.c:    RNAalifold_cmdl.h

RNAalifold_cmdl.h  RNAalifold_cmdl.c: RNAalifold.ggo
	gengetopt -i RNAalifold.ggo

RNAalifold_SOURCES=RNAalifold_cmdl.c RNAalifold.c

# RNA2Dfold.c:    RNA2Dfold_cmdl.h
# 
# RNA2Dfold_cmdl.h  RNA2Dfold_cmdl.c: RNA2Dfold.ggo
# 	gengetopt -i RNA2Dfold.ggo
# 
# RNA2Dfold_SOURCES=RNA2Dfold_cmdl.c RNA2Dfold.c




# RNAalifold.c:    RNAalifold_cmdl.h
# 
# RNAalifold_cmdl.h  RNAalifold_cmdl.c: RNAalifold.ggo
# 	gengetopt -i RNAalifold.ggo
# 
# RNAalifold_SOURCES=RNAalifold_cmdl.c RNAalifold.c \
#                    eval_microstate.cc eval_microstate.hh rtlib/string.cc rtlib/string.hh
# RNAalihishapes_CPPFLAGS = $(AM_CPPFLAGS)
# RNAalihishapes_LDFLAGS = 
# RNAalihishapes_LDADD = librna/librna.la
# 
# EXTRA_DIST = ShowTree ShowTree.opt.rb tree.rb LibChangeLog Extensions/alifold.hh    Extensions/bpfilter.hh  Extensions/mfesubopt.hh \
#              Extensions/pknot_enforce.hh  Extensions/pknot_shape.hh  Extensions/probabilities.hh        Extensions/rnaoptions.hh \
#              Extensions/shapes.hh      Extensions/typesRNAfolding.hh  Extensions/alignment.hh  Extensions/evalfold.hh  Extensions/outside.hh \
#              Extensions/pknot_extras.hh   Extensions/pknot_stems.hh  Extensions/rnaoptions_defaults.hh  Extensions/rules.hh  Extensions/singlefold.hh \
#                    rtlib/adp.hh          rtlib/bigint.hh   rtlib/erase.hh         rtlib/hashlist.hh     rtlib/output.hh     rtlib/ref.hh    rtlib/sample.hh            rtlib/table.hh \
#                    rtlib/algebra.hh      rtlib/bitops.hh   rtlib/filter.hh        rtlib/hash_stats.hh     rtlib/pool.hh       rtlib/rna.hh    rtlib/sequence.hh        rtlib/terminal.hh \
#                    rtlib/asymptotics.hh  rtlib/cm_alph.hh  rtlib/hashtng.hh     rtlib/map_pool.hh   rtlib/push_back.hh  rtlib/rope.hh   rtlib/shape_alph.hh        rtlib/treelist.hh \
#                    rtlib/backtrack.hh    rtlib/cstr.h      rtlib/generic_opts.hh  rtlib/list.hh        rtlib/move.hh       rtlib/range.hh      rtlib/shape.hh       rtlib/subopt.hh      rtlib/vector_sparse.hh \
#                    rtlib/bench.hh        rtlib/empty.hh    rtlib/hash.hh          rtlib/multipool.hh  rtlib/rational.hh   rtlib/rules.hh  rtlib/singleton.hh   rtlib/subsequence.hh
