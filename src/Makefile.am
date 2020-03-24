#libs
SUBDIRS = H lib librna

#   eval_microstate_orig eval_microstate
bin_PROGRAMS = RNAalihishapes RNAalihipath 
AM_CPPFLAGS = -I$(srcdir)/rtlib -I$(srcdir)/librna


#filter_add.hh mfe_answer_v4.hh hishape_filter.hh
RNAalihishapes_SOURCES = rnaalihishapes.cc \
                     pfall_microstate.cc pfall_microstate.hh \
                     rep_consensus.cc rep_consensus.hh \
                     probs_h_microstate.cc probs_h_microstate.hh \
                     probs_hplus_microstate.cc probs_hplus_microstate.hh \
                     probs_m_microstate.cc probs_m_microstate.hh \
                     probs_a_microstate.cc probs_a_microstate.hh \
                     rtlib/string.cc rtlib/string.hh
# line1: CFLAGS for includes, line2: LDFLAGS shared library, line3: LIBS static library
RNAalihishapes_CPPFLAGS = $(BOOST_CFLAGS) $(AM_CPPFLAGS)
RNAalihishapes_LDFLAGS = #$(BOOST_LIBS)
RNAalihishapes_LDADD = -lboost_date_time -lboost_program_options $(srcdir)/lib/libRNA.a librna/librna.la 


# RNAHeliCes_SOURCES = main.cc filter_add.hh \
#                      mfe_answer_v4.hh hishape_filter.hh \
#                      hishapeh_mfeV4_pp.cc hishapeh_mfeV4_pp.hh \
#                      hishapehplus_mfeV4_pp.cc hishapehplus_mfeV4_pp.hh \
#                      hishapem_mfeV4_pp.cc hishapem_mfeV4_pp.hh \
#                      hishapeb_mfeV4_pp.cc hishapeb_mfeV4_pp.hh \
#                      pf_answer.hh pf_filter.hh mfe_answer.hh p_func.cc p_func.hh \
#                      hishapeh_mfe_pfx.cc hishapeh_mfe_pfx.hh \
#                      hishapehplus_mfe_pfx.cc hishapehplus_mfe_pfx.hh \
#                      hishapem_mfe_pfx.cc hishapem_mfe_pfx.hh \
#                      hishapeb_mfe_pfx.cc hishapeb_mfe_pfx.hh
# 
# RNAHeliCes_CPPFLAGS = $(BOOST_CPPFLAGS) $(AM_CPPFLAGS)
# RNAHeliCes_LDFLAGS = $(BOOST_DATE_TIME_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
# RNAHeliCes_LDADD = $(BOOST_DATE_TIME_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) librna/librna.la rtlib/rtlib.la comp/comp.la
# 

#hited_util.hh hited_util.cc mfe_answer.hh
RNAalihipath_SOURCES = rnaalihipath.cc \
                     eval_microstate.cc eval_microstate.hh \
                     probs_h_microstate.cc probs_h_microstate.hh \
                     probs_hplus_microstate.cc probs_hplus_microstate.hh \
                     probs_m_microstate.cc probs_m_microstate.hh \
                     probs_a_microstate.cc probs_a_microstate.hh \
                     rtlib/string.cc rtlib/string.hh
#-I$(srcdir)/rtlib -I$(srcdir)/librna
RNAalihipath_CPPFLAGS = $(BOOST_CFLAGS) $(AM_CPPFLAGS) -DWITH_RNAOPTIONS
RNAalihipath_LDFLAGS = #$(BOOST_LIBS)
RNAalihipath_LDADD = -lboost_date_time -lboost_program_options $(srcdir)/lib/libRNA.a librna/librna.la -lm


# eval_microstate_SOURCES = eval_microstate_main.cc \
#                      eval_microstate.cc eval_microstate.hh \
#                      rtlib/string.cc rtlib/string.hh
# #$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -DWITH_RNAOPTIONS $< -c -o $@ \
# #         && $(SED) -e 's/[^ ]\+boost[^ \n]\+//' $*.d > _t && mv _t $*.d
# eval_microstate_CPPFLAGS = $(AM_CPPFLAGS) -DWITH_RNAOPTIONS
# eval_microstate_LDFLAGS = 
# eval_microstate_LDADD = librna/librna.la -lm


# eval_microstate_orig_SOURCES = eval_microstate_main_orig.cc \
#                      eval_microstate_orig.cc eval_microstate_orig.hh \
#                      rtlib/string.cc rtlib/string.hh
# eval_microstate_orig_CPPFLAGS = $(AM_CPPFLAGS) -DWITH_RNAOPTIONS
# eval_microstate_orig_LDFLAGS = 
# eval_microstate_orig_LDADD = librna/librna.la -lm


# 
# 
# BFSPath_SOURCES = bfspath.cc mfe_answer.hh \
#                      hishapeh_mfe_pp.cc hishapeh_mfe_pp.hh \
#                      hishapehplus_mfe_pp.cc hishapehplus_mfe_pp.hh \
#                      hishapem_mfe_pp.cc hishapem_mfe_pp.hh \
#                      hishapeb_mfe_pp.cc hishapeb_mfe_pp.hh
# BFSPath_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(srcdir)/rtlib -I$(srcdir)/librna
# BFSPath_LDFLAGS = $(BOOST_DATE_TIME_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
# BFSPath_LDADD = $(BOOST_DATE_TIME_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(srcdir)/libs/libRNA.a librna/librna.la rtlib/rtlib.la


#      Extensions/shapes.hh  
EXTRA_DIST = ShowTree ShowTree.opt.rb tree.rb LibChangeLog Extensions/alifold.hh    Extensions/bpfilter.hh  Extensions/mfesubopt.hh \
             Extensions/pknot_enforce.hh  Extensions/pknot_shape.hh  Extensions/probabilities.hh        Extensions/rnaoptions.hh \
                   Extensions/typesRNAfolding.hh  Extensions/alignment.hh  Extensions/evalfold.hh  Extensions/outside.hh \
             Extensions/pknot_extras.hh   Extensions/pknot_stems.hh  Extensions/rnaoptions_defaults.hh  Extensions/rules.hh  Extensions/singlefold.hh \
                   rtlib/adp.hh          rtlib/bigint.hh   rtlib/erase.hh         rtlib/hashlist.hh     rtlib/output.hh     rtlib/ref.hh    rtlib/sample.hh            rtlib/table.hh \
                   rtlib/algebra.hh      rtlib/bitops.hh   rtlib/filter.hh        rtlib/hash_stats.hh     rtlib/pool.hh       rtlib/rna.hh    rtlib/sequence.hh        rtlib/terminal.hh \
                   rtlib/asymptotics.hh  rtlib/cm_alph.hh  rtlib/hashtng.hh     rtlib/map_pool.hh   rtlib/push_back.hh  rtlib/rope.hh   rtlib/shape.hh  rtlib/shape_alph.hh        rtlib/treelist.hh \
                   rtlib/backtrack.hh    rtlib/cstr.h      rtlib/generic_opts.hh  rtlib/list.hh        rtlib/move.hh       rtlib/range.hh     rtlib/subopt.hh      rtlib/vector_sparse.hh \
                   rtlib/bench.hh        rtlib/empty.hh    rtlib/hash.hh          rtlib/multipool.hh  rtlib/rational.hh   rtlib/rules.hh  rtlib/singleton.hh   rtlib/subsequence.hh
#             librna/rnalib.c librna/rnalib.h librna/vienna/config.h \
#             librna/vienna/energy_par.c  librna/vienna/fold_vars.h  librna/vienna/intl21dH.h  librna/vienna/intl22.h    librna/vienna/params.h      librna/vienna/rna_andronescu2007.par  librna/vienna/utils.c \
#             librna/vienna/dna_mathews2004.par  librna/vienna/energy_par.h  librna/vienna/intl11dH.h   librna/vienna/intl21.h    librna/vienna/pair_mat.h  librna/vienna/read_epars.c  librna/vienna/rna_turner1999.par      librna/vienna/utils.h \
#             librna/vienna/data_structures.h  librna/vienna/energy_const.h       librna/vienna/fold_vars.c   librna/vienna/intl11.h     librna/vienna/intl22dH.h  librna/vienna/params.c    librna/vienna/read_epars.h  librna/vienna/rna_turner2004.par \
#                    librna/readme librna/vienna/COPYING

