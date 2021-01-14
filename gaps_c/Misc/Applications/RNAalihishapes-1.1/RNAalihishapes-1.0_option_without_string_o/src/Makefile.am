#libs
SUBDIRS = libs rtlib librna

bin_PROGRAMS =    RNAalihishapes
AM_CPPFLAGS = -I$(srcdir)/rtlib -I$(srcdir)/librna


RNAalihishapes_SOURCES = main.cc filter_add.hh \
                     mfe_answer_v4.hh hishape_filter.hh \
                     pfall_microstate.cc pfall_microstate.hh \
                     probs_h_microstate.cc probs_h_microstate.hh \
                     probs_hplus_microstate.cc probs_hplus_microstate.hh \
                     probs_m_microstate.cc probs_m_microstate.hh \
                     probs_a_microstate.cc probs_a_microstate.hh
RNAalihishapes_CPPFLAGS = $(BOOST_CPPFLAGS) $(AM_CPPFLAGS)
RNAalihishapes_LDFLAGS = $(BOOST_DATE_TIME_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
RNAalihishapes_LDADD = $(BOOST_DATE_TIME_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) rtlib/rtlib.la librna/librna.so


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
# 
# 
# HiPath_SOURCES = hipath.cc hited_util.hh hited_util.cc mfe_answer.hh \
#                      hishapeh_mfe_pp.cc hishapeh_mfe_pp.hh \
#                      hishapehplus_mfe_pp.cc hishapehplus_mfe_pp.hh \
#                      hishapem_mfe_pp.cc hishapem_mfe_pp.hh \
#                      hishapeb_mfe_pp.cc hishapeb_mfe_pp.hh
# HiPath_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(srcdir)/rtlib -I$(srcdir)/librna
# HiPath_LDFLAGS = $(BOOST_DATE_TIME_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
# HiPath_LDADD = $(BOOST_DATE_TIME_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(srcdir)/libs/libRNA.a librna/librna.la rtlib/rtlib.la comp/comp.la -lm
# 
# 
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


EXTRA_DIST = ShowTree ShowTree.opt.rb tree.rb LibChangeLog Extensions/alifold.hh    Extensions/bpfilter.hh  Extensions/mfesubopt.hh \
             Extensions/pknot_enforce.hh  Extensions/pknot_shape.hh  Extensions/probabilities.hh        Extensions/rnaoptions.hh \
             Extensions/shapes.hh      Extensions/typesRNAfolding.hh  Extensions/alignment.hh  Extensions/evalfold.hh  Extensions/outside.hh \
             Extensions/pknot_extras.hh   Extensions/pknot_stems.hh  Extensions/rnaoptions_defaults.hh  Extensions/rules.hh  Extensions/singlefold.hh
#                   rtlib/adp.hh          rtlib/bigint.hh   rtlib/erase.hh         rtlib/hashlist.hh     rtlib/output.hh     rtlib/ref.hh    rtlib/sample.hh      rtlib/string.cc       rtlib/table.hh \
#                   rtlib/algebra.hh      rtlib/bitops.hh   rtlib/filter.hh        rtlib/hash_stats.hh     rtlib/pool.hh       rtlib/rna.hh    rtlib/sequence.hh    rtlib/string.hh       rtlib/terminal.hh \
#                   rtlib/asymptotics.hh  rtlib/cm_alph.hh  rtlib/hashtng.hh     rtlib/map_pool.hh   rtlib/push_back.hh  rtlib/rope.hh   rtlib/shape_alph.hh  rtlib/string.lo       rtlib/treelist.hh \
#                   rtlib/backtrack.hh    rtlib/cstr.h      rtlib/generic_opts.hh  rtlib/list.hh        rtlib/move.hh       rtlib/range.hh      rtlib/shape.hh       rtlib/subopt.hh      rtlib/vector_sparse.hh \
#                   rtlib/bench.hh        rtlib/empty.hh    rtlib/hash.hh          rtlib/multipool.hh  rtlib/rational.hh   rtlib/rules.hh  rtlib/singleton.hh   rtlib/subsequence.hh

