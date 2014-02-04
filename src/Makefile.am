lib_libdir = $(top_builddir)/lib
gnulib_libdir = $(top_builddir)/gnulib
AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_srcdir)/gnulib -Wall
AM_LDFLAGS = -L$(top_builddir)/lib -L$(top_builddir)/gnulib

if HAVE_SSE
AM_CPPFLAGS += -msse
#AM_CPPFLAGS += -mfpmath=sse
endif

if HAVE_SSE2
AM_CPPFLAGS += -msse2
endif

if HAVE_SSE3
AM_CPPFLAGS += -msse3
endif

if HAVE_SSSE3
AM_CPPFLAGS += -mssse3
endif

if HAVE_SSE41
AM_CPPFLAGS += -msse4.1
endif

AM_CPPFLAGS += $(OPENMP_CFLAGS)


bin_PROGRAMS =  generate_tiltangles \
								convex_comb         \
								ai_single_axis      \
								forward_op          \
								forward_op_atonce   \
								landweber_single_axis

# CAUTION: order of libraries matters! Order from "most dependent" to "least dependent"
# (or list each library twice)

generate_tiltangles_SOURCES = $(srcdir)/generate_tiltangles.c
generate_tiltangles_LDADD = $(lib_libdir)/libtilt.a \
														$(lib_libdir)/libutil.a \
														$(lib_libdir)/libCException.a \
														$(gnulib_libdir)/libgnu.a 

convex_comb_SOURCES = $(srcdir)/convex_comb.c
convex_comb_LDADD = $(lib_libdir)/libgfunc3.a \
										$(lib_libdir)/libutil.a \
										$(lib_libdir)/libCException.a \
										$(gnulib_libdir)/libgnu.a

ai_single_axis_SOURCES = $(srcdir)/ai_single_axis.c $(top_srcdir)/lib/ETreco.h
ai_single_axis_LDADD = $(lib_libdir)/libai.a \
											 $(lib_libdir)/libet.a \
											 $(lib_libdir)/libgfunc3.a \
											 $(lib_libdir)/libtilt.a \
											 $(lib_libdir)/libutil.a \
											 $(lib_libdir)/libiniparser.a \
											 $(lib_libdir)/libCException.a \
											 $(gnulib_libdir)/libgnu.a 
												

forward_op_SOURCES = $(srcdir)/forward_op.c $(top_srcdir)/lib/ETreco.h
forward_op_LDADD = $(lib_libdir)/libfwd_op.a \
									 $(lib_libdir)/libet.a \
					 				 $(lib_libdir)/libgfunc3.a \
									 $(lib_libdir)/libtilt.a \
									 $(lib_libdir)/libutil.a \
									 $(lib_libdir)/libiniparser.a \
									 $(lib_libdir)/libCException.a \
									 $(gnulib_libdir)/libgnu.a 

forward_op_atonce_SOURCES = $(srcdir)/forward_op_atonce.c $(top_srcdir)/lib/ETreco.h
forward_op_atonce_LDADD = $(lib_libdir)/libfwd_op.a \
													$(lib_libdir)/libet.a \
													$(lib_libdir)/libgfunc3.a \
													$(lib_libdir)/libtilt.a \
													$(lib_libdir)/libutil.a \
													$(lib_libdir)/libiniparser.a \
													$(lib_libdir)/libCException.a \
													$(gnulib_libdir)/libgnu.a 

landweber_single_axis_SOURCES = $(srcdir)/landweber_single_axis.c $(top_srcdir)/lib/ETreco.h
landweber_single_axis_LDADD = $(lib_libdir)/liblandw.a \
															$(lib_libdir)/libet.a \
															$(lib_libdir)/libgfunc3.a \
															$(lib_libdir)/libtilt.a \
															$(lib_libdir)/libutil.a \
															$(lib_libdir)/libiniparser.a \
															$(lib_libdir)/libCException.a \
															$(gnulib_libdir)/libgnu.a 

