lib_srcdir = $(top_srcdir)/lib
gnulib_srcdir = $(top_srcdir)/gnulib
AM_CPPFLAGS = -I$(top_srcdir)/gnulib -Wall
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

noinst_LIBRARIES = libCException.a libiniparser.a libutil.a libgfunc3.a libetai.a

libCException_a_SOURCES = $(lib_srcdir)/CException.c $(lib_srcdir)/CException.h \
													$(lib_srcdir)/CExceptionConfig.h

libiniparser_a_SOURCES = $(lib_srcdir)/dictionary.c  	$(lib_srcdir)/dictionary.h \
												 $(lib_srcdir)/iniparser.c  	$(lib_srcdir)/iniparser.h

libutil_a_SOURCES = $(lib_srcdir)/misc.c 		$(lib_srcdir)/misc.h \
										$(lib_srcdir)/simd.h \
										$(lib_srcdir)/vec3.c		$(lib_srcdir)/vec3.h

libgfunc3_a_SOURCES = $(lib_srcdir)/fft.c  				$(lib_srcdir)/fft.h     $(lib_srcdir)/fft_private.h\
											$(lib_srcdir)/gfunc3.c  		$(lib_srcdir)/gfunc3.h \
											$(lib_srcdir)/mrc.c  				$(lib_srcdir)/mrc.h \
											$(lib_srcdir)/operators.c 	$(lib_srcdir)/operators.h \
											$(lib_srcdir)/vfunc.h											

libetai_a_SOURCES = $(lib_srcdir)/ai_opts.c 			$(lib_srcdir)/ai_opts.h \
										$(lib_srcdir)/ai_params.c  		$(lib_srcdir)/ai_params.h \
										$(lib_srcdir)/et_operators.c 	$(lib_srcdir)/et_operators.h  $(lib_srcdir)/et_operators_private.h\
										$(lib_srcdir)/et_params.c  		$(lib_srcdir)/et_params.h \
										$(lib_srcdir)/et_vfuncs.c  		$(lib_srcdir)/et_vfuncs.h  $(lib_srcdir)/et_vfuncs_private.h \
										$(lib_srcdir)/fwd_op_params.c $(lib_srcdir)/fwd_op_params.h \
										$(lib_srcdir)/fwd_op_opts.c 	$(lib_srcdir)/fwd_op_opts.h \
										$(lib_srcdir)/tiltangles.c  	$(lib_srcdir)/tiltangles.h \
										$(lib_srcdir)/vfunc.h

#lib??_a_LDFLAGS += $(OPENMP_CFLAGS)

