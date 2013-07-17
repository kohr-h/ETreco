AM_CPPFLAGS = -I$(srcdir)/../gnulib -Wall

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

noinst_LIBRARIES = libetreco.a libiniparser.a
libetreco_a_SOURCES = fft.c  fft.h \
											gfunc3.c  gfunc3.h \
											mrc.c  mrc.h \
											operators.c operators.h \
											params.c  params.h \
											tiltangles.c  tiltangles.h \
											et_vfuncs.c  et_vfuncs.h  et_vfuncs_private.h\
											vfunc.h \
											misc.c misc.h \
											options.c options.h \
											matvec3.c matvec3.h \
											simd.h \
											CException.c CException.h CExceptionConfig.h
#libetreco_a_LDFLAGS += $(OPENMP_CFLAGS)

libiniparser_a_SOURCES = iniparser.c  iniparser.h \
												 dictionary.c  dictionary.h
