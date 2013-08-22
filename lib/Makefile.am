libdir = $(top_srcdir)/lib
gnulibdir = $(top_srcdir)/gnulib
AM_CPPFLAGS = -I$(top_srcdir)/gnulib -Wall

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
libetreco_a_SOURCES = $(libdir)/fft.c  				$(libdir)/fft.h \
											$(libdir)/gfunc3.c  		$(libdir)/gfunc3.h \
											$(libdir)/mrc.c  				$(libdir)/mrc.h \
											$(libdir)/operators.c 	$(libdir)/operators.h \
											$(libdir)/params.c  		$(libdir)/params.h \
											$(libdir)/tiltangles.c  $(libdir)/tiltangles.h \
											$(libdir)/et_vfuncs.c  	$(libdir)/et_vfuncs.h  	$(libdir)/et_vfuncs_private.h\
											$(libdir)/vfunc.h \
											$(libdir)/misc.c 				$(libdir)/misc.h \
											$(libdir)/options.c 		$(libdir)/options.h \
											$(libdir)/matvec3.c 		$(libdir)/matvec3.h \
											$(libdir)/simd.h \
											$(libdir)/CException.c 	$(libdir)/CException.h 	$(libdir)/CExceptionConfig.h
#libetreco_a_LDFLAGS += $(OPENMP_CFLAGS)

libiniparser_a_SOURCES = $(libdir)/iniparser.c  	$(libdir)/iniparser.h \
												 $(libdir)/dictionary.c  	$(libdir)/dictionary.h
