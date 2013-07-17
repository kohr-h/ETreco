AM_CPPFLAGS = -I$(srcdir)/../lib -I$(srcdir)/../gnulib

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

bin_PROGRAMS =  ai_single_axis generate_tiltangles convex_comb
ai_single_axis_SOURCES = ai_single_axis.c $(srcdir)/../lib/ETreco.h
ai_single_axis_LDFLAGS = -lm
ai_single_axis_LDADD = ../lib/libetreco.a \
											 ../lib/libiniparser.a \
											 ../gnulib/libgnu.a \
											 $(OPENMP_CFLAGS)

generate_tiltangles_SOURCES = generate_tiltangles.c
generate_tiltangles_LDADD = ../lib/libetreco.a \
											      ../gnulib/libgnu.a

convex_comb_SOURCES = convex_comb.c
convex_comb_LDADD = ../lib/libetreco.a \
										../gnulib/libgnu.a
