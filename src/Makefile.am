libdir = $(top_builddir)/lib
gnulibdir = $(top_builddir)/gnulib
AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_srcdir)/gnulib -Wall

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
ai_single_axis_SOURCES = $(srcdir)/ai_single_axis.c $(top_srcdir)/lib/ETreco.h
ai_single_axis_LDFLAGS = -lm
ai_single_axis_LDADD = $(libdir)/libetreco.a \
											 $(libdir)/libiniparser.a \
											 $(gnulibdir)/libgnu.a \
											 $(OPENMP_CFLAGS)

generate_tiltangles_SOURCES = $(srcdir)/generate_tiltangles.c $(top_srcdir)/lib/ETreco.h
generate_tiltangles_LDADD = $(libdir)/libetreco.a \
											      $(gnulibdir)/libgnu.a

convex_comb_SOURCES = $(srcdir)/convex_comb.c $(top_srcdir)/lib/ETreco.h
convex_comb_LDADD = $(libdir)/libetreco.a \
										$(gnulibdir)/libgnu.a
