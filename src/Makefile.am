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


bin_PROGRAMS =  ai_single_axis generate_tiltangles convex_comb

# CAUTION: order of libraries matters! Order from "most dependent" to "least dependent"
# or list each library twice

ai_single_axis_SOURCES = $(srcdir)/ai_single_axis.c $(top_srcdir)/lib/ETreco.h
ai_single_axis_LDADD = 	-letai -lgfunc3 -lgnu -lutil -liniparser -lCException -lm
												

generate_tiltangles_SOURCES = $(srcdir)/generate_tiltangles.c
generate_tiltangles_LDADD = -letai -lutil -lCException -lgnu

convex_comb_SOURCES = $(srcdir)/convex_comb.c
convex_comb_LDADD = -lgfunc3 -lutil -lCException -lgnu
