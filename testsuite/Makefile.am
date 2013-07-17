AM_CPPFLAGS = -I$(srcdir)/../lib -I$(srcdir)/../gnulib
AUTOMAKE_OPTIONS = dejagnu

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

bin_PROGRAMS =  gf3_test_grid cexceptions

gf3_test_grid_SOURCES = gf3_test_grid.c
gf3_test_grid_LDFLAGS = -lm
gf3_test_grid_LDADD = ../lib/libetreco.a \
											../gnulib/libgnu.a

cexceptions_SOURCES = cexceptions.c
cexceptions_LDADD = ../lib/libetreco.a \
									  ../lib/libiniparser.a \
  								  ../gnulib/libgnu.a
