libdir = $(top_builddir)/lib
gnulibdir = $(top_builddir)/gnulib
AM_CPPFLAGS = -I$(top_srcdir)/lib -I$(top_srcdir)/gnulib -Wall

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

bin_PROGRAMS =  gf3_test_grid cexceptions ft_shift_test

gf3_test_grid_SOURCES = $(srcdir)/gf3_test_grid.c
gf3_test_grid_LDFLAGS = -lm
gf3_test_grid_LDADD = $(libdir)/libetreco.a \
											$(gnulibdir)/libgnu.a

cexceptions_SOURCES = $(srcdir)/cexceptions.c
cexceptions_LDADD = $(libdir)/libetreco.a \
									  $(libdir)/libiniparser.a \
  								  $(gnulibdir)/libgnu.a

ft_shift_test_SOURCES = $(srcdir)/ft_shift_test.c
ft_shift_test_LDADD = $(libdir)/libetreco.a \
											$(libdir)/libiniparser.a \
											$(gnulibdir)/libgnu.a
