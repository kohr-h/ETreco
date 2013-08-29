# Preliminaries
AC_INIT([ETreco], [20130828], [kohr(at)num.uni-sb.de], [], [https://github.com/kohr-h/ETreco])
AC_CONFIG_SRCDIR([src/ai_single_axis.c])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Authors needed for `--version' or `--help' output
AC_DEFINE([AUTHORS], ["Holger Kohr"], [The authors of this software.])

# Checks for programs
AC_PROG_CC

# Gnulib macro (should come right after AC_PROG_CC)
gl_EARLY

# Checks for programs, part 2
AC_PROG_RANLIB
AM_PROG_AR

# Another Gnulib macro
gl_INIT

AC_LANG([C])

# Check for OpenMP
AC_OPENMP

# Activate config file for CException
AC_DEFINE([CEXCEPTION_USE_CONFIG_FILE], [ ], [Define to use a config file for CException.])

# Check for SSE
AC_CACHE_CHECK([whether SSE instructions are available], [ac_cv_cpu_have_sse],
  [ac_cv_cpu_have_sse=no
   AC_RUN_IFELSE(
   [
     AC_LANG_PROGRAM(
       [[\
          void cpuid(int CPUInfo[4],int InfoType){
              __asm__ __volatile__ (
                  "cpuid":
                  "=a" (CPUInfo[0]),
                  "=b" (CPUInfo[1]),
                  "=c" (CPUInfo[2]),
                  "=d" (CPUInfo[3]) :
                  "a" (InfoType)
              );}
        ]],
        [[\
          int SSE = 0;  int info[4];
          cpuid(info, 0);
          int nIds = info[0];
          if (nIds >= 1){
              cpuid(info,0x00000001);
              SSE   = (info[3] & ((int)1 << 25)) != 0;
            }
          return !SSE;  /* 0 for success (SSE) and 1 for failure (no SSE) */
        ]])
   ],
   [ac_cv_cpu_have_sse=yes],
   []
  )]
)

if test "x$ac_cv_cpu_have_sse" = xyes ; then
  AC_DEFINE([HAVE_SSE], [1], [Define to 1 if SSE instructions are available.])
else
  AC_DEFINE([HAVE_SSE], [0], [Define to 1 if SSE instructions are available.])
fi
AM_CONDITIONAL([HAVE_SSE], [test "x$ac_cv_cpu_have_sse" = xyes])

# Check for SSE2
AC_CACHE_CHECK([whether SSE2 instructions are available], [ac_cv_cpu_have_sse2],
  [ac_cv_cpu_have_sse2=no
   AC_RUN_IFELSE(
   [
     AC_LANG_PROGRAM(
       [[\
          void cpuid(int CPUInfo[4],int InfoType){
              __asm__ __volatile__ (
                  "cpuid":
                  "=a" (CPUInfo[0]),
                  "=b" (CPUInfo[1]),
                  "=c" (CPUInfo[2]),
                  "=d" (CPUInfo[3]) :
                  "a" (InfoType)
              );}
        ]],
        [[\
          int SSE2 = 0;  int info[4];
          cpuid(info, 0);
          int nIds = info[0];
          if (nIds >= 1){
              cpuid(info,0x00000001);
              SSE2  = (info[3] & ((int)1 << 26)) != 0;
            }
          return !SSE2;  /* 0 for success (SSE2) and 1 for failure (no SSE2) */
        ]])
   ],
   [ac_cv_cpu_have_sse2=yes],
   []
  )]
)

if test "x$ac_cv_cpu_have_sse2" = xyes ; then
  AC_DEFINE([HAVE_SSE2], [1], [Define to 1 if SSE2 instructions are available.])
else
  AC_DEFINE([HAVE_SSE2], [0], [Define to 1 if SSE2 instructions are available.])
fi
AM_CONDITIONAL([HAVE_SSE2], [test "x$ac_cv_cpu_have_sse2" = xyes])

# Check for SSE3
AC_CACHE_CHECK([whether SSE3 instructions are available], [ac_cv_cpu_have_sse3],
  [ac_cv_cpu_have_sse3=no
   AC_RUN_IFELSE(
   [
     AC_LANG_PROGRAM(
       [[\
          void cpuid(int CPUInfo[4],int InfoType){
              __asm__ __volatile__ (
                  "cpuid":
                  "=a" (CPUInfo[0]),
                  "=b" (CPUInfo[1]),
                  "=c" (CPUInfo[2]),
                  "=d" (CPUInfo[3]) :
                  "a" (InfoType)
              );}
        ]],
        [[\
          int SSE3 = 0;  int info[4];
          cpuid(info, 0);
          int nIds = info[0];
          if (nIds >= 1){
              cpuid(info,0x00000001);
              SSE3  = (info[2] & ((int)1 <<  0)) != 0;
            }
          return !SSE3;  /* 0 for success (SSE3) and 1 for failure (no SSE3) */
        ]])
   ],
   [ac_cv_cpu_have_sse3=yes],
   []
  )]
)

if test "x$ac_cv_cpu_have_sse3" = xyes ; then
  AC_DEFINE([HAVE_SSE3], [1], [Define to 1 if SSE3 instructions are available.])
else
  AC_DEFINE([HAVE_SSE3], [0], [Define to 1 if SSE3 instructions are available.])
fi
AM_CONDITIONAL([HAVE_SSE3], [test "x$ac_cv_cpu_have_sse3" = xyes])

# Check for SSSE3
AC_CACHE_CHECK([whether SSSE3 instructions are available], [ac_cv_cpu_have_ssse3],
  [ac_cv_cpu_have_ssse3=no
   AC_RUN_IFELSE(
   [
     AC_LANG_PROGRAM(
       [[\
          void cpuid(int CPUInfo[4],int InfoType){
              __asm__ __volatile__ (
                  "cpuid":
                  "=a" (CPUInfo[0]),
                  "=b" (CPUInfo[1]),
                  "=c" (CPUInfo[2]),
                  "=d" (CPUInfo[3]) :
                  "a" (InfoType)
              );}
        ]],
        [[\
          int SSSE3 = 0;  int info[4];
          cpuid(info, 0);
          int nIds = info[0];
          if (nIds >= 1){
              cpuid(info,0x00000001);
              SSSE3 = (info[2] & ((int)1 <<  9)) != 0;
            }
          return !SSSE3;  /* 0 for success (SSSE3) and 1 for failure (no SSSE3) */
        ]])
   ],
   [ac_cv_cpu_have_ssse3=yes],
   []
  )]
)

if test "x$ac_cv_cpu_have_ssse3" = xyes ; then
  AC_DEFINE([HAVE_SSSE3], [1], [Define to 1 if SSSE3 instructions are available.])
else
  AC_DEFINE([HAVE_SSSE3], [0], [Define to 1 if SSSE3 instructions are available.])
fi
AM_CONDITIONAL([HAVE_SSSE3], [test "x$ac_cv_cpu_have_ssse3" = xyes])

# Check for SSE4.1
AC_CACHE_CHECK([whether SSE4.1 instructions are available], [ac_cv_cpu_have_sse41],
  [ac_cv_cpu_have_sse41=no
   AC_RUN_IFELSE(
   [
     AC_LANG_PROGRAM(
       [[\
          void cpuid(int CPUInfo[4],int InfoType){
              __asm__ __volatile__ (
                  "cpuid":
                  "=a" (CPUInfo[0]),
                  "=b" (CPUInfo[1]),
                  "=c" (CPUInfo[2]),
                  "=d" (CPUInfo[3]) :
                  "a" (InfoType)
              );}
        ]],
        [[\
          int SSE41 = 0;  int info[4];
          cpuid(info, 0);
          int nIds = info[0];
          if (nIds >= 1){
              cpuid(info,0x00000001);
              SSE41 = (info[2] & ((int)1 << 19)) != 0;
            }
          return !SSE41;  /* 0 for success (SSE4.1) and 1 for failure (no SSE4.1) */
        ]])
   ],
   [ac_cv_cpu_have_sse41=yes],
   []
  )]
)

if test "x$ac_cv_cpu_have_sse41" = xyes ; then
  AC_DEFINE([HAVE_SSE41], [1], [Define to 1 if SSE4.1 instructions are available.])
else
  AC_DEFINE([HAVE_SSE41], [0], [Define to 1 if SSE4.1 instructions are available.])
fi
AM_CONDITIONAL([HAVE_SSE41], [test "x$ac_cv_cpu_have_sse41" = xyes])


# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T


# Inlining is used by GSL, too
if test "x$ac_cv_c_inline" != xno ; then
  AC_DEFINE([HAVE_INLINE], [1], [Define to 1 if inline is available (for GSL).])
  AC_SUBST([HAVE_INLINE])
fi

# Range checking
# Disabled by `./configure --disable-range-checking'
AC_ARG_ENABLE(
  [range-checking],
  [AS_HELP_STRING([--disable-range-checking], 
                  [disable range checking of array indices (USE WITH CAUTION!)])],
  [AC_DEFINE([RANGE_CHECKING], [0], [Define to 1 to enable range checking of array indices.])
   AC_SUBST([RANGE_CHECKING])],
  [AC_DEFINE([RANGE_CHECKING], [1], [Define to 1 to enable range checking of array indices.])
   AC_SUBST([RANGE_CHECKING])])

# Function argument checking
# Disabled by `./configure --disable-arg-checking'
AC_ARG_ENABLE(
  [arg-checking],
  [AS_HELP_STRING([--disable-arg-checking],
                  [disable checking of function arguments (USE WITH CAUTION!)])],
  [AC_DEFINE([ARG_CHECKING], [0], [Define to 1 to enable checking of function arguments.])
   AC_SUBST([ARG_CHECKING])],
  [AC_DEFINE([ARG_CHECKING], [1], [Define to 1 to enable checking of function arguments.])
   AC_SUBST([ARG_CHECKING])])

# Debugging output and checks (likely to slow down execution)
# Enabled by `./configure --enable-debug' 
AC_ARG_ENABLE(
  [debugging],
  [AS_HELP_STRING([--enable-debugging],
                  [enable debugging checks and output])],
  [AC_DEFINE([DEBUGGING], [1], [Define to 1 for debugging output and argument checks.])
   AC_SUBST([DEBUGGING])],
  [AC_DEFINE([DEBUGGING], [0], [Define to 1 for debugging output and argument checks.])
   AC_SUBST([DEBUGGING])])


# Checks for library functions
AC_FUNC_MALLOC
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNC([powf], [], AC_CHECK_LIB([m], [powf]))
AC_CHECK_FUNC([sqrtf], [], AC_CHECK_LIB([m], [sqrtf]))
AC_CHECK_FUNCS([strchr strrchr strtol strdup memset mkdir])


# Checks for headers
AC_CHECK_HEADERS([stdlib.h string.h error.h float.h stddef.h limits.h])

# BLAS related checks
AC_SEARCH_LIBS([cblas_scopy], \
               [atlas cblas goto-blas gslcblas], \
               [AC_DEFINE([HAVE_CBLAS], [1], [Define to 1 if CBLAS is available.]) \
                AC_SUBST([HAVE_CBLAS])], \
               [AC_MSG_WARN([No CBLAS implementation found - using naive methods instead.])])

# GSL related checks
AC_CHECK_LIB([gsl], [gsl_sf_bessel_J0], [],\
                [AC_MSG_ERROR([GSL not found.])])
AC_CHECK_HEADER([gsl/gsl_math.h], [], \
                [AC_MSG_ERROR([GSL headers not found.])])
AC_CHECK_HEADER([gsl/gsl_sf_bessel.h], [], \
                [AC_MSG_ERROR([GSL headers not found.])])
AC_CHECK_HEADER([gsl/gsl_blas.h], [], \
                [AC_MSG_ERROR([GSL headers not found.])])

# FFWT related checks
AC_CHECK_LIB([fftw3f], [fftwf_malloc], [], \
                [AC_MSG_ERROR([FFTW3 single-precision library not found.])])
AC_CHECK_HEADERS([fftw3.h], [], \
                [AC_MSG_ERROR([FFTW3 header not found.])])

# Targets to configure
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 lib/Makefile \
                 gnulib/Makefile \
                 testsuite/Makefile \
                 ])
                 
# Output
AC_OUTPUT
