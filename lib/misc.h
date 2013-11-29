/*
 * misc.h -- miscellaneous helper functions
 * 
 * Copyright 2013 Holger Kohr <kohr@num.uni-sb.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef __MISC_H__
#define __MISC_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if GNULIB_DIRNAME
#include "dirname.h"
#endif

#include <stdlib.h>
#include <errno.h>
#include <float.h>

#include "vec3.h"


/*-------------------------------------------------------------------------------------------------*
 * Constants
 *-------------------------------------------------------------------------------------------------*/

#define EPS_DENOM 1E-14
#define EPS_GRID  2E-04
#define ONE_DEGREE 0.017453292519943295F  // [rad]
#define M_SQRT2PI 2.5066282746310002F

#define TRUE  1
#define FALSE 0

static float const c_zero[2] = {0.0f, 0.0f};
static float const c_one[2]  = {1.0f, 0.0f};


/*-------------------------------------------------------------------------------------------------*
 * Debugging 
 *-------------------------------------------------------------------------------------------------*/

#if DEBUGGING

#define PRINT_AT_LINE(_fmt, ...) \
  fprintf (stderr, "[file %s, line %d]: " _fmt, \
  __FILE__, __LINE__,  ##__VA_ARGS__)

#define DEBUG(X) X

#else
#define PRINT_AT_LINE(_fmt, ...)  /* empty */
#define DEBUG(X)  /* empty */

#endif /* DEBUGGING */


/*-------------------------------------------------------------------------------------------------*
 * Function argument checking
 *-------------------------------------------------------------------------------------------------*/

#if ARG_CHECKING

#define CAPTURE_NULL(X, _retval) \
do { \
  if ((X) == NULL) { EXC_THROW_PRINT (EXC_NULL); return _retval; }\
} while (0)

#define CAPTURE_NULL_VOID(X) \
do { \
  if ((X) == NULL) { EXC_THROW_PRINT (EXC_NULL); return; }\
} while (0)

#else
#define CAPTURE_NULL(X, _retval)  /* empty */
#define CAPTURE_NULL_VOID(X)  /* empty */

#endif  /* ARG_CHECKING */


/*-------------------------------------------------------------------------------------------------*
 * Verbosity-dependent output
 *-------------------------------------------------------------------------------------------------*/

extern int verbosity_level;

#define VERB_LEVEL_QUIET    -1
#define VERB_LEVEL_NORMAL    0
#define VERB_LEVEL_VERBOSE   1

/* Print if VERBOSITY_LEVEL is larger than _VLEVEL */
#define PRINT_VERB_LEVEL(_fmt, _vlevel, ...) \
do { \
  if (verbosity_level >= _vlevel) \
    printf (_fmt, ##__VA_ARGS__); \
} while (0)

#define PRINT_QUIET(_fmt, ...)   PRINT_VERB_LEVEL(_fmt, VERB_LEVEL_QUIET,   ##__VA_ARGS__)
#define PRINT_NORMAL(_fmt, ...)  PRINT_VERB_LEVEL(_fmt, VERB_LEVEL_NORMAL,  ##__VA_ARGS__)
#define PRINT_VERBOSE(_fmt, ...) PRINT_VERB_LEVEL(_fmt, VERB_LEVEL_VERBOSE, ##__VA_ARGS__)

/*-------------------------------------------------------------------------------------------------*
 * Helper functions
 *-------------------------------------------------------------------------------------------------*/

/* Allocate a 16-byte aligned array of size NBYTES. Returns NULL on failure. 
 * 
 * Thrown exceptions:
 * - EXC_NOMEMORY
 */
void *
ali16_malloc (size_t nbytes) __attribute__((alloc_size(1)));


/* A simple replacement for the base_name function in Gnulib.  Returns a pointer to the position in
 * FILENAME where the base name (without path) starts.
 */
#if !GNULIB_DIRNAME
char const *
base_name (char const *filename);
#endif


/* Print program name, version, authors, year, license and some other interesting stuff. */
void 
print_version_etc (char const *progname);


/* Compute the rotated basis vectors OM_X, OM_Y, OM_Z in the standard basis.  They correspond to the
 * columns in the rotation matrix defined by the Euler angles in ANGLES_DEG (in degrees).  The 
 * "X" convention as defined in https://de.wikipedia.org/wiki/Eulersche_Winkel is used.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 */
void
compute_rotated_basis (vec3 const angles_deg, vec3 om_x, vec3 om_y, vec3 om_z);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __MISC_H__ */
