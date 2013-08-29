/*
 * misc.c -- miscellaneous helper functions
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "CException.h"

#include "misc.h"

#include "mrc.h"


/*-------------------------------------------------------------------------------------------------*/

void *
ali16_malloc (size_t nbytes)
{
  void *alloc_arr;
  
  if (posix_memalign ((void **) &alloc_arr, 16, nbytes))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_NOMEMORY, "Unable to allocate %lu bytes of memory.", nbytes);
      return NULL;
    }
    
  return alloc_arr;
}

/*-------------------------------------------------------------------------------------------------*/

#if !GNULIB_DIRNAME
char const *
base_name (char const *filename)
{
  char const *s;
  
  if ((s = strrchr (filename, '/')) == NULL)
    s = filename;
  else
    s++;
    
  return s;
}
#endif

/*-------------------------------------------------------------------------------------------------*/

enum { COPYRIGHT_YEAR = 2013 };

void 
print_version_etc (char const *progname)
{
  printf ("%s (%s) %s\n", progname, PACKAGE_NAME, PACKAGE_VERSION);
  printf ("Copyright (C) %d %s\n", COPYRIGHT_YEAR, AUTHORS);
  
  puts ("\
\n\
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.");
  
  printf ("\nReport bugs to: %s\n", PACKAGE_BUGREPORT);
  #ifdef PACKAGE_URL
  printf ("%s home page: <%s>\n", PACKAGE_NAME, PACKAGE_URL);
  #endif

  return;
}

/*-------------------------------------------------------------------------------------------------*/

#define SHIFT_PSI 0.0
#define SHIFT_THETA 0.0
#define SHIFT_PHI 0.0

void
compute_rotated_basis (vec3 const angles_deg, vec3 om_x, vec3 om_y, vec3 om_z)
{
  float alpha, beta, gamma;
  
  CAPTURE_NULL (angles_deg);
  CAPTURE_NULL (om_x);
  CAPTURE_NULL (om_y);
  CAPTURE_NULL (om_z);

  alpha = ONE_DEGREE * angles_deg[2] + SHIFT_PHI;
  beta  = ONE_DEGREE * angles_deg[1] + SHIFT_THETA;
  gamma = ONE_DEGREE * angles_deg[0] + SHIFT_PSI;

  // We follow the "x" convention for the Euler angles; https://de.wikipedia.org/wiki/Eulersche_Winkel

  om_x[0] = cosf (alpha) * cosf (gamma) - sinf (alpha) * cosf (beta) * sinf (gamma);
  om_x[1] = sinf (alpha) * cosf (gamma) + cosf (alpha) * cosf (beta) * sinf (gamma);
  om_x[2] = sinf (beta)  * sinf (gamma);

  om_y[0] = -sinf (alpha) * cosf (gamma) - sinf (alpha) * cosf (beta) * cosf (gamma);
  om_y[1] = -sinf (alpha) * sinf (gamma) + cosf (alpha) * cosf (beta) * cosf (gamma);
  om_y[2] =  sinf (beta)  * cosf (gamma);

  om_z[0] = sinf (alpha) * sinf (beta);
  om_z[1] = cosf (alpha) * sinf (beta);
  om_z[2] = cosf (beta);

  return;
}

/*-------------------------------------------------------------------------------------------------*/
