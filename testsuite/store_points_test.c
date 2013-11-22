/*
 * store_points_test.c
 * 
 * Copyright 2013 Holger Kohr <kohr@num.uni-sb.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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

#include "ETreco.h"


int main(int argc, char **argv)
{
  int i;
  vec3 x0 = {0.2, 0.0, 0.0};
  vec3 cs = {0.1, 45.0 * ONE_DEGREE, 60.0 * ONE_DEGREE};
  idx3 shp = {5, 4, 6};

  float *pts = NULL;
  
  gfunc3 *gf = new_gfunc3 ();
  
  gfunc3_init (gf, x0, cs, shp, VOLUME);
  
  // gfunc3_print_grid (gf, "Polar grid:");
  
  pts = gfunc3_grid_points (gf, EUCLIDEAN);
  
  for (i = 0; i < gf->ntotal; i++)
    printf ("%g %g %g\n", pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
  
  free (pts);
  gfunc3_free (&gf);
  
  return 0;
}

