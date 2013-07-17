/*
 * gf3_test_grid.c -- test some of the gf3 functions related to grid
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


#include <stdio.h>
#include <stdlib.h>

#include "ETreco.h"

int test_subgrid_indices (void)
{
  size_t i, *idcs;
  
  idx3 shp = {11, 8, 9}, shp_sub = {8, 2, 3};
  vec3 cs = {1.0, 0.3, 2.7}, cs_sub = {1.0, 0.9, 5.4};
  vec3 x0 = {0.0, -6.3, 2.4}, x0_sub = {0.0 , -5.4, 7.8};
  
  gfunc3 *gf, *gf_sub;
  
  gf = new_gfunc3 ();
  gfunc3_init (gf, x0, cs, shp, REAL);
  
  gf_sub = new_gfunc3 ();
  gfunc3_init (gf_sub, x0_sub, cs_sub, shp_sub, REAL);
  
  gfunc3_print_grid (gf, "gf grid");
  gfunc3_print_grid (gf_sub, "gf_sub grid");
  
  idcs = gfunc3_subgrid_flatidcs (gf, gf_sub);
  
  if (idcs == NULL)
    return -1;
  
  for (i = 0; i < gf_sub->ntotal; i++)
    printf ("%lu\n", idcs[i]);
    
  free (idcs);
  
  return 0;
}

int test_zeropad (void)
{
  idx3 shp = {11, 8, 9};
  vec3 cs = {1.0, 0.3, 2.7};
  vec3 x0 = {0.0, -6.3, 2.4};
  
  idx3 padding = {1, 2, 0};
  
  float one = 1.0;
  
  gfunc3 *gf;
  
  gf = new_gfunc3 ();
  gfunc3_init (gf, x0, cs, shp, REAL);
  gfunc3_set_all (gf, &one);

  gfunc3_print_grid (gf, "gf grid");
  
  gfunc3_zeropad (gf, padding);
  gfunc3_print_grid (gf, "after zero-padding");

  gfunc3_to_mrc (gf, "padded.mrc");

  gfunc3_unpad (gf, padding);
  gfunc3_print_grid (gf, "after unpadding");
  
  gfunc3_free (&gf);
  
  return 0;
}

int main(int argc, char **argv)
{
  // if (test_subgrid_indices() != 0)
    // puts ("Test failed.");
    
  if (test_zeropad () != 0)
    puts ("Test failed.");
  return 0;
}

