/*
 * convex_comb.c -- compute a convex combination of 2 MRC files
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


int verbosity_level = 0;


int 
main (int argc, char **argv)
{
  CEXCEPTION_T e = EXC_NONE;
  char *fname1, *fname2, *fname_out;
  char const *progname = base_name (argv[0]);
  float lambda;
  float min1, min2, max1, max2;
  
  gfunc3 *gf1, *gf2;
  
  Try
  {
    if (argc != 5)
      {
        print_version_etc (progname);
        fprintf (stderr, "\nUsage: %s <file1> <file2> <lambda> <outfile>\n\n", progname);
        return 0;
      }
    
    fname1 = argv[1];
    fname2 = argv[2];
    lambda = atof (argv[3]);
    fname_out = argv[4];
    
    gf1 = new_gfunc3 ();
    gf2 = new_gfunc3 ();
    
    gfunc3_init_mrc (gf1, fname1, NULL, NULL, VOLUME);
    gfunc3_init_mrc (gf2, fname2, NULL, NULL, VOLUME);

    /* Shift and scale both to [0,1] and add (1-lambda)*f1 + lambda*f2 */
    min1 = gfunc3_min (gf1);
    max1 = gfunc3_max (gf1);
    PRINT_AT_LINE ("Function 1 min: %e\n", min1);
    PRINT_AT_LINE ("Function 1 max: %e\n", max1);
    gfunc3_add_constant (gf1, -min1);
    gfunc3_scale (gf1, 1.0 / (max1 - min1));
    
    min2 = gfunc3_min (gf2);
    max2 = gfunc3_max (gf2);
    PRINT_AT_LINE ("Function 2 min: %e\n", min2);
    PRINT_AT_LINE ("Function 2 max: %e\n", max2);
    gfunc3_add_constant (gf2, -min2);
    gfunc3_scale (gf2, 1.0 / (max2 - min2));
    
    gfunc3_scale (gf2, lambda);
    gfunc3_axpy (1.0 - lambda, gf1, gf2);
    
    gfunc3_to_mrc (gf1, fname_out);

    gfunc3_free (&gf1);
    gfunc3_free (&gf2);
  }
  Catch (e)
  {
    EXC_REPRINT;
  }
    
  return 0;
}

