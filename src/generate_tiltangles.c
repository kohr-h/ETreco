/*
 * generate_tiltangles.c -- generate a tiltangles.txt file from starting
 *                          angle and increment
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

#if GNULIB_DIRNAME
#include "dirname.h"
#endif

int main (int argc, char **argv)
{
  CEXCEPTION_T e = EXC_NONE;
  int i;
  
  int ntilts;
  float theta_st, theta_inc;
  char *fname_out;
  tiltangles *ta = new_tiltangles ();
  
  if (argc != 5)
    {
      print_version_etc (base_name (argv[0]));
      printf ("\nUsage: %s <start_angle> <increment> <num_tilts> <outfile>\n\n", base_name(argv[0]));
      return 0;
    }
  
  theta_st  = atof (argv[1]);
  theta_inc = atof (argv[2]);
  ntilts    = atoi (argv[3]);
  fname_out = argv[4];
  
  Try
  {
    tiltangles_init (ta, ntilts, 3);
    
    for (i = 0; i < ntilts; i++)
      ta->angles_deg[i][1] = theta_st + i * theta_inc;
    
    printf ("Writing tiltangles to %s\n", fname_out);
    tiltangles_to_file (ta, fname_out);
    
    tiltangles_free (&ta);
  }
  Catch (e)
  {
    EXC_REPRINT;
  }
    
  return 0;
}

