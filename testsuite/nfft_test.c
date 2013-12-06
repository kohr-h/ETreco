/*
 * nfft_test.c
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
#include <complex.h>

#include "ETreco.h"

int verbosity_level = VERB_LEVEL_NORMAL;
int autocenter_vol_flag = TRUE;
int fft_padding = 0;

int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;
  vec3 x0 = {0.1, 0, 0}, cs = {1.0, 1.0, 1.0}, cs_plane = {0.2, 0.2, 1.0}, angles = {0.0, 45.0, 45.0};
  idx3 shp = {10, 7, 3}, padding = {2, 2, 2}, shp_plane = {10, 9, 1};

  float *freqs = NULL, wave_number = 45.0;
  
  gfunc3 *vol = new_gfunc3 (), *perp_plane = new_gfunc3 (), *perp_plane_re = NULL;
  
  gfunc3_init (vol, x0, cs, shp, COMPLEX);
  gfunc3_init (perp_plane, NULL, cs_plane, shp_plane, COMPLEX);

  Try {
    gfunc3_set_all (vol, 1.0 - 1.5 * I);
    gfunc3_zeropad (vol, padding); 
    gfunc3_to_mrc (vol, "temp/vol.mrc");
  } CATCH_EXIT_FAIL (_e);

  gfunc3_print_grid (vol, "Volume grid:");
  gfunc3_print_grid (perp_plane, "Perp plane reciprocal grid:");

  Try { 
    fft_forward (vol);
    gfunc3_to_mrc (vol, "temp/vol_ft.mrc");
    fft_backward (vol);
    gfunc3_to_mrc (vol, "temp/vol_ift.mrc");
  } CATCH_EXIT_FAIL (_e);

  Try { freqs = perp_plane_freqs (perp_plane, angles); } CATCH_EXIT_FAIL (_e);

  Try { 
    nfft3_transform (vol, freqs, perp_plane->ntotal, (float complex *) perp_plane->fvals);
  } CATCH_EXIT_FAIL (_e);
  
  Try { gfunc3_to_mrc (perp_plane, "temp/plane_ft.mrc"); } CATCH_EXIT_FAIL (_e);
  
  Try { fft_backward (perp_plane); } CATCH_EXIT_FAIL (_e);
  gfunc3_print_grid (perp_plane, "Perp plane grid:");
  Try { gfunc3_to_mrc (perp_plane, "temp/plane.mrc"); } CATCH_EXIT_FAIL (_e);
  
  Try { perp_plane_re = gfunc3_realpart (perp_plane, NULL); } CATCH_EXIT_FAIL (_e);
  gfunc3_print_grid (perp_plane_re, "Real part plane grid:");
  Try { gfunc3_to_mrc (perp_plane_re, "temp/plane_re.mrc"); } CATCH_EXIT_FAIL (_e);
  Try { gfunc3_real2complex (perp_plane_re); } CATCH_EXIT_FAIL (_e);
  Try { gfunc3_to_mrc (perp_plane_re, "temp/plane_re_ascomplex.mrc"); } CATCH_EXIT_FAIL (_e);
  
  free (freqs);
  gfunc3_free (&vol);
  gfunc3_free (&perp_plane);
  gfunc3_free (&perp_plane_re);
  
  return 0;
}

