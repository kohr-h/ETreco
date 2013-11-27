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

#include "ETreco.h"

int verbosity_level = VERB_LEVEL_NORMAL;
int autocenter_vol_flag = TRUE;
int fft_padding = 0;

int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i;
  vec3 x0 = {0.1, 0, 0}, cs = {0.2, 0.3, 0.5}, cs_plane = {0.1, 0.1, 1.0}, angles = {0.0, 45.0, 45.0};
  idx3 shp = {10, 7, 3}, padding = {2, 2, 2}, shp_plane = {25, 24, 1};

  float *freqs = NULL, wave_number = 45.0;
  
  gfunc3 *vol = new_gfunc3 (), *perp_plane_re = new_gfunc3 (), *perp_plane_im = new_gfunc3 ();
  
  gfunc3_init (vol, x0, cs, shp, VOLUME);
  gfunc3_init (perp_plane_re, NULL, cs_plane, shp_plane, VOLUME);
  gfunc3_init (perp_plane_im, NULL, cs_plane, shp_plane, VOLUME);

  gfunc3_set_all (vol, c_one);
  gfunc3_zeropad (vol, padding);
  gfunc3_to_mrc (vol, "temp/vol.mrc");

  gfunc3_grid_fwd_reciprocal (perp_plane_re);
  gfunc3_grid_fwd_reciprocal (perp_plane_im);

  Try { freqs = perp_plane_freqs (perp_plane_re, angles); } CATCH_EXIT_FAIL (_e);
  nfft3_transform (vol, NULL, freqs, perp_plane_re->ntotal, 
    perp_plane_re->fvals, perp_plane_im->fvals);
  
  gfunc3_to_mrc (perp_plane_re, "temp/plane_re.mrc");
  gfunc3_to_mrc (perp_plane_re, "temp/plane_im.mrc");

  
  free (freqs);
  gfunc3_free (&vol);
  gfunc3_free (&perp_plane_re);
  gfunc3_free (&perp_plane_im);
  
  return 0;
}

