/*
 * forward_op_test.c
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

#include "ETreco.h"

int verbosity_level = VERB_LEVEL_NORMAL;
int autocenter_vol_flag = TRUE;
int fft_padding = 0;

int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  char const *phantom_name = "../test_data/Testdata/Balls/phantom_balls.mrc";
  vec3 cs_proj = {8.0, 8.0, 1.0}, x0_proj = {0.0, 5.0, 0.0};
  idx3 shp_proj = {210, 250, 1};
  vec3 angles = {90.0, 90.0, 90.0};
  gfunc3 *phantom = new_gfunc3 (), *xr_proj = new_gfunc3 (), *proj_r = new_gfunc3();
  

  
  
  Try { 
    gfunc3_init_mrc (phantom, phantom_name, NULL, NULL, VOLUME); 
    gfunc3_init (xr_proj, x0_proj, cs_proj, shp_proj, COMPLEX); 
  } CATCH_EXIT_FAIL (_e);
    
  gfunc3_print_grid (phantom, "Phantom grid:");
  gfunc3_print_grid (xr_proj, "Plane grid:");
  
  Try { 
    gfunc3_real2complex (phantom);
    xray_projection (phantom, angles, xr_proj); 
  } CATCH_EXIT_FAIL (_e);
  
  Try { gfunc3_to_mrc (xr_proj, "xray.mrc"); }  CATCH_EXIT_FAIL (_e);
  
  Try { 
    proj_r = gfunc3_realpart (xr_proj, NULL); 
    gfunc3_to_mrc (proj_r, "xray_re.mrc");
  } CATCH_EXIT_FAIL (_e);
  
  Try { 
    proj_r = gfunc3_imagpart (xr_proj, proj_r); 
    gfunc3_to_mrc (proj_r, "xray_im.mrc");
  } CATCH_EXIT_FAIL (_e);
  
  gfunc3_free (&phantom);
  gfunc3_free (&xr_proj);
  gfunc3_free (&proj_r);
  
  return EXIT_SUCCESS;
}

