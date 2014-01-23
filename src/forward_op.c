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


int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  gfunc3 *phantom = new_gfunc3 (), *proj_img = new_gfunc3 (), *proj_r = new_gfunc3();
  OptionData *opt_data = new_OptionData ();
  RecParams *rec_p = new_RecParams ();
  
  /* WIP
   * TODO: write fwd_op_options!!
   */
  Try { OptionData_assign_from_args (opt_data, argc, argv); }  CATCH_EXIT_FAIL (_e);
  Try { 
    RecParams_assign_from_OptionData (rec_p, opt_data);
    gfunc3_init_mrc (phantom, phantom_name, NULL, NULL, VOLUME); 
    gfunc3_init (xr_proj, x0_proj, cs_proj, shp_proj, COMPLEX); 
  } CATCH_EXIT_FAIL (_e);
    
  gfunc3_print_grid (phantom, "Phantom grid:");
  
  gfunc3_print_grid (xr_proj, "Plane grid:");
  
  Try { 
    gfunc3_real2complex (phantom);
    // xray_projection (phantom, angles, xr_proj); 
    et_scattering_projection (phantom, angles, rec_p, xr_proj, PROJ_ASSUMPTION); 
  } CATCH_EXIT_FAIL (_e);

  Try { gfunc3_to_mrc (xr_proj, "xray.mrc"); }  CATCH_EXIT_FAIL (_e);
  
  Try { 
    proj_r = gfunc3_realpart (xr_proj, NULL); 
    gfunc3_to_mrc (proj_r, "born_re.mrc");
  } CATCH_EXIT_FAIL (_e);
  
  Try { 
    proj_r = gfunc3_imagpart (xr_proj, proj_r); 
    gfunc3_to_mrc (proj_r, "born_im.mrc");
  } CATCH_EXIT_FAIL (_e);
  
  gfunc3_free (&phantom);
  gfunc3_free (&xr_proj);
  gfunc3_free (&proj_r);
  RecParams_free (&rec_p);
  
  return EXIT_SUCCESS;
}

