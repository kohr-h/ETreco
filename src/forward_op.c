/*
 * forward_op_test.c
 * 
 * Copyright 2014 Holger Kohr <kohr@num.uni-sb.de>
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

  int last_index;

  FwdOpts *opts = NULL;
  EtParams *et_params = NULL;
  FwdParams *fwd_params = NULL;
  tiltangles *tilts = NULL;
  gfunc3 *volume = NULL, *tiltseries = NULL;
  
  Try { 
    opts = new_FwdOpts (); 
    et_params = new_EtParams (); 
    fwd_params = new_FwdParams (); 
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    tiltseries = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);


  /* Assign options, parameters, tiltangles and init gfunc structures */
  Try { FwdOpts_assign_from_args (opts, argc, argv); } CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    FwdParams_assign_from_file (fwd_params, opts->fname_params); 
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
    gfunc3_init_mrc (volume, opts->fname_in, NULL, NULL, VOLUME);
    // gfunc3_init (tiltseries, NULL, opts)
  }  CATCH_EXIT_FAIL (_e);
  
  
  /* Scale and shift the volume if necessary */
  Try {
    FwdParams_apply_to_volume (fwd_params, volume);
    gfunc3_scale_grid (volume, et_params->magnification); 
  }  CATCH_EXIT_FAIL (_e);


  /* Set some parameters and resolve conflicts */
  if (opts->num_images == 0)  /* Not specified by option, thus taken from data */
    opts->num_images = tilts->ntilts - opts->start_index;

  last_index = opts->start_index + opts->num_images - 1;

  Try {
    if (last_index >= tilts->ntilts)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Image indices (max: %d) must be "
      "smaller\n than the number of tiltangles (%d).\n", last_index, tilts->ntilts);
  } CATCH_EXIT_FAIL (_e);


  /* Print a summary of everything before starting */
  FwdOpts_print (opts);
  EtParams_print (et_params);
  FwdParams_print (fwd_params);
  gfunc3_print_grid (volume, "volume grid:");
  
  FwdOpts_free (&opts);
  EtParams_free (&et_params);
  FwdParams_free (&fwd_params);
  tiltangles_free (&tilts);
  gfunc3_free (&volume);
  
  return EXIT_SUCCESS;
  
  /*
  gfunc3 *phantom = new_gfunc3 (), *proj_img = new_gfunc3 (), *proj_r = new_gfunc3();
  OptionData *opt_data = new_OptionData ();
  RecParams *rec_p = new_RecParams ();
  
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
  */
}

