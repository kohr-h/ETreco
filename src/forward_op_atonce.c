/*
 * forward_op_atonce.c -- Implementation of the ET forward operator
 * creating the image stack in one go
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
#include <complex.h>

#include "ETreco.h"


int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;

  FwdOpts *opts = NULL;
  EtParams *et_params = NULL;
  FwdParams *fwd_params = NULL;
  tiltangles *tilts = NULL;
  gfunc3 *volume = NULL, *proj_stack = NULL, *proj_stack_imag = NULL;
  
  Try { 
    opts = new_FwdOpts (); 
    et_params = new_EtParams (); 
    fwd_params = new_FwdParams (); 
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    proj_stack = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);


  /* Assign options, parameters, tiltangles and init gfunc structures */
  Try { FwdOpts_assign_from_args (opts, argc, argv); } CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    FwdParams_assign_from_file (fwd_params, opts->fname_params); 
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
    gfunc3_init_mrc (volume, opts->fname_in, NULL, NULL);
  }  CATCH_EXIT_FAIL (_e);

  /* Scale and shift the volume if necessary; transfer to complex */
  Try {
    FwdParams_apply_to_volume (fwd_params, volume);
    gfunc3_scale_grid (volume, et_params->magnification);
    gfunc3_real2complex (volume);
  }  CATCH_EXIT_FAIL (_e);


  /* Initialize the stack */
  Try {
    idx3_copy (proj_stack->shape, fwd_params->detector_shape);
    proj_stack->shape[2] = tilts->ntilts;
    gfunc3_init (proj_stack, NULL, fwd_params->detector_px_size, proj_stack->shape, COMPLEX);
  } CATCH_EXIT_FAIL (_e);

  
  /* Print a summary of everything before starting */
  FwdOpts_print (opts);
  EtParams_print (et_params);
  FwdParams_print (fwd_params);
  gfunc3_print_grid (volume, "volume grid:");
  gfunc3_print_grid (proj_stack, "Projection stack grid:");

  
  /* Compute the projections */
  Try {
    et_scattering_projection_atonce (volume, tilts, et_params, proj_stack, opts->model);
  } CATCH_EXIT_FAIL (_e);
      
  /* Get the imaginary part and write it to disk */
  Try { proj_stack_imag = gfunc3_imagpart (proj_stack, NULL); }  CATCH_EXIT_FAIL (_e);
  Try { gfunc3_to_mrc (proj_stack_imag, opts->fname_out, NULL); }  CATCH_EXIT_FAIL (_e);

  printf ("\n\nProjections written to %s.\n\n", opts->fname_out);
  
  /* Delete all objects before exiting */
  FwdOpts_free (&opts);
  EtParams_free (&et_params);
  FwdParams_free (&fwd_params);
  tiltangles_free (&tilts);
  gfunc3_free (&volume);
  gfunc3_free (&proj_stack);
  gfunc3_free (&proj_stack_imag);
  
  return EXIT_SUCCESS;
}
