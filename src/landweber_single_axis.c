/*
 * landweber_single_axis.c -- Implementation of a Landweber iteration for 
 * single-axis ET
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

  int i;

  LandwOpts *opts = NULL;
  EtParams *et_params = NULL;
  FwdParams *fwd_params = NULL;
  tiltangles *tilts = NULL;
  gfunc3 *volume = NULL, *proj_stack = NULL, *proj_stack_imag = NULL, *cur_img = NULL;
  
  Try { 
    opts = new_LandwOpts (); 
    et_params = new_EtParams (); 
    fwd_params = new_FwdParams (); 
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    proj_stack = new_gfunc3 ();
    cur_img = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);


  /* Assign options, parameters, tiltangles and init gfunc structures */
  Try { LandwOpts_assign_from_args (opts, argc, argv); } CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    FwdParams_assign_from_file (fwd_params, opts->fname_params); 
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
    gfunc3_init_mrc (proj_stack, opts->fname_in, NULL, NULL);
  }  CATCH_EXIT_FAIL (_e);

  Try {
    if (tilts->ntilts != proj_stack->shape[2])
      {
        EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Number of tilts (%d) must agree with stack "
        "size (%d).", tilts->ntilts, proj_stack->shape[2]);
      }
  } CATCH_EXIT_FAIL (_e);
    
  if (normalize_flag)
    {
      for (i = 0; i < tilts->ntilts; i++)
        {
          Try { gfunc3_set_stack_pointer (cur_img, proj_stack, i); }  CATCH_EXIT_FAIL (_e);
          Try { 
            histogram_normalization (cur_img, opts->bg_patch_ix0, opts->bg_patch_shape);
          } CATCH_EXIT_FAIL (_e);
        }
    }
  

  /* TODO: continue here */
  
  /* Further steps:
   * - Init volume
   * - Init residual
   * - Get initial guess
   * - back-project
   * - Start iteration
   */
  
// 
// 
  // /* Scale and shift the volume if necessary; transfer to complex and add imaginary part */
  // Try {
    // FwdParams_apply_to_volume (fwd_params, volume);
    // gfunc3_scale_grid (volume, et_params->magnification);
    // if (volume->type == REAL)
      // gfunc3_scale (volume, 1.0 + 0.1 * I);  /* TODO: use acr from config file */
  // }  CATCH_EXIT_FAIL (_e);
// 
// 
  // /* Initialize the stack */
  // Try {
    // idx3_copy (proj_stack->shape, fwd_params->detector_shape);
    // proj_stack->shape[2] = tilts->ntilts;
    // gfunc3_init (proj_stack, NULL, fwd_params->detector_px_size, proj_stack->shape, COMPLEX);
  // } CATCH_EXIT_FAIL (_e);
// 
// 
  // /* Print a summary of everything before starting */
  // LandwOpts_print (opts);
  // EtParams_print (et_params);
  // FwdParams_print (fwd_params);
  // gfunc3_print_grid (volume, "volume grid:");
  // gfunc3_print_grid (proj_stack, "Projection stack grid:");
// 
  // 
  // /* Compute the projections */
  // Try {
    // et_scattering_projection_atonce (volume, tilts, et_params, proj_stack, opts->model);
  // } CATCH_EXIT_FAIL (_e);
// 
  // /* Get the imaginary part, invert contrast if necessary and write it to disk */
  // Try { proj_stack_imag = gfunc3_imagpart (proj_stack, NULL); }  CATCH_EXIT_FAIL (_e);
  // if (invert_contrast_flag)
    // gfunc3_scale (proj_stack_imag, -1.0);
  // Try { gfunc3_to_mrc (proj_stack_imag, opts->fname_out, NULL); }  CATCH_EXIT_FAIL (_e);
// 
  // printf ("\n\nProjections written to %s\n\n", opts->fname_out);
  // 
  // /* Delete all objects before exiting */
  LandwOpts_free (&opts);
  EtParams_free (&et_params);
  FwdParams_free (&fwd_params);
  tiltangles_free (&tilts);
  // gfunc3_free (&volume);
  gfunc3_free (&proj_stack);
  // gfunc3_free (&proj_stack_imag);
  
  return EXIT_SUCCESS;
}
