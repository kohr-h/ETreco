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
#include <math.h>
#include <complex.h>

#include "ETreco.h"


int main(int argc, char **argv)
{
  CEXCEPTION_T _e = EXC_NONE;

  int i, k;
  float theta_last = 0.0, theta_cur = 0.0, theta_next = 0.0, *weights = NULL;

  idx3 res_shp;
  vec3 angles;

  LandwOpts *opts = NULL;
  EtParams *et_params = NULL;
  LandwParams *landw_params = NULL;
  tiltangles *tilts = NULL;
  gfunc3 *volume = NULL, *volume_r = NULL, *proj_stack = NULL, *cur_img = NULL, *residual = NULL;
  
  Try { 
    opts = new_LandwOpts (); 
    et_params = new_EtParams (); 
    landw_params = new_LandwParams (); 
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    volume_r = new_gfunc3 ();
    proj_stack = new_gfunc3 ();
    cur_img = new_gfunc3 ();
    residual = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);


  /* Assign options, parameters, tiltangles and init gfunc structures */
  Try { LandwOpts_assign_from_args (opts, argc, argv); } CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    LandwParams_assign_from_file (landw_params, opts->fname_params); 
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
    gfunc3_init_mrc (proj_stack, opts->fname_in, NULL, NULL);
  }  CATCH_EXIT_FAIL (_e);


  /* Check consistency in number of tilts */
  Try {
    if (tilts->ntilts != proj_stack->shape[2])
      {
        EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Number of tilts (%d) must agree with stack "
        "size (%d).", tilts->ntilts, proj_stack->shape[2]);
      }
  } CATCH_EXIT_FAIL (_e);
    
  /* Apply parameters to stack (fix detector pixel size) */
  Try { LandwParams_apply_to_proj_stack (landw_params, proj_stack); }  CATCH_EXIT_FAIL (_e);

  /* Initialize, shift and scale the volume if necessary */
  Try { 
    gfunc3_init (volume, NULL, landw_params->vol_csize, landw_params->vol_shape, COMPLEX);
  }  CATCH_EXIT_FAIL (_e);
  Try {
    LandwParams_apply_to_volume (landw_params, volume);
    gfunc3_scale_grid (volume, et_params->magnification);
  }  CATCH_EXIT_FAIL (_e);


  /* Initialize the residual (only one 2D image) */
  idx3_copy (res_shp, proj_stack->shape);
  res_shp[2] = 1;
  Try { 
    gfunc3_init (residual, proj_stack->x0, proj_stack->csize, res_shp, REAL); 
  } CATCH_EXIT_FAIL (_e);


  /* Print a summary befor starting the real work */
  LandwOpts_print (opts);
  EtParams_print (et_params);
  LandwParams_print (landw_params);
  gfunc3_print_grid (volume, "volume grid:");
  gfunc3_print_grid (proj_stack, "Projection stack grid:");


  /* Normalize the stack and make it complex */
  if (normalize_flag)
    {
      printf ("\n\nNormalizing...\n\n");
      for (i = 0; i < tilts->ntilts; i++)
        {
          Try { gfunc3_set_stack_pointer (cur_img, proj_stack, i); }  CATCH_EXIT_FAIL (_e);
          Try { 
            histogram_normalization (cur_img, opts->bg_patch_ix0, opts->bg_patch_shape);
          } CATCH_EXIT_FAIL (_e);
        }
    }
  Try { gfunc3_real2complex (proj_stack); }  CATCH_EXIT_FAIL (_e);


  /* Compute the integration weights and store them (no endpoint correction) */
  Try { weights = (float *) ali16_malloc (tilts->ntilts * sizeof (float)); }  CATCH_EXIT_FAIL (_e);
  Try { tiltangles_get_angles (tilts, angles, 0); }  CATCH_EXIT_FAIL (_e);
  theta_cur = angles[1];
  for (i = 0; i < tilts->ntilts; i++)
    {
      Try { gfunc3_set_stack_pointer (cur_img, proj_stack, i); }  CATCH_EXIT_FAIL (_e);

      /* Get next tilt angles, compute integration weight and scale image accordingly */
      if (i != tilts->ntilts - 1)
        {
          Try { tiltangles_get_angles (tilts, angles, i + 1); }  CATCH_EXIT_FAIL (_e);
          theta_next = angles[1];
        }

      if (i == 0)
        weights[i] = (theta_next - theta_cur) / M_PI;
      else if (i == tilts->ntilts - 1)
        weights[i] = (theta_cur - theta_last) / M_PI;
      else
        weights[i] = (theta_next - theta_last) / (2.0 * M_PI);
      
      if (invert_contrast_flag)
        weights[i] = -weights[i];
    }


  /* Compute the starting point for the iteration */
  printf ("\n\nComputing iteration starting point...\n\n");
  for (i = 0; i < tilts->ntilts; i++)
    {
      Try { gfunc3_set_stack_pointer (cur_img, proj_stack, i); }  CATCH_EXIT_FAIL (_e);
      Try { tiltangles_get_angles (tilts, angles, i); }  CATCH_EXIT_FAIL (_e);        
      Try { 
        et_scattering_adjoint_single_axis (cur_img, angles[1], landw_params->tilt_axis, 
          et_params, volume, opts->model, weights[i]);
      } CATCH_EXIT_FAIL (_e);
    }


  volume_r = gfunc3_realpart (volume, NULL);
  temp_mrc_out (volume_r, "backproj_r", 0);

  volume_r = gfunc3_imagpart (volume, volume_r);
  temp_mrc_out (volume_r, "backproj_i", 0);

  /* TODO: continue here */
  
  /* Further steps:
   * - Get initial guess
   * - back-project
   * - Start iteration
   */
  
// 
// 
  // /* Scale and shift the volume if necessary; transfer to complex and add imaginary part */
  // Try {
    // LandwParams_apply_to_volume (landw_params, volume);
    // gfunc3_scale_grid (volume, et_params->magnification);
    // if (volume->type == REAL)
      // gfunc3_scale (volume, 1.0 + 0.1 * I);  /* TODO: use acr from config file */
  // }  CATCH_EXIT_FAIL (_e);
// 
// 
  // /* Initialize the stack */
  // Try {
    // idx3_copy (proj_stack->shape, landw_params->detector_shape);
    // proj_stack->shape[2] = tilts->ntilts;
    // gfunc3_init (proj_stack, NULL, landw_params->detector_px_size, proj_stack->shape, COMPLEX);
  // } CATCH_EXIT_FAIL (_e);
// 
// 
  // /* Print a summary of everything before starting */
  // LandwOpts_print (opts);
  // EtParams_print (et_params);
  // LandwParams_print (landw_params);
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
  LandwParams_free (&landw_params);
  tiltangles_free (&tilts);
  gfunc3_free (&volume);
  gfunc3_free (&proj_stack);
  gfunc3_free (&residual);
  free (weights);
  
  return EXIT_SUCCESS;
}
