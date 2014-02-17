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

  vec3 angles;

  LandwOpts *opts = NULL;
  EtParams *et_params = NULL;
  LandwParams *landw_params = NULL;
  tiltangles *tilts = NULL;
  gfunc3 *volume = NULL, *proj_stack = NULL, *residual_stack = NULL;
  gfunc3 *cur_img1 = NULL, *cur_img2 = NULL;
  
  Try { 
    opts = new_LandwOpts (); 
    et_params = new_EtParams (); 
    landw_params = new_LandwParams (); 
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    proj_stack = new_gfunc3 ();
    cur_img1 = new_gfunc3 ();
    cur_img2 = new_gfunc3 ();
    residual_stack = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);


  /* Assign options, parameters, tiltangles and init gfunc structures */
  Try { LandwOpts_assign_from_args (opts, argc, argv); } CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    LandwParams_assign_from_file (landw_params, opts->fname_params); 
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
    gfunc3_init_mrc (proj_stack, opts->fname_in, NULL, NULL);
  } CATCH_EXIT_FAIL (_e);


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
  } CATCH_EXIT_FAIL (_e);


  /* Print a summary befor starting the real work */
  LandwOpts_print (opts);
  EtParams_print (et_params);
  LandwParams_print (landw_params);
  gfunc3_print_grid (volume, "volume grid:");
  gfunc3_print_grid (proj_stack, "Projection stack grid:");


  /* Normalize the stack */
  if (normalize_flag)
    {
      printf ("\n\nNormalizing...\n\n");
      for (i = 0; i < tilts->ntilts; i++)
        {
          Try { gfunc3_set_stack_pointer (cur_img1, proj_stack, i); }  CATCH_EXIT_FAIL (_e);
          Try { 
            histogram_normalization (cur_img1, opts->bg_patch_ix0, opts->bg_patch_shape);
          } CATCH_EXIT_FAIL (_e);
        }
      Try { gfunc3_real2complex (proj_stack); }  CATCH_EXIT_FAIL (_e);
    }


  /* Initialize the residual stack */
  Try { gfunc3_init_from_foreign_grid (residual_stack, proj_stack); }  CATCH_EXIT_FAIL (_e);



  /* Compute the integration weights and store them (no endpoint correction) */
  Try { weights = (float *) ali16_malloc (tilts->ntilts * sizeof (float)); }  CATCH_EXIT_FAIL (_e);
  Try { tiltangles_get_angles (tilts, angles, 0); }  CATCH_EXIT_FAIL (_e);
  theta_cur = angles[1];
  for (i = 0; i < tilts->ntilts; i++)
    {
      Try { gfunc3_set_stack_pointer (cur_img1, proj_stack, i); }  CATCH_EXIT_FAIL (_e);

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
      Try { gfunc3_set_stack_pointer (cur_img1, proj_stack, i); }  CATCH_EXIT_FAIL (_e);
      Try { tiltangles_get_angles (tilts, angles, i); }  CATCH_EXIT_FAIL (_e);        
      Try { 
        et_scattering_adjoint_single_axis (cur_img1, angles[1], landw_params->tilt_axis, 
          et_params, volume, opts->model, weights[i]);
      } CATCH_EXIT_FAIL (_e);
    }


  /* Iterate */
  
  for (k = 0; k < opts->max_iter; k++)
    {
      /* Project current guess */
      Try { 
        et_scattering_projection_atonce (volume, tilts, et_params, residual_stack, PROJ_ASSUMPTION);
      } CATCH_EXIT_FAIL (_e);
      
      if (DEBUGGING)
        temp_mrc_out (residual_stack, "reproj_stk_", k + 1);
      
      /* Now image per image */
      
      for (i = 0; i < tilts->ntilts; i++)
        {
          /* Compute residual */
          Try { 
            gfunc3_set_stack_pointer (cur_img1, proj_stack, i);
            gfunc3_set_stack_pointer (cur_img2, residual_stack, i); 
          } CATCH_EXIT_FAIL (_e);
          
          Try { gfunc3_axpy (-1, cur_img2, cur_img1); }  CATCH_EXIT_FAIL (_e);
          
          /* Backproject and add to volume */
          Try { 
            et_scattering_adjoint_single_axis (cur_img2, angles[1], landw_params->tilt_axis, 
              et_params, volume, opts->model, opts->relax_param * weights[i]);
          } CATCH_EXIT_FAIL (_e);
        }
      
      if (DEBUGGING)
        {
          temp_mrc_out (residual_stack, "residual_stk_", k + 1);
          temp_mrc_out (volume, "volume_", k + 1);
        }
      
      /* TODO: Introduce stopping criterion */
    }
  

  /* Delete all objects before exiting */
  LandwOpts_free (&opts);
  EtParams_free (&et_params);
  LandwParams_free (&landw_params);
  tiltangles_free (&tilts);
  gfunc3_free (&volume);
  gfunc3_free (&proj_stack);
  gfunc3_free (&residual_stack);
  free (weights);
  
  return EXIT_SUCCESS;
}
