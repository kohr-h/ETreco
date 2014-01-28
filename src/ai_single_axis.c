/*
 * ai_single_axis.c -- approximate inverse for single-axis ET data
 * 
 * Copyright 2014 Holger Kohr <kohr@num.uni-sb.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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

/* Configuration header */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* System headers */
#include <stdio.h>
#include <math.h>
#include <omp.h>


/* All ETreco headers */
#include "ETreco.h"

#define ONE_MICROMETER  1E3   // [nm]

int kernel_varies_with_tilt = 0;

int
main (int argc, char *argv[])
{
  CEXCEPTION_T _e = EXC_NONE;
  AiOpts *opts;

  int i, nz, last_index;
  float weight, theta_last = 0.0, theta_cur = 0.0, theta_next = 0.0;
  vec3 angles, ones = {1.0, 1.0, 1.0};
  FILE *fp = NULL;
  
  vfunc vf_rk, vf_lambda, vf_recip_mtf;
  EtParams *et_params;
  AiParams *ai_params;
  EtAiParWrapper pwrapper;
  gfunc3 *proj_image, *rk, *recip_mtf, *volume;
  tiltangles *tilts;

  /* Create structures */
  Try 
  { 
    opts = new_AiOpts (); 
    ai_params = new_AiParams ();
    et_params = new_EtParams ();
    proj_image = new_gfunc3 ();
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    rk  = new_gfunc3 ();
    recip_mtf = new_gfunc3 ();
  } CATCH_EXIT_FAIL (_e);

  /* Initialize command-line options, reconstruction parameters, first projection image, 
   * tilt angles and volume
   */
  Try { AiOpts_assign_from_args (opts, argc, argv); }  CATCH_EXIT_FAIL (_e);
  Try { 
    EtParams_assign_from_file (et_params, opts->fname_params); 
    AiParams_assign_from_AiOpts (ai_params, opts);
    }  CATCH_EXIT_FAIL (_e);
  Try { 
    AiParams_assign_ctftrunc_from_EtParams (ai_params, et_params);
    gfunc3_init_mrc (proj_image, opts->fname_in, &fp, &nz);
    tiltangles_assign_from_file (tilts, opts->fname_tiltangles);
  } CATCH_EXIT_FAIL (_e);

  /* Wrap the parameters */
  pwrapper.et_params = et_params;
  pwrapper.ai_params = ai_params;
  
  Try { 
    gfunc3_init (volume, NULL, ai_params->vol_csize, ai_params->vol_shape, REAL); 
  }  CATCH_EXIT_FAIL (_e);


  /* Scale and shift volume according to the parameters */
  Try { 
    AiParams_apply_to_volume (ai_params, volume); 
    gfunc3_scale_grid (volume, et_params->magnification); 
  } CATCH_EXIT_FAIL (_e);


  /* Set some parameters and resolve conflicts */
  if (opts->num_images == 0)  /* Not specified by option, thus taken from data */
    opts->num_images = nz - opts->start_index;

  last_index = opts->start_index + opts->num_images - 1;

  Try {
    if (last_index >= tilts->ntilts)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Image indices (max: %d) from input stack must be "
      "smaller\n than the number of tiltangles (%d).\n", last_index, tilts->ntilts);
  } CATCH_EXIT_FAIL (_e);

  /* If zero, skip computation of kernel from the second image on.
   * Add conditions here that cause the reco kernel to be dependent on tilt. 
   */
  kernel_varies_with_tilt = use_lambda_flag;


  /* Print summary of everything before the real work starts */
  AiOpts_print (opts);
  EtParams_print (et_params);
  AiParams_print (ai_params);
  if (verbosity_level >= VERB_LEVEL_NORMAL)
    gfunc3_print_grid (volume, "Volume grid:");

  /* Initialize current theta for the loop */
  Try { tiltangles_get_angles (tilts, angles, 0); }  CATCH_EXIT_FAIL (_e);
  theta_cur = angles[1];

  
  for (i = opts->start_index; i <= last_index; i++)
    {
      /* Initialize image and normalize if desired */
      Try { 
        gfunc3_read_from_stack (proj_image, fp, i);
        AiParams_apply_to_proj_image (ai_params, proj_image); 
      } CATCH_EXIT_FAIL (_e);

      if ((verbosity_level >= VERB_LEVEL_VERBOSE) && (i == opts->start_index))
        gfunc3_print_grid (proj_image, "Data grid:");

      if (verbosity_level >= VERB_LEVEL_NORMAL)
        printf ("Image %3d of %3d\n", i - opts->start_index + 1, opts->num_images);
     
      if (normalize_flag)
        {
          Try {
            histogram_normalization (proj_image, opts->bg_patch_ix0, opts->bg_patch_shape);
          } CATCH_EXIT_FAIL (_e);
          if (DEBUGGING)
            {
              Try { temp_mrc_out (proj_image, "normalized_", i + 1); }  CATCH_EXIT_FAIL (_e);
          // probability_normalization (proj_image);
          // temp_mrc_out (proj_image, "normalized2_", i + 1);
            }
        }

      /* Get next tilt angles, compute integration weight and scale image accordingly */
      if (i != last_index)
        {
          Try { tiltangles_get_angles (tilts, angles, i + 1); }  CATCH_EXIT_FAIL (_e);
          theta_next = angles[1];
        }

      if (i == opts->start_index)
        weight = (theta_next - theta_cur) / 2.0;
      else if (i == last_index)
        weight = (theta_cur - theta_last) / 2.0;
      else
        weight = (theta_next - theta_last) / 2.0;

      if (verbosity_level >= VERB_LEVEL_VERBOSE)
        {
          printf ("theta = %f degrees\n", theta_cur);
          printf ("weight = %f\n", weight);
        }

      if (invert_contrast_flag)
        weight = -weight;
        
      Try { gfunc3_scale (proj_image, weight); }  CATCH_EXIT_FAIL (_e);

      
      /* Rotate image to align tilt axis with x or y axis */
      Try { image_rotation (proj_image, -ai_params->tilt_axis_rotation); }  CATCH_EXIT_FAIL (_e);

      if (DEBUGGING)
        {
          Try { temp_mrc_out (proj_image, "rotated_", i + 1); }  CATCH_EXIT_FAIL (_e);
        }

      /* Fourier transform the image */
      Try { fft_forward (proj_image); }   CATCH_EXIT_FAIL (_e);
      if ((verbosity_level >= VERB_LEVEL_VERBOSE) && (i == opts->start_index))
        gfunc3_print_grid (proj_image, "Image FT grid");
      if (DEBUGGING)
        {
          Try { temp_mrc_out (proj_image, "ft_image_", i + 1); }  CATCH_EXIT_FAIL (_e);
        }

      /* Compute reco kernel (only if necessary) */
      if ((i == opts->start_index) || kernel_varies_with_tilt )
        {
          Try { gfunc3_init_from_foreign_grid (rk, proj_image); }  CATCH_EXIT_FAIL (_e);
          gfunc3_set_all (rk, 0.0);

          Try { vfunc_init_ft_rk_single_axis (&vf_rk, &pwrapper); }  CATCH_EXIT_FAIL (_e);
          Try { gfunc3_assign_fvals_from_vfunc (rk, &vf_rk); }  CATCH_EXIT_FAIL (_e);

          if (DEBUGGING)
            {
              Try { temp_mrc_out (rk, "ft_rk_", i + 1); }  CATCH_EXIT_FAIL (_e);
            }
        }

      /* Apply kernel */
      Try { gfunc3_mul (proj_image, rk); }  CATCH_EXIT_FAIL (_e);
      if (DEBUGGING)
        {
          Try { temp_mrc_out (proj_image, "multiplied_ft_rk_", i + 1); }  CATCH_EXIT_FAIL (_e);
        }

      /* Initialize reciprocal detector MTF (only once) */
      if ((i == opts->start_index) && use_mtf_flag)
        {
          Try { gfunc3_init_from_foreign_grid (recip_mtf, proj_image); }  CATCH_EXIT_FAIL (_e);
          gfunc3_set_csize (recip_mtf, ones);
          gfunc3_compute_xmin_xmax (recip_mtf);
          gfunc3_set_all (recip_mtf, 0.0);

          Try { vfunc_init_detector_recip_mtf (&vf_recip_mtf, et_params); }  CATCH_EXIT_FAIL (_e);
          Try { gfunc3_assign_fvals_from_vfunc (recip_mtf, &vf_recip_mtf); }  CATCH_EXIT_FAIL (_e);
          gfunc3_set_csize (recip_mtf, proj_image->csize);

          if (DEBUGGING)
            {
              Try { temp_mrc_out (recip_mtf, "recip_mtf_", i + 1); }  CATCH_EXIT_FAIL (_e);
            }
        }

      /* Apply reciprocal MTF */
      if (use_mtf_flag)
        {
          Try { gfunc3_mul (proj_image, recip_mtf); }  CATCH_EXIT_FAIL (_e);
          if (DEBUGGING)
            {
              Try { temp_mrc_out (proj_image, "multiplied_mtf_", i + 1); }  CATCH_EXIT_FAIL (_e);
            }
        }  

      /* TODO: implement and distinguish lambda_x and lambda_y */
      /* Apply lambda if desired */
      if (use_lambda_flag)
        {
          Try { vfunc_init_ft_lambda (&vf_lambda, &opts->lambda_pow); }  CATCH_EXIT_FAIL (_e);
          Try { gfunc3_mul_vfunc (proj_image, &vf_lambda); }  CATCH_EXIT_FAIL (_e);
        }

      /* Fourier transform back */
      Try { fft_backward (proj_image); }  CATCH_EXIT_FAIL (_e);
      
      if (DEBUGGING)
        {
          Try { temp_mrc_out (proj_image, "filtered_", i + 1); }  CATCH_EXIT_FAIL (_e);
        }

      /* Compute backprojection */
      Try { 
        xray_backprojection_single_axis (proj_image, theta_cur, ai_params->tilt_axis,
          ai_params->tilt_axis_par_shift_px, volume); 
      } CATCH_EXIT_FAIL (_e);
      // Try { tiltangles_get_angles (tilts, angles, i); }  CATCH_EXIT_FAIL (_e);
      // Try 
      // { 
        // xray_backprojection (proj_image, angles, volume); 
      // } CATCH_EXIT_FAIL (_e);
    
      
      /* Update thetas */
      theta_last = theta_cur;
      theta_cur  = theta_next;
    }

  Try { gfunc3_to_mrc (volume, opts->fname_out, NULL); }  CATCH_EXIT_FAIL (_e);
  printf ("Reconstructed volume written to %s\n", opts->fname_out);

  if (fp)
    fclose (fp);
  AiOpts_free (&opts);
  EtParams_free (&et_params);
  AiParams_free (&ai_params);
  tiltangles_free (&tilts);
  gfunc3_free (&proj_image);
  gfunc3_free (&rk);
  gfunc3_free (&recip_mtf);
  gfunc3_free (&volume);
  
  return 0;
}
