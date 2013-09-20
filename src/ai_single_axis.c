/*
 * ai_single_axis.c -- approximate inverse for single-axis ET data
 * 
 * Copyright 2013 Holger Kohr <kohr@num.uni-sb.de>
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
  CEXCEPTION_T e = EXC_NONE;
  OptionData *opt_data;

  int i, nz, last_index;
  float weight, theta_last = 0.0, theta_cur = 0.0, theta_next = 0.0;
  vec3 angles, ones = {1.0, 1.0, 1.0};
  FILE *fp = NULL;
  
  vfunc vf_rk, vf_lambda, vf_recip_mtf;
  RecParams *rec_p;
  gfunc3 *proj_image, *rk, *recip_mtf, *volume;
  tiltangles *tilts;

  /* Create structures */
  Try 
  { 
    opt_data = new_OptionData (); 
    rec_p = new_RecParams ();
    proj_image = new_gfunc3 ();
    tilts = new_tiltangles ();
    volume = new_gfunc3 ();
    rk  = new_gfunc3 ();
    recip_mtf = new_gfunc3 ();
  }  CATCH_EXIT_FAIL (e);

  /* Initialize command-line options, reconstruction parameters, first projection image, 
   * tilt angles and volume
   */
  Try 
  { 
    OptionData_assign_from_args (opt_data, argc, argv); 
    RecParams_assign_from_OptionData (rec_p, opt_data);
    gfunc3_init_mrc (proj_image, opt_data->fname_in, &fp, &nz, STACK);
    tiltangles_assign_from_file (tilts, opt_data->fname_tiltangles);
    gfunc3_init (volume, NULL, rec_p->vol_csize, rec_p->vol_shape, REAL);
  } CATCH_EXIT_FAIL (e);

  /* Scale and shift volume according to the reco parameters */
  Try 
  { 
    gfunc3_scale_grid (volume, rec_p->magnification); 
    RecParams_apply_to_volume (rec_p, volume); 
  }  CATCH_EXIT_FAIL (e);



  /* Set some parameters and resolve conflicts */
  if (opt_data->num_images == 0)  /* Not specified by option, so taken from data */
    opt_data->num_images = nz - opt_data->start_index;

  last_index = opt_data->start_index + opt_data->num_images - 1;

  Try {
  if (last_index >= tilts->ntilts)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Image indices (max: %d) from input stack must be "
    "smaller\n than the number of tiltangles (%d).\n", last_index, tilts->ntilts);
  }  CATCH_EXIT_FAIL (e);

  /* If zero, skip computation of kernel from the second image on.
   * Add conditions here that cause the reco kernel to be dependent on tilt. 
   */
  kernel_varies_with_tilt = use_lambda_flag;


  /* Print summary of everything before the real work starts */
  OptionData_print (opt_data);
  RecParams_print (rec_p);
  if (verbosity_level >= VERB_LEVEL_NORMAL)
    gfunc3_print_grid (volume, "Volume grid:");


  /* Initialize current theta for the loop */
  Try { tiltangles_get_angles (tilts, angles, 0); }  CATCH_EXIT_FAIL (e);
  theta_cur = angles[1];

  
  for (i = opt_data->start_index; i <= last_index; i++)
    {
      /* Initialize image and normalize if desired */
      Try 
      { 
        gfunc3_read_from_stack (proj_image, fp, i);
        RecParams_apply_to_proj_image (rec_p, proj_image); 
      }  CATCH_EXIT_FAIL (e);

      if ((verbosity_level >= VERB_LEVEL_VERBOSE) && (i == opt_data->start_index))
        gfunc3_print_grid (proj_image, "Data grid:");

      if (verbosity_level >= VERB_LEVEL_NORMAL)
        printf ("Image %3d of %3d\n", i - opt_data->start_index + 1, opt_data->num_images);
     
      if (normalize_flag)
        {
          Try {
          histogram_normalization (proj_image, opt_data->bg_patch_ix0, opt_data->bg_patch_shape);
          }  CATCH_EXIT_FAIL (e);
          if (DEBUGGING)
            temp_mrc_out (proj_image, "normalized_", i + 1);
          // probability_normalization (proj_image);
          // gfunc3_to_mrc (proj_image, "normalized2_.mrc");
        }

      /* Get next tilt angles, compute integration weight and scale image accordingly */
      if (i != last_index)
        {
          Try { tiltangles_get_angles (tilts, angles, i + 1); }  CATCH_EXIT_FAIL (e);
          theta_next = angles[1];
        }

      if (i == opt_data->start_index)
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
        
      Try { gfunc3_scale (proj_image, weight); }  CATCH_EXIT_FAIL (e);

      
      /* Rotate image to align tilt axis with x axis */
      Try { image_rotation (proj_image, -rec_p->tilt_axis_rotation); }  CATCH_EXIT_FAIL (e);

      if (DEBUGGING)
        temp_mrc_out (proj_image, "rotated_", i + 1);

      /* Fourier transform the image */
      Try { fft_forward (proj_image); }   CATCH_EXIT_FAIL (e);
      if ((verbosity_level >= VERB_LEVEL_VERBOSE) && (i == opt_data->start_index))
        gfunc3_print_grid (proj_image, "Image FT grid");
      if (DEBUGGING)
        temp_mrc_out (proj_image, "ft_image_", i + 1);

      /* Compute reco kernel (only if necessary) */
      if ((i == opt_data->start_index) || kernel_varies_with_tilt )
        {
          float nul[2] = {0.0f, 0.0f};
          Try { gfunc3_init_from_foreign_grid (rk, proj_image); }  CATCH_EXIT_FAIL (e);
          gfunc3_set_all (rk, nul);

          Try { vfunc_init_ft_rk_single_axis_x (&vf_rk, rec_p); }  CATCH_EXIT_FAIL (e);
          Try { gfunc3_assign_fvals_from_vfunc (rk, &vf_rk); }  CATCH_EXIT_FAIL (e);

          if (DEBUGGING)
            temp_mrc_out (rk, "ft_rk_", i + 1);
        }

      /* Apply kernel */
      Try { gfunc3_mul (proj_image, rk); }  CATCH_EXIT_FAIL (e);
      if (DEBUGGING)
        temp_mrc_out (proj_image, "multiplied_ft_rk_", i + 1);

      /* Initialize reciprocal detector MTF (only once) */
      if ((i == opt_data->start_index) && use_mtf_flag)
        {
          float nul[2] = {0.0f, 0.0f};
          Try { gfunc3_init_from_foreign_grid (recip_mtf, proj_image); }  CATCH_EXIT_FAIL (e);
          gfunc3_set_csize (recip_mtf, ones);
          gfunc3_compute_xmin_xmax (recip_mtf);
          gfunc3_set_all (recip_mtf, nul);

          Try { vfunc_init_detector_recip_mtf (&vf_recip_mtf, rec_p); }  CATCH_EXIT_FAIL (e);
          Try { gfunc3_assign_fvals_from_vfunc (recip_mtf, &vf_recip_mtf); }  CATCH_EXIT_FAIL (e);
          gfunc3_set_csize (recip_mtf, proj_image->csize);

          if (DEBUGGING)
            temp_mrc_out (recip_mtf, "recip_mtf_", i + 1);
        }

      /* Apply reciprocal MTF */
      if (use_mtf_flag)
        {
          Try { gfunc3_mul (proj_image, recip_mtf); }  CATCH_EXIT_FAIL (e);
          if (DEBUGGING)
            temp_mrc_out (proj_image, "multiplied_mtf_", i + 1);
        }  

      /* Apply lambda if desired */
      if (use_lambda_flag)
        {
          Try { vfunc_init_ft_lambda (&vf_lambda, &opt_data->lambda_pow); }  CATCH_EXIT_FAIL (e);
          Try { gfunc3_mul_vfunc (proj_image, &vf_lambda); }  CATCH_EXIT_FAIL (e);
        }

      /* Fourier transform back */
      Try { fft_backward (proj_image); }  CATCH_EXIT_FAIL (e);
      
      if (DEBUGGING)
        temp_mrc_out (proj_image, "filtered_", i + 1);
      
      /* Compute backprojection */
      Try { xray_backprojection_sax (proj_image, theta_cur, volume); }  CATCH_EXIT_FAIL (e);
      // Try { tiltangles_get_angles (tilts, angles, i); }  CATCH_EXIT_FAIL (e);
      // Try { xray_backprojection (proj_image, angles, volume); }  CATCH_EXIT_FAIL (e);
      
      /* Update thetas */
      theta_last = theta_cur;
      theta_cur  = theta_next;
    }

  Try { gfunc3_to_mrc (volume, opt_data->fname_out); }  CATCH_EXIT_FAIL (e);
  printf ("Reconstructed volume written to %s\n", opt_data->fname_out);

  if (fp)
    fclose (fp);
  OptionData_free (&opt_data);
  RecParams_free (&rec_p);
  tiltangles_free (&tilts);
  gfunc3_free (&proj_image);
  gfunc3_free (&rk);
  gfunc3_free (&recip_mtf);
  gfunc3_free (&volume);
  
  return 0;
}
