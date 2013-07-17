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

/* Gnulib headers */


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
  vec3 angles, ones = {1.0, 1.0, 1.0};
  FILE *fp = NULL;
  
  vfunc vf_rk, vf_lambda, vf_recip_mtf;
  RecParams *rec_p;
  gfunc3 *proj_image, *rk, *recip_mtf, *volume;
  tiltangles *tilts;

  Try
  {
    /* Initialize options from command line */
    opt_data = new_OptionData ();
    OptionData_assign_from_args (opt_data, argc, argv);
    OptionData_print (opt_data);
  
    /* Initialize reconstruction parameters */
    rec_p = new_RecParams ();
    RecParams_assign_from_OptionData (rec_p, opt_data);
    RecParams_print (rec_p);
  
    proj_image = new_gfunc3 ();
    gfunc3_init_mrc (proj_image, opt_data->fname_in, &fp, &nz, STACK);
  
    if (rec_p->detector_px_size[0] != 0.0)  /* Detector pixel size is set in config file */
      gfunc3_set_csize (proj_image, rec_p->detector_px_size);
    
    if (opt_data->num_images == 0)  /* Not specified by option, so taken from data */
      opt_data->num_images = nz - opt_data->start_index;
  
    last_index = opt_data->start_index + opt_data->num_images - 1;
  
    if (verbosity_level >= VERB_LEVEL_NORMAL)
      gfunc3_print_grid (proj_image, "Data grid");
    
    /* Initialize tilt angles */
    tilts = new_tiltangles ();
    tiltangles_assign_from_file (tilts, opt_data->fname_tiltangles);
    
    if (last_index >= tilts->ntilts)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Image indices (max: %d) from input stack must be "
      "smaller\n than the number of tiltangles (%d).\n", last_index, tilts->ntilts);
  
    volume = new_gfunc3 ();
    gfunc3_init (volume, NULL, rec_p->vol_csize, rec_p->vol_shape, REAL);
    gfunc3_scale_grid (volume, rec_p->magnification);
    if (verbosity_level >= VERB_LEVEL_NORMAL)
      gfunc3_print_grid (volume, "volume");
  
    rk  = new_gfunc3 ();
    recip_mtf = new_gfunc3 ();
  
    /* If zero, skip computation of kernel from the second image on.
     * Add conditions here that cause the reco kernel to be dependent on tilt. 
     */
    kernel_varies_with_tilt = use_lambda_flag;
    
    for (i = opt_data->start_index; i <= last_index; i++)
      {
        printf ("Image %3d of %3d\n", i - opt_data->start_index + 1, opt_data->num_images);
        
        /* Initialize and rotate image to align tilt axis with x axis */
        gfunc3_read_from_stack (proj_image, fp, i);
        image_rotation (proj_image, -rec_p->tilt_axis);
        
        if (DEBUGGING)
          temp_mrc_out (proj_image, "rotated_", i + 1);
        
        /* Normalization of the current image */
        if (normalize_flag)
          {
            histogram_normalization (proj_image, opt_data->bg_patch_ix0, opt_data->bg_patch_shape);
            if (DEBUGGING)
              temp_mrc_out (proj_image, "normalized_", i + 1);
            // probability_normalization (proj_image);
            // gfunc3_to_mrc (proj_image, "normalized2_.mrc");
          }
        
        if (invert_contrast_flag)
          gfunc3_scale (proj_image, -1.0);
        
        /* Fourier transform the image */
        fft_forward (proj_image);
        if (verbosity_level >= VERB_LEVEL_VERBOSE)
          gfunc3_print_grid (proj_image, "Image FT grid");
        if (DEBUGGING)
          temp_mrc_out (proj_image, "ft_image_", i + 1);

        /* Compute reco kernel (only if necessary) */
        if ((i == opt_data->start_index) || kernel_varies_with_tilt )
          {
            float nul[2] = {0.0f, 0.0f};
            gfunc3_init_from_foreign_grid (rk, proj_image);
            gfunc3_set_all (rk, nul);
            if (DEBUGGING)
              temp_mrc_out (rk, "ft_rk_zero_", i + 1);
            
            vfunc_init_ft_rk_single_axis_x (&vf_rk, rec_p);
            gfunc3_assign_fvals_from_vfunc (rk, &vf_rk);
  
            if (DEBUGGING)
              temp_mrc_out (rk, "ft_rk_", i + 1);
          }
  
        /* Apply kernel */
        gfunc3_mul (proj_image, rk);
        if (DEBUGGING)
          temp_mrc_out (proj_image, "multiplied_ft_rk_", i + 1);
  
        /* Initialize reciprocal detector MTF (only once) */
        if ((i == opt_data->start_index) && use_mtf_flag)
          {
            gfunc3_init_from_foreign_grid (recip_mtf, proj_image);
            gfunc3_set_csize (recip_mtf, ones);
            gfunc3_compute_xmin_xmax (recip_mtf);
  
            vfunc_init_detector_recip_mtf (&vf_recip_mtf, rec_p);
            gfunc3_assign_fvals_from_vfunc (recip_mtf, &vf_recip_mtf);
            gfunc3_set_csize (recip_mtf, proj_image->csize);
  
            if (DEBUGGING)
              temp_mrc_out (recip_mtf, "recip_mtf_", i + 1);
          }
  
        /* Apply reciprocal MTF */
        if (use_mtf_flag)
          gfunc3_mul (proj_image, recip_mtf);
  
  
        /* Apply lambda if desired */
        if (use_lambda_flag)
          {
            vfunc_init_ft_lambda (&vf_lambda, &opt_data->lambda_pow);
            gfunc3_mul_vfunc (proj_image, &vf_lambda);
          }
  
        /* Fourier transform back */
        fft_backward (proj_image);
        
        if (DEBUGGING)
          temp_mrc_out (proj_image, "filtered_", i + 1);
        
        /* Get current tilt angles and compute backprojection; Since we assume single axis geometry, 
         * use only angles[1] (theta).
         */
        tiltangles_get_angles (tilts, angles, i);
        if (verbosity_level >= VERB_LEVEL_VERBOSE)
          printf ("theta = %f degrees\n", angles[1]);
        
        xray_backprojection_sax (proj_image, angles[1], volume);
      }
  
    // gfunc3_make_nonneg (volume);
    gfunc3_to_mrc (volume, opt_data->fname_out);
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
  }
  Catch (e)
  {
    EXC_REPRINT;
  }
  
  return 0;
}
