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

  int i;
  int last_index;
  FILE *fp = NULL;
  vec3 angles_deg;

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
    proj_stack_imag = new_gfunc3 ();
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


  /* Set some parameters and resolve conflicts */
  if (opts->num_images == 0)  /* Not specified by option, thus taken from data */
    opts->num_images = tilts->ntilts - opts->start_index;

  last_index = opts->start_index + opts->num_images - 1;

  Try {
    if (last_index >= tilts->ntilts)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Image indices (max: %d) must be "
      "smaller\n than the number of tiltangles (%d).\n", last_index, tilts->ntilts);
  } CATCH_EXIT_FAIL (_e);

  /* Initialize the stack */
  Try {
    idx3_copy (proj_stack->shape, fwd_params->detector_shape);
    proj_stack->shape[2] = opts->num_images;
    gfunc3_init (proj_img, NULL, fwd_params->detector_px_size, proj_stack->shape, COMPLEX);
  } CATCH_EXIT_FAIL (_e);

  
  /* TODO: continue here */

  /* Initialize the stack on the disk */
  Try { proj_img_imag = gfunc3_imagpart (proj_img, NULL); }  CATCH_EXIT_FAIL (_e);
  Try { gfunc3_to_mrc (proj_img_imag, opts->fname_out, &fp); } CATCH_EXIT_FAIL (_e);


  /* Swap x and z axes of volume due to different ordering in NFFT */
  // Try { gfunc3_swapxz (volume); }  CATCH_EXIT_FAIL (_e);


  /* Print a summary of everything before starting */
  FwdOpts_print (opts);
  EtParams_print (et_params);
  FwdParams_print (fwd_params);
  gfunc3_print_grid (volume, "volume grid:");
  gfunc3_print_grid (proj_img, "Projection image grid:");


  for (i = opts->start_index; i <= last_index; i++)
    {
      if (verbosity_level >= VERB_LEVEL_NORMAL)
        printf ("Image %3d of %3d\n", i - opts->start_index + 1, opts->num_images);

      tiltangles_get_angles (tilts, angles_deg, i);
      Try {
        et_scattering_projection (volume, angles_deg, et_params, proj_img, opts->model);
      } CATCH_EXIT_FAIL (_e);
      
      Try { proj_img_imag = gfunc3_imagpart (proj_img, proj_img_imag); }  CATCH_EXIT_FAIL (_e);
      
      Try { gfunc3_write_to_stack (proj_img_imag, fp, i - opts->start_index); }  CATCH_EXIT_FAIL (_e);
    }

  printf ("\n\nProjections written to %s.\n\n", opts->fname_out);
  
  /* Delete all objects before exiting */
  FwdOpts_free (&opts);
  EtParams_free (&et_params);
  FwdParams_free (&fwd_params);
  tiltangles_free (&tilts);
  gfunc3_free (&volume);
  gfunc3_free (&proj_img);
  fclose (fp);
  
  return EXIT_SUCCESS;
}