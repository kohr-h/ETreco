/*
 * et_operators.c -- operators specific for ET
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include "CException.h"

#include "vec3.h"
#include "misc.h"

#include "tiltangles.h"

#include "fft.h"
#include "fft_private.h"
#include "gfunc3.h"
#include "gfunc3_private.h"
#include "mrc.h"
#include "operators_private.h"

#include "et_params.h"
#include "et_operators.h"
#include "et_vfuncs.h"

/*-------------------------------------------------------------------------------------------------*/

float *
ewald_sphere_freqs (gfunc3 const *ft_proj_img_grid, vec3 const normal_angles_deg, float wave_number)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  float *freqs = NULL;

  vec3 p;
  int ix, iy;
  float sin_psi, cos_psi, sin_theta, cos_theta, sin_phi, cos_phi;
  float *cur_freq;
  double tmp, v;

  Try { 
    freqs = (float *) ali16_malloc (3 * ft_proj_img_grid->ntotal * sizeof (float));
  } CATCH_RETURN (_e, NULL);
  
  sin_psi   = sinf (normal_angles_deg[0] * ONE_DEGREE);
  cos_psi   = cosf (normal_angles_deg[0] * ONE_DEGREE);
  sin_theta = sinf (normal_angles_deg[1] * ONE_DEGREE);
  cos_theta = cosf (normal_angles_deg[1] * ONE_DEGREE);
  sin_phi   = sinf (normal_angles_deg[2] * ONE_DEGREE);
  cos_phi   = cosf (normal_angles_deg[2] * ONE_DEGREE);
  
  vec3_copy (p, ft_proj_img_grid->xmin);

  for (iy = 0, cur_freq = freqs; iy < ft_proj_img_grid->shape[1]; iy++)
    {
      for (ix = 0; ix < ft_proj_img_grid->shape[0]; ix++, cur_freq += 3)
        {
          /* The same as in perp_plane_freqs, but with an additional summand consisting of
           * (k-v(\xi))\omega, where k is the wave number, v(xi) = sqrt(k^2-v^2), and \omega
           * is the third column of the transpose of the rotation matrix.
           */
          tmp = p[0] * p[0] + p[1] * p[1];
          v = sqrt (wave_number * wave_number - tmp);
          tmp = wave_number - v;
          
          cur_freq[0] =   p[0] * ( cos_phi * cos_psi - sin_phi * cos_theta * sin_psi)
                        + p[1] * (-cos_phi * sin_psi - sin_phi * cos_theta * cos_psi)
                        + tmp  * ( sin_phi * sin_theta);
          cur_freq[1] =   p[0] * ( sin_phi * cos_psi + cos_phi * cos_theta * sin_psi)
                        + p[1] * (-sin_phi * sin_psi + cos_phi * cos_theta * cos_psi)
                        + tmp  * (-cos_phi * sin_theta);
          cur_freq[2] =   p[0] * ( sin_theta * sin_psi)
                        + p[1] * ( sin_theta * cos_psi)
                        + tmp  * ( cos_theta);
          
          p[0] += ft_proj_img_grid->csize[0];
        }
      p[0]  = ft_proj_img_grid->xmin[0];
      p[1] += ft_proj_img_grid->csize[1];
    }

  return freqs;
}

/*-------------------------------------------------------------------------------------------------*/

float *
ewald_sphere_stack_freqs (gfunc3 const *ft_proj_img_grid, tiltangles const *tilts, float wave_number)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int i;
  size_t j;
  
  vec3 angles_deg;
  float *freqs = NULL, *allfreqs = NULL, *cur_freq;
  
  CAPTURE_NULL (ft_proj_img_grid, NULL);
  CAPTURE_NULL (tilts, NULL);
  
  /* Alloc space for the huge frequency array */
  Try { 
    allfreqs = (float *) ali16_malloc (3 * ft_proj_img_grid->ntotal * tilts->ntilts * sizeof (float));
  } CATCH_RETURN (_e, NULL);
  
  cur_freq = allfreqs;
  for (i = 0; i < tilts->ntilts; i++)
    {
      /* Compute frequencies for the i'th sphere and copy them to the huge array. */
      Try { tiltangles_get_angles (tilts, angles_deg, i); }  CATCH_RETURN (_e, NULL);
      Try { 
        freqs = ewald_sphere_freqs (ft_proj_img_grid, angles_deg, wave_number);
      } CATCH_RETURN (_e, NULL);
      
      for (j = 0; j < 3 * ft_proj_img_grid->ntotal; j++)
        *(cur_freq++) = freqs[j];
    }
  
  free (freqs);
  
  return allfreqs;
}

/*-------------------------------------------------------------------------------------------------*/

void
et_scattering_projection (gfunc3 const *scatterer, vec3 const angles_deg, EtParams const *params, 
                          gfunc3 *proj_img, scattering_model sct_model)
{
  CEXCEPTION_T _e = EXC_NONE;

  float *freqs = NULL;
  vfunc ctf;
  
  CAPTURE_NULL_VOID (scatterer);
  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (angles_deg);
  if (sct_model != PROJ_ASSUMPTION)
    CAPTURE_NULL_VOID (params);
  
  GFUNC_CAPTURE_UNINIT_VOID (scatterer);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);

  if (proj_img->type == REAL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Projection image must be of COMPLEX type.");
      return;
    }

  if ( (sct_model == BORN_APPROX) && (params->wave_number == 0.0) )
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Wave number must not be 0 in BORN_APPROX model.");
      return;
    }

  /* Take the projection image to the reciprocal space */
  Try { gfunc3_grid_fwd_reciprocal (proj_img); }  CATCH_RETURN_VOID (_e);

  /* Evaluate the FT of the SCATTERER at frequencies according to the chosen model. For the
   * PROJ_ASSUMPTION model, the frequencies lie on the 2D plane perpendicular to the unit vector 
   * defined by ANGLES_DEG. For BORN_APPROX, they lie on the corresponding Ewald sphere.
   */
  if (sct_model == PROJ_ASSUMPTION)
    {
      Try { freqs = perp_plane_freqs (proj_img, angles_deg); }  CATCH_RETURN_VOID (_e);

      Try { 
        nfft_transform (scatterer, freqs, proj_img->ntotal, (float complex **) (&proj_img->fvals)); 
      } CATCH_RETURN_VOID (_e);
   }
  else if (sct_model == BORN_APPROX)
    {
      Try { 
        freqs = ewald_sphere_freqs (proj_img, angles_deg, params->wave_number); 
      } CATCH_RETURN_VOID (_e);

      Try { 
        nfft_transform (scatterer, freqs, proj_img->ntotal, (float complex **) (&proj_img->fvals)); 
      } CATCH_RETURN_VOID (_e);
    }
  
  /* FIXXME: there is something wrong with the factors. Check with math! */
  /* Multiply with CTF and scaling factors. Finally apply inverse FT */
  Try { vfunc_init_ctf (&ctf, params); } CATCH_RETURN_VOID (_e);
  Try {
    gfunc3_mul_vfunc (proj_img, &ctf);
    gfunc3_scale (proj_img, 2 * M_PI * M_SQRT2PI);
  } CATCH_RETURN_VOID (_e);

  Try { fft_backward (proj_img); } CATCH_RETURN_VOID (_e);
  
  free (freqs);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
et_scattering_projection_atonce (gfunc3 const *scatterer, tiltangles const *tilts, 
                                 EtParams const *params, gfunc3 *proj_stack, 
                                 scattering_model sct_model)
{
  CEXCEPTION_T _e = EXC_NONE;

  int i;
  float *freqs = NULL;
  idx3 proj_shp;
  vfunc ctf;
  gfunc3 *gf_tmp = new_gfunc3 ();
  
  CAPTURE_NULL_VOID (scatterer);
  CAPTURE_NULL_VOID (proj_stack);
  CAPTURE_NULL_VOID (tilts);
  if (sct_model != PROJ_ASSUMPTION)
    CAPTURE_NULL_VOID (params);
  
  GFUNC_CAPTURE_UNINIT_VOID (scatterer);
  GFUNC_CAPTURE_UNINIT_VOID (proj_stack);

  if (proj_stack->type != COMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Projection stack must be of COMPLEX type.");
      return;
    }

  if (proj_stack->shape[2] < tilts->ntilts)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Projection stack z shape too small.");
      return;
    }

  if ( (sct_model == BORN_APPROX) && (params->wave_number == 0.0) )
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Wave number must not be 0 in BORN_APPROX model.");
      return;
    }

  /* Create the reciprocal grid for one single projection */
  idx3_copy (proj_shp, proj_stack->shape);
  proj_shp[2] = 1;
  gfunc3_init_gridonly (gf_tmp, proj_stack->x0, proj_stack->csize, proj_shp, COMPLEX);
  Try { gfunc3_grid_fwd_reciprocal (gf_tmp); }  CATCH_RETURN_VOID (_e);

  /* Evaluate the FT of the SCATTERER at frequencies according to the chosen model. For the
   * PROJ_ASSUMPTION model, the frequencies lie on the union of 2D planes perpendicular to the 
   * unit vectors defined by TILTS. For BORN_APPROX, they lie on the corresponding Ewald spheres.
   */
  if (sct_model == PROJ_ASSUMPTION)
    {
      Try { freqs = perp_plane_stack_freqs (gf_tmp, tilts); }  CATCH_RETURN_VOID (_e);
    }
  else if (sct_model == BORN_APPROX)
    {
      Try { 
        freqs = ewald_sphere_stack_freqs (gf_tmp, tilts, params->wave_number); 
      } CATCH_RETURN_VOID (_e);
    }
  
  Try { 
    nfft_transform (scatterer, freqs, proj_stack->ntotal, (float complex **) (&proj_stack->fvals)); 
  } CATCH_RETURN_VOID (_e);

  Try { vfunc_init_ctf (&ctf, params); } CATCH_RETURN_VOID (_e);
  gf_tmp->is_initialized = TRUE;
  for (i = 0; i < tilts->ntilts; i++)
    {
      /* Have GF_TMP->FVALS point to the current image in the stack */
      gf_tmp->fvals = &proj_stack->fvals[i * 2 * gf_tmp->ntotal];
      
      /* Image per image:
       * - Multiply with CTF
       * - FT backwards
       * - Multiply with the other factors
       */
      Try { 
        gfunc3_mul_vfunc (gf_tmp, &ctf); 
        fft_backward (gf_tmp);
        gfunc3_scale (gf_tmp, 2 * M_PI * M_SQRT2PI);
      }  CATCH_RETURN_VOID (_e);

      /* Return to reciprocal grid */
      Try { gfunc3_grid_fwd_reciprocal (gf_tmp); }  CATCH_RETURN_VOID (_e);  
    }

  gf_tmp->fvals = NULL;
  gfunc3_free (&gf_tmp);

  free (freqs);
  return;
}

/*-------------------------------------------------------------------------------------------------*/
