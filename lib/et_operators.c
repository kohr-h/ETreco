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

#include "fft.h"
#include "gfunc3.h"
#include "mrc.h"
#include "operators.h"
#include "operators_private.h"

#include "et_params.h"

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
    freqs = (float *) ali16_malloc ((ft_proj_img_grid->ntotal * 3) * sizeof (float));
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

void
et_scattering_projection (gfunc3 const *scatterer, vec3 const angles_deg, RecParams const *rec_p, 
                          gfunc3 *proj_img, scattering_model sct_model)
{
  CEXCEPTION_T _e = EXC_NONE;

  int real_output = FALSE;
  float *freqs = NULL;
  vfunc ctf;
  
  CAPTURE_NULL_VOID (scatterer);
  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (angles_deg);
  if (sct_model != PROJ_ASSUMPTION)
    CAPTURE_NULL_VOID (rec_p);
  
  GFUNC_CAPTURE_UNINIT_VOID (scatterer);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);

  if (proj_img->type == REAL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Projection image must be of COMPLEX type.");
      return;
    }

  if ( (sct_model == BORN_APPROX) && (rec_p->wave_number == 0.0) )
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
        nfft3_transform (scatterer, freqs, proj_img->ntotal, (float complex *) proj_img->fvals); 
      } CATCH_RETURN_VOID (_e);
   }
  else if (sct_model == BORN_APPROX)
    {
      Try { 
        freqs = ewald_sphere_freqs (proj_img, angles_deg, rec_p->wave_number); 
      } CATCH_RETURN_VOID (_e);

      Try { 
        nfft3_transform (scatterer, freqs, proj_img->ntotal, (float complex *) proj_img->fvals); 
      } CATCH_RETURN_VOID (_e);
    }
  
  /* Multiply with CTF and scaling factors. Finally apply inverse FT */
  Try { vfunc_init_ctf (ctf, rec_p); } CATCH_RETURN_VOID (_e);
  Try {
    gfunc3_mul_vfunc (proj_img, ctf);
    gfunc3_scale (proj_img, 2 * M_PI * M_SQRT2PI);
  } CATCH_RETURN_VOID (_e);

  Try { fft_backward (proj_img); } CATCH_RETURN_VOID (_e);
  
  free (freqs);
  return;
}

/*-------------------------------------------------------------------------------------------------*/
