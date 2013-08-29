/*
 * et_vfuncs.c -- vector functions related to ET
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#include "CException.h"

#include "misc.h"
#include "simd.h"

#include "vfunc.h"

#include "et_params.h"
#include "et_vfuncs.h"

/*-------------------------------------------------------------------------------------------------*/

float
ft_delta (float const *xi, vec3 const gamma)
{
  return 1.0 / (2 * M_PI * M_SQRT2PI);
}

/*-------------------------------------------------------------------------------------------------*/

float
ft_gaussian (float const *xi, vec3 const gamma)
{
  /* width of gaussian with gamma^2=1 is 3, thus compensate by 6.0 instead of 2.0 in denominator */
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) xi, *p2 = (__m128 const *) gamma; 
  __m128 mgxi2 = _mm_mul_ps (*p1, *p2);
  mgxi2 = _mm_mul_ps (mgxi2, mgxi2);
  float *pmgxi2 = (float *) &mgxi2;
  return expf (-(pmgxi2[0] + pmgxi2[1] + pmgxi2[2]) / 6.0)  / (2 * M_PI * M_SQRT2PI);

  #else
  float gxi0 = xi[0] * gamma[0], gxi1 = xi[1] * gamma[1], gxi2 = xi[2] * gamma[2];
  return expf (-(gxi0 * gxi0 + gxi1 * gxi1 + gxi2 * gxi2) / 6.0) / (2 * M_PI * M_SQRT2PI);

  #endif
}

/*-------------------------------------------------------------------------------------------------*/

float
envelope_energy_spread (float t, RecParams const *rec_p)
{
  float ct = rec_p->energy_spread * rec_p->cc1 * t / (4 * rec_p->wave_number);

  return exp (-ct * ct);
}

/*-------------------------------------------------------------------------------------------------*/

float
envelope_source_size (float t, RecParams const *rec_p)
{
  float c = rec_p->wave_number * rec_p->cond_ap_angle / (2 * rec_p->focal_length);
  float b = rec_p->cs * t / (rec_p->focal_length * rec_p->focal_length) - rec_p->defocus_nominal;

  return exp (-c * c * b * b);
}

/*-------------------------------------------------------------------------------------------------*/

float
oscillating_part_angle_function (float t, RecParams const *rec_p)
{
//   float ls = rec_p->magnification / (rec_p->magnification - 1) * rec_p->focal_length;

//  // The following term is probably close to zero (CAUTION: possible OVERFLOW)
//   float b = - ls * ( t*t*t*t / (8*rec_p->wave_number*rec_p->wave_number*rec_p->wave_number) );

  float zp = -t / (4 * rec_p->wave_number) 
    * (rec_p->cs * t / (rec_p->wave_number * rec_p->wave_number) - 2 * rec_p->defocus_nominal);

  return /*b + */ zp + atan (rec_p->acr);
}

/*-------------------------------------------------------------------------------------------------*/

float
ctf_scaling_function (float t, RecParams const *rec_p)
{
  // float c = rec_p->magnification * rec_p->magnification * rec_p->wave_number * rec_p->wave_number
    // * sqrtf (1 + rec_p->acr * rec_p->acr);
  float c = rec_p->wave_number;

  float a = sqrtf (rec_p->wave_number * rec_p->wave_number - t);

  return c / a;
}

/*-------------------------------------------------------------------------------------------------*/

float
ctf_unscaled_radial (float t, RecParams const *rec_p)
{
  float r = rec_p->aper_cutoff;

  if (t > r * r)
    return 0.0;

  return sinf (oscillating_part_angle_function (t, rec_p)) 
    * envelope_source_size (r, rec_p) * envelope_energy_spread (t, rec_p);
}

/*-------------------------------------------------------------------------------------------------*/

float
detector_mtf_radial (float t, RecParams const *rec_p)
{
  float fa = 1.0 + rec_p->mtf_alpha * t;
  float fb = 1.0 + rec_p->mtf_beta * t;
  float pp = powf (2, rec_p->mtf_p) / (2 * powf (fa, rec_p->mtf_p) + powf (2, rec_p->mtf_p) - 2.0);
  float pq = powf (2, rec_p->mtf_q) / (2 * powf (fb, rec_p->mtf_q) + powf (2, rec_p->mtf_q) - 2.0);

  return rec_p->mtf_a * pp + rec_p->mtf_b * pq + rec_p->mtf_c;
}

/*-------------------------------------------------------------------------------------------------*/

void
detector_mtf (float const *grid_xi, float *zp, void const *params)
{
  RecParams *rec_p = (RecParams *) params;
  float absxi2 = grid_xi[0] * grid_xi[0] + grid_xi[1] * grid_xi[1];

  *zp = detector_mtf_radial (absxi2, rec_p);

  return;
}

void
vfunc_init_detector_mtf (vfunc *vf, RecParams const *rec_p)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (rec_p);
  
  vf->f = detector_mtf;
  vf->params = rec_p;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
detector_recip_mtf (float const *grid_xi, float *zp, void const *params)
{
  RecParams *rec_p = (RecParams *) params;
  float absxi2 = grid_xi[0] * grid_xi[0] + grid_xi[1] * grid_xi[1];

  *zp = 1.0 / detector_mtf_radial (absxi2, rec_p);

  return;
}

void
vfunc_init_detector_recip_mtf (vfunc *vf, RecParams const *rec_p)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (rec_p);
  
  vf->f = detector_recip_mtf;
  vf->params = rec_p;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ctf (float const *xi, float *zp, void const *params)
{
  RecParams *rec_p = (RecParams *) params;
  float mxi2;

  mxi2 = xi[0] * xi[0] + xi[1] * xi[1];
  mxi2 *= rec_p->magnification * rec_p->magnification;

  *zp = ctf_scaling_function (mxi2, rec_p) * ctf_unscaled_radial (mxi2, rec_p);

  return;
}

void
vfunc_init_ctf (vfunc *vf, RecParams const *rec_p)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (rec_p);
    
  vf->f = ctf;
  vf->params = rec_p;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

float
ramp_filter_sax (float const *xi)
{
  return fabsf (xi[1]) / M_SQRT2PI;
}

/*-------------------------------------------------------------------------------------------------*/

float
recip_ctf_unscaled_radial (float t, RecParams const *rec_p)
{
  int i;
  float r = rec_p->aper_cutoff;
  float t0, t1, dt, a0, a1, a2, a3;

  t1 = rec_p->xover_pts[NUM_CTF_LOBES - 1][1]; // 0 from last zero on
  if ((t > t1) || (t > r * r))
    return 0.0;

  for (i = 0; i < NUM_CTF_LOBES; i++)
    {
      t0 = rec_p->xover_pts[i][0];
      t1 = rec_p->xover_pts[i][1];
      if ((t > t0) && (t < t1))
        {
          a0 = rec_p->xover_cspline_coeff[i][0];
          a1 = rec_p->xover_cspline_coeff[i][1];
          a2 = rec_p->xover_cspline_coeff[i][2];
          a3 = rec_p->xover_cspline_coeff[i][3];
          dt = t - t0;

          return a0 + dt * (a1 + dt * (a2 + dt * a3));
        }
    }

  return 1.0 / ctf_unscaled_radial (t, rec_p);
}

/*-------------------------------------------------------------------------------------------------*/

/* TODO: insert the correct constants */

void ft_rk_single_axis_x (float const *xi, float *zp, void const *params)
{
  RecParams *rec_p = (RecParams *) params;
  float mxi2;
  
  mxi2 = xi[0] * xi[0] + xi[1] * xi[1];
  mxi2 *= rec_p->magnification * rec_p->magnification;

  if (sqrtf (mxi2) >= rec_p->aper_cutoff)
    *zp = 0.0;

  else
    {
      *zp = recip_ctf_unscaled_radial (mxi2, rec_p) / ctf_scaling_function (mxi2, rec_p) 
        * ramp_filter_sax (xi) * rec_p->moll_ft (xi, rec_p->gamma);
    }
        
  return;
}

void ft_rk_single_axis_x_noctf (float const *xi, float *zp, void const *params)
{
  RecParams *rec_p = (RecParams *)params;

  *zp = ramp_filter_sax (xi) * rec_p->moll_ft (xi, rec_p->gamma);
  
  return;
}

void
vfunc_init_ft_rk_single_axis_x (vfunc *vf, RecParams const *rec_p)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (rec_p);
    
  if (use_ctf_flag)
    vf->f = ft_rk_single_axis_x;
  else
    vf->f = ft_rk_single_axis_x_noctf;
    
  vf->params = rec_p;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ft_charfun_ball3 (float const *xi, float *zp, void const *params)
{
  float rad = *((float *) params);
  float tmp, r_absxi = rad * sqrtf (xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);

  if (r_absxi < 0.01)
    {
      /* asymptotic form for small argument: 
       * 1/2^{3/2}*(1/Gamma(5/2) - 1/Gamma(7/2) * z/2)
       */
      tmp = 2.0 / (3 * M_SQRT2PI) * (1 - r_absxi / 5.0);
    }
  else if (r_absxi > 10.0)
    {
      /* asymptotic form for large argument: 
       * sqrt(2/(pi*z))/z^{3/2}*(cos(z-pi) - sin(z-pi)/z)
       */
      tmp = sqrtf (2 / M_PI) / (r_absxi * r_absxi) 
        * (cosf (r_absxi - M_PI) - sinf (r_absxi - M_PI) / r_absxi);
    }
  else
    {
      /* regular form: J_(3/2)(z) / z^(3/2) */
      tmp = gsl_sf_bessel_Jnu (1.5, r_absxi) / (r_absxi * sqrtf (r_absxi));
    }

  *zp = tmp * rad * rad * rad;

  return;
}

void
vfunc_init_ft_charfun_ball3 (vfunc *vf, float const *pradius)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (pradius);
  
  vf->f = ft_charfun_ball3;
  vf->params = pradius;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ft_charfun_cyl3 (float const *xi, float *zp, void const *params)
{
  // TODO: check and explain
  float *p = (float *) params;
  float len = p[0], rad = p[1];

  float argx, argyz, sinc, bessel;

  argx = len * xi[0];
  argyz = rad * sqrtf (xi[1] * xi[1] + xi[2] * xi[2]);

  sinc = (fabsf (argx) > 1E-4) ? sinf (argx) / argx : 1.0 - argx * argx / 3.0;

  if (fabsf (argyz) < 1E-2)
    bessel = 0.5 * (1 - argyz / 4.0);
  else if (fabsf (argyz) > 10.0)
    bessel = sqrtf (2.0 / (M_PI * argyz)) / argyz
      * (cosf (argyz - 3 * M_PI / 4.0) - sinf (argyz - 3 * M_PI / 4.0) * 3.0 / (8.0 * argyz));
  else
    bessel = gsl_sf_bessel_J1 (argyz) / argyz;

  *zp = len * sinc * rad * rad * bessel / M_SQRT2PI;

  return;
}

void
vfunc_init_ft_charfun_cyl3 (vfunc *vf, float const *plength_radius)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (plength_radius);
    
  vf->f = ft_charfun_cyl3;
  vf->params = plength_radius;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ft_lambda (float const *xi, float *zp, void const *params)
{
  float a = *((float const *) params);
  float absxi = sqrtf (xi[0] * xi[0] + xi[1] * xi[1]);

  if (a >= 0.0)
    *zp = powf (absxi, a);
  else
  *zp = (absxi > EPS_DENOM) ? powf (absxi, a) : powf (EPS_DENOM, a);

  return;
}

void
vfunc_init_ft_lambda (vfunc *vf, float const *ppow)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (ppow);
    
  vf->f = ft_lambda;
  vf->params = ppow;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ft_lambda_v (float const *xi, float *zp, void const *params)
{
  float a = *((float *) params);
  float absxi = fabs (xi[1]);

  if (a >= 0.0)
    *zp = powf (absxi, a);
  else
  *zp = (absxi > EPS_DENOM) ? powf (absxi, a) : powf (EPS_DENOM, a);

  return;
}

void
vfunc_init_ft_lambda_v (vfunc *vf, float const *ppow)
{
  CAPTURE_NULL (vf);
  CAPTURE_NULL (ppow);
    
  vf->f = ft_lambda_v;
  vf->params = ppow;

  return;
}

/*-------------------------------------------------------------------------------------------------*/
