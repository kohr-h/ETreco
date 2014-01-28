/*
 * et_vfuncs.c -- vector functions related to ET
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
envelope_energy_spread (float t, EtParams const *params)
{
  float ct = params->energy_spread * params->cc1 * t / (4 * params->wave_number);

  return exp (-ct * ct);
}

/*-------------------------------------------------------------------------------------------------*/

float
envelope_source_size (float t, EtParams const *params)
{
  float c = params->wave_number * params->cond_ap_angle / (2 * params->focal_length);
  float b = params->cs * t / (params->focal_length * params->focal_length) - params->defocus_nominal;

  return exp (-c * c * b * b);
}

/*-------------------------------------------------------------------------------------------------*/

float
oscillating_part_angle_function (float t, EtParams const *params)
{
  return -t / (4 * params->wave_number) 
    * (params->cs * t / (params->wave_number * params->wave_number) - 2 * params->defocus_nominal);
}

/*-------------------------------------------------------------------------------------------------*/

float
oscillating_part_angle_function_acr (float t, EtParams const *params)
{
  float zp = -t / (4 * params->wave_number) 
    * (params->cs * t / (params->wave_number * params->wave_number) - 2 * params->defocus_nominal);

  return zp + atan (params->acr);
}

/*-------------------------------------------------------------------------------------------------*/

float
ctf_scaling_function (float t, EtParams const *params)
{
  /* Account for volume dilation factor here: sqrt(m) -> m^2 */
  float c = params->magnification * params->magnification * params->wave_number * params->wave_number;

  double a = sqrt ( (double) params->wave_number * params->wave_number - t);

  return (float) (c / a);
}

/*-------------------------------------------------------------------------------------------------*/

float
ctf_acr_scaling_function (float t, EtParams const *params)
{
  /* Account for volume dilation factor here: sqrt(m) -> m^2 */
  float c = params->magnification * params->magnification * params->wave_number * params->wave_number
    * sqrtf (1 + params->acr * params->acr);

  double a = sqrt ( (double) params->wave_number * params->wave_number - t);

  return (float) (c / a);
}

/*-------------------------------------------------------------------------------------------------*/

float complex
ctf_unscaled_radial (float t, EtParams const *params)
{
  float r = params->aper_cutoff;

  if (t > r * r)
    return 0.0;

  return cexpf (I * oscillating_part_angle_function (t, params)) 
    * envelope_source_size (t, params) * envelope_energy_spread (t, params);
}

/*-------------------------------------------------------------------------------------------------*/

float
ctf_acr_unscaled_radial (float t, EtParams const *params)
{
  float r = params->aper_cutoff;

  if (t > r * r)
    return 0.0;

  return sinf (oscillating_part_angle_function_acr (t, params)) 
    * envelope_source_size (t, params) * envelope_energy_spread (t, params);
}

/*-------------------------------------------------------------------------------------------------*/

void
ctf_acr (float const *xi, float *zp, void const *par)
{
  EtParams *params = (EtParams *) par;
  float mxi2;

  mxi2 = xi[0] * xi[0] + xi[1] * xi[1];
  mxi2 *= params->magnification * params->magnification;

  *zp = ctf_acr_scaling_function (mxi2, params) * ctf_acr_unscaled_radial (mxi2, params) / (2 * M_PI);

  return;
}

void
vfunc_init_ctf_acr (vfunc *vf, EtParams const *params)
{
  CAPTURE_NULL_VOID (vf);
  CAPTURE_NULL_VOID (params);
    
  vf->f = ctf_acr;
  vf->params = params;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
ctf (float const *xi, float *zp, void const *par)
{
  EtParams *params = (EtParams *) par;
  float mxi2;
  float complex *z = (float complex *) zp;

  mxi2 = xi[0] * xi[0] + xi[1] * xi[1];
  mxi2 *= params->magnification * params->magnification;

  *z = ctf_scaling_function (mxi2, params) * ctf_unscaled_radial (mxi2, params) / (2 * M_PI);

  return;
}

void
vfunc_init_ctf (vfunc *vf, EtParams const *params)
{
  CAPTURE_NULL_VOID (vf);
  CAPTURE_NULL_VOID (params);
    
  vf->f = ctf;
  vf->params = params;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

float
detector_mtf_radial (float t, EtParams const *params)
{
  float fa = 1.0 + params->mtf_alpha * t;
  float fb = 1.0 + params->mtf_beta * t;
  float pp = powf (2, params->mtf_p) / (2 * powf (fa, params->mtf_p) + powf (2, params->mtf_p) - 2.0);
  float pq = powf (2, params->mtf_q) / (2 * powf (fb, params->mtf_q) + powf (2, params->mtf_q) - 2.0);

  return params->mtf_a * pp + params->mtf_b * pq + params->mtf_c;
}

/*-------------------------------------------------------------------------------------------------*/

void
detector_mtf (float const *grid_xi, float *zp, void const *par)
{
  EtParams *params = (EtParams *) par;
  float absxi2 = grid_xi[0] * grid_xi[0] + grid_xi[1] * grid_xi[1];

  *zp = detector_mtf_radial (absxi2, params);

  return;
}

void
vfunc_init_detector_mtf (vfunc *vf, EtParams const *params)
{
  CAPTURE_NULL_VOID (vf);
  CAPTURE_NULL_VOID (params);
  
  vf->f = detector_mtf;
  vf->params = params;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
detector_recip_mtf (float const *grid_xi, float *zp, void const *par)
{
  EtParams *params = (EtParams *) par;
  float absxi2 = grid_xi[0] * grid_xi[0] + grid_xi[1] * grid_xi[1];

  *zp = 1.0 / detector_mtf_radial (absxi2, params);

  return;
}

void
vfunc_init_detector_recip_mtf (vfunc *vf, EtParams const *params)
{
  CAPTURE_NULL_VOID (vf);
  CAPTURE_NULL_VOID (params);
  
  vf->f = detector_recip_mtf;
  vf->params = params;

  return;
}

/*-------------------------------------------------------------------------------------------------*/
