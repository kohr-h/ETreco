/*
 * et_params.c -- functions to handle ET input parameter structures
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
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "CException.h"

#include "dictionary.h"
#include "iniparser.h"

#include "vec3.h"
#include "misc.h"

#include "gfunc3.h"

#include "et_params.h"
#include "et_vfuncs_private.h"

#include "ai_opts.h"
#include "ai_params.h"
#include "ai_vfuncs.h"


/*-------------------------------------------------------------------------------------------------*/

#define ONE_MICROMETER  1E3   // [nm]

mollifier_ft_function moll_types[] = {NULL, ft_delta, ft_gaussian, NULL};

/*-------------------------------------------------------------------------------------------------*/

AiParams *
new_AiParams (void)
{
  CEXCEPTION_T e = EXC_NONE;
  AiParams *params = NULL;
  
  Try { params = (AiParams *) ali16_malloc (sizeof (AiParams)); }  CATCH_RETURN (e, NULL);
    
  params->tilt_axis               = 0;
  params->tilt_axis_rotation      = 0.0;
  params->tilt_axis_par_shift_px  = 0.0;
  params->ctf_trunc               = 0.0;
  params->moll_ft                 = NULL;

  return params;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_free (AiParams **pparams)
{
  if (pparams == NULL)
    return;
  
  free (*pparams);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_assign_from_AiOpts (AiParams *params, const AiOpts *opts)
{
  int itmp;
  float ta_sx, ta_sy;
  double dtmp;
  dictionary *dict;

  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (opts);

  if ((dict = iniparser_load (opts->fname_params)) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read parameters from %s.", opts->fname_params);
      return;
    }

  /* Regularization */

  /* TODO: implement anisotropic mollifier (3 gamma's) */
  vec3_set_all (params->gamma, opts->gamma * ONE_MICROMETER);
  params->moll_ft = moll_types[opts->moll_type];


  /* VOLUME */

  if ((itmp = iniparser_getint (dict, "volume:nx", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nx' not found in %s.", opts->fname_params);
      return;
    }

  params->vol_shape[0] = itmp;
  
  if ((itmp = iniparser_getint (dict, "volume:ny", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'ny' not found in %s.", opts->fname_params);
      return;
    }

  params->vol_shape[1] = itmp;

  if ((itmp = iniparser_getint (dict, "volume:nz", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nz' not found in %s.", opts->fname_params);
      return;
    }

  params->vol_shape[2] = itmp;

  if ((dtmp = iniparser_getdouble (dict, "volume:voxel_size", FLT_MAX)) == FLT_MAX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'voxel_size' not found in %s.", opts->fname_params);
      return;
    }

  vec3_set_all (params->vol_csize, (float) dtmp);

  /* Overrides MRC header */
  params->vol_shift_px[0] = (float) iniparser_getdouble (dict, "volume:shift_x", FLT_MAX);
  params->vol_shift_px[1] = (float) iniparser_getdouble (dict, "volume:shift_y", FLT_MAX);
  params->vol_shift_px[2] = (float) iniparser_getdouble (dict, "volume:shift_z", FLT_MAX);

  
  /* GEOMETRY */
  
  /* Single axis: parallel tilt axis shift */
  dtmp = iniparser_getdouble (dict, "geometry:tilt_axis", 0.0);
  if (fabsf (dtmp) < 45.0)  /* Use "x" backprojection variant */
    {
      params->tilt_axis = 0;
      params->tilt_axis_rotation = (float) dtmp;
    }
  else
    {
      params->tilt_axis = 1;
      /* What's missing to +- 90 degrees */
      params->tilt_axis_rotation = (dtmp > 0) ? (float)(dtmp - 90.0) : (float)(-dtmp + 90.0);
    }

  ta_sx = (float) iniparser_getdouble (dict, "geometry:axis_shift_x", 0.0);
  ta_sy = (float) iniparser_getdouble (dict, "geometry:axis_shift_y", 0.0);

  params->tilt_axis_par_shift_px = ta_sx * sinf (params->tilt_axis_rotation * ONE_DEGREE) 
    + ta_sy * cosf (params->tilt_axis_rotation * ONE_DEGREE);


  
  /* DETECTOR */
  
  /* Overrides MRC header */
  dtmp = iniparser_getdouble (dict, "detector:pixel_size", 0.0);
  if (dtmp != 0.0)
    {
      params->detector_px_size[0] = (float) dtmp * ONE_MICROMETER;
      params->detector_px_size[1] = (float) dtmp * ONE_MICROMETER;
      params->detector_px_size[2] = 1.0;
    }

  iniparser_freedict (dict);


  /* REGULARIZATION 2 */

  if (use_ctf_flag && truncate_ctf_flag)
    params->ctf_trunc = opts->ctf_trunc;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

#define MAX_SECANT_STEPS    15

float
approx_ctf_xval (float yval, EtParams const *params, float xstart0, float xstart1, float *pdf)
{
  int n;
  float xn_1 = xstart0, xn = xstart1, fval, dfval, update;

  for (n = 0; n < MAX_SECANT_STEPS; n++)
    {
      fval = ctf_acr_unscaled_radial (xn, params) - yval;
      dfval = (ctf_acr_unscaled_radial (xn, params) - 
        ctf_acr_unscaled_radial (xn_1, params)) / (xn - xn_1);

      update = -fval / dfval;

      xn_1 = xn;
      xn += update;

      if (fabsf (update / (xstart1 - xstart0)) < 1E-3)
        break;
    }

  if (n == MAX_SECANT_STEPS)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_COMPUTE, "Error finding cutoff value: no convergence after"
        " %d secant steps (last update = %e)!\n", n, update);
      return FLT_MAX;
    }
    
  if (pdf != NULL)
    *pdf = dfval;

  return xn;
}

/*-------------------------------------------------------------------------------------------------*/

void
compute_xover_spline_params (AiParams *ai_params, EtParams const *et_params)
{
  /* This function computes the coefficients of the interpolating splines between the points of the 
   * function graph of 1/CTF where the absolute value of this function reaches the value 
   * ai_params->ctf_trunc. This way, instead of truncating, a smooth transition between the function 
   * values is guaranteed.
   * A secant algorithm is used here to approximate these values. Starting from the zero of the 
   * CTF (can be computed analytically using a and b below), the x value is approximated where 
   * the CTF absolute value is 1/(ai_params->ctf_trunc). 
   * This procedure is repeated for NUM_CTF_LOBES zeros of the CTF. The sign change depends on the
   * parity of the zero.
   */

  CEXCEPTION_T _e = EXC_NONE;

  int i;
  float a, b, c;

  float x0 = 0.0, x1, dx, x0start, x1start, df0, df1, f1, f0, yval = 1.0 / ai_params->ctf_trunc;

  float u0, v0, w0;

  a = et_params->cs / (4 * et_params->wave_number * et_params->wave_number * et_params->wave_number);
  b = et_params->defocus_nominal / (2 * et_params->wave_number);

  for (i = 0; i < NUM_CTF_LOBES; i++)
    {
      c = atan (et_params->acr) - (i + 1) * M_PI;
      /* (i+1)'th zero = x0start = starting point of secant iteration */
      x0start = (b - sqrtf (b * b + 4 * a * c)) / (2 * a);
      /* x1start = second starting point, a little bit into the correct direction */
      x1start = (i % 2) ? 1.005 * x0start : 0.995 * x0start;  
      /* sign of target value depends on the parity of the zero */
      yval = (i % 2) ? -1.0 / ai_params->ctf_trunc : 1.0 / ai_params->ctf_trunc;  

      /* Compute x0 with CTF(x0) = yval; df0 ~ CTF'(x0) */
      Try { x0 = approx_ctf_xval (yval, et_params, x0start, x1start, &df0); } CATCH_RETURN_VOID (_e);
      /* Make interval [x0,x1] symmetric around x0start;
       * The last spline decays to zero with zero derivative.
       */
      x1 = (i < NUM_CTF_LOBES - 1) ? 2 * x0start - x0 : x0start;  

      f0 = 1.0 / yval;
      f1 = (i < NUM_CTF_LOBES - 1) ? 1.0 / ctf_acr_unscaled_radial (x1, et_params) : 0.0;

      df0 = -df0 / (yval * yval); // (f^-1)'(x) = -f'(x)/f(x)^2
      df1 = (i < NUM_CTF_LOBES - 1) ?
        (1.0 / ctf_acr_unscaled_radial (        x1, et_params) - 
         1.0 / ctf_acr_unscaled_radial (0.999 * x1, et_params)) / (0.001 * x1) 
        : 0.0;

      /* Compute coefficients of the spline
       * p(x) = a0 + a1 * (x - x0) + a2 * (x - x0)^2 + a3 * (x - x0)^3, x0 <= x <= x1 
       */
      dx = x1 - x0;
      u0 = (f1 - f0) / dx;
      v0 = (df1 - df0) / dx;
      w0 = (u0 - df0) / dx;
      ai_params->xover_cspline_coeff[i][3] = (v0 - 2 * w0) / dx;  /* a3 */
      ai_params->xover_cspline_coeff[i][2] = 3 * w0 - v0;         /* a2 */
      ai_params->xover_cspline_coeff[i][1] = df0;                 /* a1 */
      ai_params->xover_cspline_coeff[i][0] = f0;                  /* a0 */
      ai_params->xover_pts[i][0] = x0;                            /* x0 */
      ai_params->xover_pts[i][1] = x1;                            /* x1 */
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_assign_ctftrunc_from_EtParams (AiParams *ai_params, EtParams const *et_params)
{
  CEXCEPTION_T _e = EXC_NONE;

  CAPTURE_NULL_VOID (ai_params);
  CAPTURE_NULL_VOID (et_params);
  
  if (ai_params->ctf_trunc > 0.0)
    {
      Try { compute_xover_spline_params (ai_params, et_params); } CATCH_RETURN_VOID (_e);
    }
  return;
}

/*-------------------------------------------------------------------------------------------------*/


void
AiParams_print (AiParams const *params)
{
  int i;
  
  CAPTURE_NULL_VOID (params);

  printf ("\n");
  puts ("AI parameters:");
  puts ("==============\n");

  puts ("Geometry:");
  puts ("---------\n");
  
  printf ("volume voxel size    : % 9.2f [nm]\n", params->vol_csize[0]);
  printf ("volume shift         : (");
  for (i = 0; i < 3; i++)
    {
      if (params->vol_shift_px[i] == FLT_MAX)  printf ("(from data)");
      else  printf ("%7.2f", params->vol_shift_px[i]);
      if (i != 2)  printf (", ");
    }
  printf (") [pixels]\n\n");

  printf ("detector shape      : ");
  if (params->detector_px_size[0] != 0.0)  printf ("% 9.2f [nm]\n", params->detector_px_size[0]);
  else  printf ("(from data)\n");
  printf ("tilt axis rotation  : % 9.2f [degrees]\n", params->tilt_axis_rotation);
  printf ("\n");

  puts ("Regularization:");
  puts ("---------------\n");
  
  printf ("gamma (scaled)       : (%.2e, %.2e, %.2e)\n", params->gamma[0], params->gamma[1], 
    params->gamma[2]);
    
  if (use_ctf_flag && truncate_ctf_flag)
    {
      printf ("recip CTF truncation : %9.2f\n", params->ctf_trunc);
      printf ("\n");
      if (verbosity_level >= VERB_LEVEL_VERBOSE)
        {
          printf ("truncation intervals and spline parameters:\n");

          for (i = 0; i < NUM_CTF_LOBES; i++)
            {
              printf ("[% 9.6f, % 9.6f] a0 = % f, a1 = % f, a2 = % f, a3 = % f\n", 
                params->xover_pts[i][0], params->xover_pts[i][1], 
                params->xover_cspline_coeff[i][0], params->xover_cspline_coeff[i][1], 
                params->xover_cspline_coeff[i][2], params->xover_cspline_coeff[i][3]);
            }
        }
    }
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_apply_to_volume (AiParams const *params, gfunc3 *vol)
{
  int i;
  
  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (vol);
  
  for (i = 0; i < 3; i++)
    {
      if (params->vol_shift_px[i] != FLT_MAX)
        vol->x0[i] = params->vol_shift_px[i] * vol->csize[i];
    }
  
  gfunc3_compute_xmin_xmax (vol);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_apply_to_proj_image (AiParams const *params, gfunc3 *proj_img)
{
  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (proj_img);
  
  if (!GFUNC_IS_2D (proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "proj_img must be 2-dimensional");
      return;
    }
  
  if (params->detector_px_size[0] != 0.0)  /* Detector pixel size is set in config file */
    gfunc3_set_csize (proj_img, params->detector_px_size);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/
