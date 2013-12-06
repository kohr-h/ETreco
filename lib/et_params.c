/*
 * et_params.c -- functions to handle ET input parameter structures
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
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "CException.h"

#include "dictionary.h"
#include "iniparser.h"

#include "vec3.h"
#include "misc.h"

#include "gfunc3.h"

#include "ai_options.h"
#include "et_params.h"
#include "et_vfuncs.h"
#include "et_vfuncs_private.h"

// Physical constants and conversion factors
#define EL_REST_ENERGY  510998.928  // [eV]
#define HC              1239.84193  // [eV*nm]

#define ONE_MILLIMETER  1E6   // [nm]
#define ONE_MICROMETER  1E3   // [nm]
#define ONE_KILOVOLT    1E3   // [V]
#define ONE_MILLIRADIAN 1E-3  // [1]

// Amplitude contrast ratio, fixed value seems to be suitable for organic specimens
#define ACR             0.2

// Parameter for secant method
#define MAX_STEPS    15

mollifier_ft_function moll_types[] = {NULL, ft_delta, ft_gaussian, NULL};

/*-------------------------------------------------------------------------------------------------*/

int use_ctf_flag = 0;
int use_mtf_flag = 0;

/*-------------------------------------------------------------------------------------------------*/

RecParams *
new_RecParams (void)
{
  CEXCEPTION_T e = EXC_NONE;
  RecParams *rp = NULL;
  
  Try { rp = (RecParams *) ali16_malloc (sizeof (RecParams)); }  CATCH_RETURN (e, NULL);
    
  rp->acc_voltage             = 0.0;
  rp->energy_spread           = 0.0;
  rp->magnification           = 0.0;
  rp->cs                      = 0.0;
  rp->cc                      = 0.0;
  rp->aperture                = 0.0;
  rp->focal_length            = 0.0;
  rp->cond_ap_angle           = 0.0;
  rp->defocus_nominal         = 0.0;
  rp->mtf_a                   = 0.0;
  rp->mtf_b                   = 0.0;
  rp->mtf_c                   = 0.0;
  rp->mtf_alpha               = 0.0;
  rp->mtf_beta                = 0.0;
  rp->mtf_p                   = 0;
  rp->mtf_q                   = 0;
  rp->acr                     = 0.0;
  rp->tilt_axis               = 0;
  rp->tilt_axis_rotation      = 0.0;
  rp->tilt_axis_par_shift_px  = 0.0;
  rp->wave_number             = 0.0;
  rp->cc1                     = 0.0;
  rp->aper_cutoff             = 0.0;
  rp->ctf_trunc               = 0.0;
  rp->moll_ft                 = NULL;

  return rp;
}

/*-------------------------------------------------------------------------------------------------*/

void
RecParams_free (RecParams **prec_p)
{
  if (prec_p == NULL)
    return;
  
  free (*prec_p);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

/* Internal declaration */
void
compute_xover_spline_params (RecParams *rec_p);


void
RecParams_assign_from_OptionData (RecParams *rec_p, const OptionData *od)
{
  int itmp;
  float ta_sx, ta_sy;
  double dtmp;
  dictionary *dict;

  CAPTURE_NULL_VOID (rec_p);
  CAPTURE_NULL_VOID (od);

  if ((dict = iniparser_load (od->fname_reco_params)) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read reco parameters from %s.", 
        od->fname_reco_params);
      return;
    }

  /* Regularization */

  /* TODO: implement anisotropic mollifier (3 gamma's) */
  vec3_set_all (rec_p->gamma, od->gamma * ONE_MICROMETER);
  rec_p->moll_ft = moll_types[od->moll_type];


  /* Geometry part */

  if ((itmp = rec_p->mtf_p = iniparser_getint (dict, "volume:nx", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nx' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->vol_shape[0] = itmp;
  
  if ((itmp = rec_p->mtf_p = iniparser_getint (dict, "volume:ny", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'ny' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->vol_shape[1] = itmp;

  if ((itmp = rec_p->mtf_p = iniparser_getint (dict, "volume:nz", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nz' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->vol_shape[2] = itmp;

  if ((dtmp = iniparser_getdouble (dict, "volume:voxel_size", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'voxel_size' not found in %s.", od->fname_reco_params);
      return;
    }

  vec3_set_all (rec_p->vol_csize, (float) dtmp);

  /* Overrides MRC header */
  rec_p->vol_shift_px[0] = (float) iniparser_getdouble (dict, "volume:shift_x", FLT_MAX);
  rec_p->vol_shift_px[1] = (float) iniparser_getdouble (dict, "volume:shift_y", FLT_MAX);
  rec_p->vol_shift_px[2] = (float) iniparser_getdouble (dict, "volume:shift_z", FLT_MAX);

  
  /* Single axis: parallel tilt axis shift */
  dtmp = iniparser_getdouble (dict, "geometry:tilt_axis", 0.0);
  if (fabsf (dtmp) < 45.0)  /* Use "x" backprojection variant */
    {
      rec_p->tilt_axis = 0;
      rec_p->tilt_axis_rotation = (float) dtmp;
    }
  else
    {
      rec_p->tilt_axis = 1;
      /* What's missing to +- 90 degrees */
      rec_p->tilt_axis_rotation = (dtmp > 0) ? (float)(dtmp - 90.0) : (float)(-dtmp + 90.0);
    }

  ta_sx = (float) iniparser_getdouble (dict, "geometry:axis_shift_x", 0.0);
  ta_sy = (float) iniparser_getdouble (dict, "geometry:axis_shift_y", 0.0);

  rec_p->tilt_axis_par_shift_px = ta_sx * sinf (rec_p->tilt_axis_rotation * ONE_DEGREE) 
    + ta_sy * cosf (rec_p->tilt_axis_rotation * ONE_DEGREE);


  if ((dtmp = iniparser_getdouble (dict, "optics:magnification", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'magnification' not found in %s.", 
        od->fname_reco_params);
      return;
    }

  rec_p->magnification = (float) dtmp;
  
  /* Overrides MRC header */
  dtmp = iniparser_getdouble (dict, "detector:pixel_size", 0.0);
  rec_p->detector_px_size[0] = (float) dtmp * ONE_MICROMETER;
  rec_p->detector_px_size[1] = (float) dtmp * ONE_MICROMETER;
  rec_p->detector_px_size[2] = 1.0;


  /* CTF part */

  /* Ignore CTF if acc_voltage is zero or not present */
  if ((dtmp = iniparser_getdouble (dict, "electronbeam:acc_voltage", 0.0)) == 0.0)
    {
      use_ctf_flag = 0;
      iniparser_freedict (dict);
      return;
    }

  use_ctf_flag = 1;
  rec_p->acc_voltage = (float) dtmp * ONE_KILOVOLT;

  if ((dtmp = iniparser_getdouble (dict, "electronbeam:energy_spread", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'energy_spread' not found in %s.", 
        od->fname_reco_params);
      return;
    }

  rec_p->energy_spread = (float) dtmp;

  if ((dtmp = iniparser_getdouble (dict, "optics:cs", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cs' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->cs = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:cc", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cc' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->cc = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:aperture", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'aperture' not found in %s.", od->fname_reco_params);
      return;
    }

  rec_p->aperture = (float) dtmp * ONE_MICROMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:focal_length", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'focal_length' not found in %s.", 
        od->fname_reco_params);
      return;
    }

  rec_p->focal_length = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:cond_ap_angle", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cond_ap_angle' not found in %s.", 
        od->fname_reco_params);
      return;
    }

  rec_p->cond_ap_angle = (float) dtmp * ONE_MILLIRADIAN;

  if ((dtmp = iniparser_getdouble (dict, "optics:defocus_nominal", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'defocus_nominal' not found in %s.", 
        od->fname_reco_params);
      return;
    }

  rec_p->defocus_nominal = (float) dtmp * ONE_MICROMETER;

  /* If not found in config file, set to predefined ACR */
  rec_p->acr = (float) iniparser_getdouble (dict, "sample:famp", ACR);

  /* The MTF parameters all have default values */

  use_mtf_flag = 1;
  
  rec_p->mtf_a     = iniparser_getdouble (dict, "detector:mtf_a", 0.0);
  rec_p->mtf_b     = iniparser_getdouble (dict, "detector:mtf_b", 0.0);
  rec_p->mtf_c     = iniparser_getdouble (dict, "detector:mtf_c", 1.0);
  rec_p->mtf_alpha = iniparser_getdouble (dict, "detector:mtf_alpha", 0.0);
  rec_p->mtf_beta  = iniparser_getdouble (dict, "detector:mtf_beta", 0.0);
  rec_p->mtf_p     = iniparser_getint (dict, "detector:mtf_p", 1);  
  rec_p->mtf_q     = iniparser_getint (dict, "detector:mtf_q", 1);  

  /* If a and b are zero, the MTF collapses to a constant. This means 'no MTF'. */
  if ((rec_p->mtf_a == 0.0) && (rec_p->mtf_b == 0.0))
    use_mtf_flag = 0;
   

  iniparser_freedict (dict);


  /* Derived parameters */

  /* Compute relativistic wave number (unit: [1/nm]) */

  /* momentum * c [eV] */
  dtmp = sqrtf (rec_p->acc_voltage * rec_p->acc_voltage + 2 * rec_p->acc_voltage * EL_REST_ENERGY); 
  rec_p->wave_number = (float) 2 * M_PI * dtmp / HC;


  /* Compute constant derived from cc (unit: [nm]) */
  dtmp = 1.0 / (2 * EL_REST_ENERGY); // some factor
  rec_p->cc1 = (float) (1 + 2 * dtmp * rec_p->acc_voltage) 
                / (rec_p->acc_voltage * (1 + dtmp * rec_p->acc_voltage)) * rec_p->cc;

  rec_p->aper_cutoff = (rec_p->wave_number * rec_p->aperture) / rec_p->focal_length;


  /* Regularization part 2 */

  if (use_ctf_flag && truncate_ctf_flag)
    rec_p->ctf_trunc = od->ctf_trunc;
    
  if (rec_p->ctf_trunc > 0.0)
    compute_xover_spline_params (rec_p);

  return;
}

/*-------------------------------------------------------------------------------------------------*/


void
RecParams_print (RecParams const *rec_p)
{
  int i;
  
  CAPTURE_NULL_VOID (rec_p);

  /* TODO: fix alignment of printout */
  /* TODO: make dependent on verbosity */
  printf ("\n");
  puts ("Reconstruction parameters:");
  puts ("==========================\n");

  if (use_ctf_flag)
    {
      puts ("CTF part:");
      puts ("---------\n");
      
      printf ("acc_voltage     : % 7.2e [V]\n", rec_p->acc_voltage);
      printf ("energy_spread   : % 7.2f [eV]\n", rec_p->energy_spread);
      printf ("magnification   : % 7.2f\n", rec_p->magnification);
      printf ("cs              : % 7.2e [nm]\n", rec_p->cs);
      printf ("cc              : % 7.2e [nm]\n", rec_p->cc);
      printf ("aperture        : % 7.2f [nm]\n", rec_p->aperture);
      printf ("focal_length    : % 7.2e [nm]\n", rec_p->focal_length);
      printf ("cond_ap_angle   : % 7.5f [rad]\n", rec_p->cond_ap_angle);
      printf ("defocus_nominal : % 7.2f [nm]\n", rec_p->defocus_nominal);
      printf ("mtf_a           : % 7.2f\n", rec_p->mtf_a);
      printf ("mtf_b           : % 7.2f\n", rec_p->mtf_b);
      printf ("mtf_c           : % 7.2f\n", rec_p->mtf_c);
      printf ("mtf_alpha       : % 7.2f\n", rec_p->mtf_alpha);
      printf ("mtf_beta        : % 7.2f\n", rec_p->mtf_beta);
      printf ("mtf_p           : % 7d\n", rec_p->mtf_p);
      printf ("mtf_q           : % 7d\n", rec_p->mtf_q);
      printf ("\n");
      printf ("acr             : % 7.2f\n", rec_p->acr);
      printf ("wave_number     : % 7.2f [1/nm]\n", rec_p->wave_number);
      printf ("cc1             : % 7.2f [nm]\n", rec_p->cc1);
      printf ("aperture cutoff : % 7.2f [1/nm]\n", rec_p->aper_cutoff);
      printf ("\n");
    }
    
  puts ("Geometry part:");
  puts ("------------\n");
  
  printf ("vol_shape    : (%d, %d, %d)\n", rec_p->vol_shape[0], rec_p->vol_shape[1], 
  rec_p->vol_shape[2]);
  printf ("vol_csize    : % 7.2f [nm]\n", rec_p->vol_csize[0]);
  printf ("vol_shift_px : (");
  for (i = 0; i < 3; i++)
    {
      if (rec_p->vol_shift_px[i] == FLT_MAX)  printf ("(from data)");
      else  printf ("%7.2f", rec_p->vol_shift_px[i]);
      if (i != 2)  printf (", ");
    }
  printf (")\n\n");

  printf ("detector_px_size  :");
  if (rec_p->detector_px_size[0] != 0.0)  printf ("% 7.2f\n", rec_p->detector_px_size[0]);
  else  printf ("(from data)\n");
  printf ("tilt_axis       : % 7.2f [degrees]\n", rec_p->tilt_axis_rotation);
  printf ("\n");
  if (!use_ctf_flag)
    {
      printf ("magnification   : % 7.2f\n", rec_p->magnification);
      printf ("\n");
    }

  puts ("Regularization part:");
  puts ("--------------------\n");
  
  printf ("gamma = (%e, %e, %e)\n", rec_p->gamma[0], rec_p->gamma[1], rec_p->gamma[2]);
  if (use_ctf_flag && truncate_ctf_flag)
    {
      printf ("reciprocal CTF truncation at M = %f\n", rec_p->ctf_trunc);
      printf ("\n");
      printf ("truncation intervals and spline parameters:\n");

      for (i = 0; i < NUM_CTF_LOBES; i++)
        {
          printf ("[% 9.6f, % 9.6f] a0 = % f, a1 = % f, a2 = % f, a3 = % f\n", rec_p->xover_pts[i][0],
            rec_p->xover_pts[i][1], rec_p->xover_cspline_coeff[i][0],
            rec_p->xover_cspline_coeff[i][1], rec_p->xover_cspline_coeff[i][2],
            rec_p->xover_cspline_coeff[i][3]);
        }
    }
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
RecParams_apply_to_volume (RecParams const *rec_p, gfunc3 *vol)
{
  int i;
  
  CAPTURE_NULL_VOID (rec_p);
  CAPTURE_NULL_VOID (vol);
  
  for (i = 0; i < 3; i++)
    {
      if (rec_p->vol_shift_px[i] == FLT_MAX)  continue;
      
      vol->x0[i] = rec_p->vol_shift_px[i] * vol->csize[i];
    }
  
  gfunc3_compute_xmin_xmax (vol);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
RecParams_apply_to_proj_image (RecParams const *rec_p, gfunc3 *proj_img)
{
  CAPTURE_NULL_VOID (rec_p);
  CAPTURE_NULL_VOID (proj_img);
  
  if (!GFUNC_IS_2D (proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "proj_img must be 2-dimensional");
      return;
    }
  
  if (rec_p->detector_px_size[0] != 0.0)  /* Detector pixel size is set in config file */
    gfunc3_set_csize (proj_img, rec_p->detector_px_size);
    
  gfunc3_compute_xmin_xmax (proj_img);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

float
approx_ctf_xval (float yval, RecParams const *rec_p, float xstart0, float xstart1, float *pdf)
{
  int n;
  float xn_1 = xstart0, xn = xstart1, fval, dfval, update;

  for (n = 0; n < MAX_STEPS; n++)
    {
      fval = ctf_unscaled_radial (xn, rec_p) - yval;
      dfval = (ctf_unscaled_radial (xn, rec_p) - ctf_unscaled_radial (xn_1, rec_p)) / (xn - xn_1);

      update = -fval / dfval;

      xn_1 = xn;
      xn += update;

      if (fabsf (update / (xstart1 - xstart0)) < 1E-3)
        break;
    }

  if (n == MAX_STEPS)
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
compute_xover_spline_params (RecParams *rec_p)
{
  /* This function computes the coefficients of the interpolating splines between the points of the 
   * function graph of 1/CTF where the absolute value of this function reaches the value 
   * rec_p->ctf_trunc. This way, instead of truncating, a smooth transition between the function 
   * values is guaranteed.
   * A secant algorithm is used here to approximate these values. Starting from the zero of the 
   * CTF (can be computed analytically using a and b below), the x value is approximated where 
   * the CTF absolute value is 1/rec_p->ctf_trunc. 
   * This procedure is repeated for NUM_CTF_LOBES zeros of the CTF. The sign change depends on the
   * parity of the zero.
   */

  int i;
  float a, b, c;

  float x0, x1, dx, x0start, x1start, df0, df1, f1, f0, yval = 1.0 / rec_p->ctf_trunc;

  float u0, v0, w0;

  a = rec_p->cs / (4 * rec_p->wave_number * rec_p->wave_number * rec_p->wave_number);
  b = rec_p->defocus_nominal / (2 * rec_p->wave_number);

  for (i = 0; i < NUM_CTF_LOBES; i++)
    {
      c = atan (rec_p->acr) - (i + 1) * M_PI;
      /* (i+1)'th zero = x0start = starting point of secant iteration */
      x0start = (b - sqrtf (b * b + 4 * a * c)) / (2 * a);
      /* x1start = second starting point, a little bit into the correct direction */
      x1start = (i % 2) ? 1.005 * x0start : 0.995 * x0start;  
      /* sign of target value depends on the parity of the zero */
      yval = (i % 2) ? -1.0 / rec_p->ctf_trunc : 1.0 / rec_p->ctf_trunc;  

      /* Compute x0 with CTF(x0) = yval; df0 ~ CTF'(x0) */
      x0 = approx_ctf_xval (yval, rec_p, x0start, x1start, &df0); 
      /* Make interval [x0,x1] symmetric around x0start;
       * The last spline decays to zero with zero derivative.
       */
      x1 = (i < NUM_CTF_LOBES - 1) ? 2 * x0start - x0 : x0start;  

      f0 = 1.0 / yval;
      f1 = (i < NUM_CTF_LOBES - 1) ? 1.0 / ctf_unscaled_radial (x1, rec_p) : 0.0;

      df0 = -df0 / (yval * yval); // (f^-1)'(x) = -f'(x)/f(x)^2
      df1 = (i < NUM_CTF_LOBES - 1) ?
        (1.0 / ctf_unscaled_radial (x1, rec_p) - 1.0 / ctf_unscaled_radial (0.999 * x1, rec_p)) 
        / (0.001 * x1) 
        : 0.0;

      /* Compute coefficients of the spline
       * p(x) = a0 + a1 * (x - x0) + a2 * (x - x0)^2 + a3 * (x - x0)^3, x0 <= x <= x1 
       */
      dx = x1 - x0;
      u0 = (f1 - f0) / dx;
      v0 = (df1 - df0) / dx;
      w0 = (u0 - df0) / dx;
      rec_p->xover_cspline_coeff[i][3] = (v0 - 2 * w0) / dx;  /* a3 */
      rec_p->xover_cspline_coeff[i][2] = 3 * w0 - v0;         /* a2 */
      rec_p->xover_cspline_coeff[i][1] = df0;                 /* a1 */
      rec_p->xover_cspline_coeff[i][0] = f0;                  /* a0 */
      rec_p->xover_pts[i][0] = x0;                            /* x0 */
      rec_p->xover_pts[i][1] = x1;                            /* x1 */
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/
