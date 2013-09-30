/*
 * operators.c -- operations on grid functions
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
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CException.h"

#include "matvec3.h"
#include "misc.h"

#include "fft.h"
#include "gfunc3.h"
#include "mrc.h"
#include "operators.h"

#include "ai_options.h"

// TODO: descriptions

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection (gfunc3 const *proj_img, vec3 const angles_deg, vec3 const axis_shift_px, 
                     gfunc3 *volume)
{
  CEXCEPTION_T e = EXC_NONE;
  int ix, iy, iz;
  size_t idx = 0;
  vec3 omega_x, omega_y, omega;
  vec3 xm;
  float Pxmin[2], Pdx[2], Pdy[2], Pdz[2], Pv[3], Pvx0[2], Pvy0[2];
  
  CAPTURE_NULL (proj_img);
  CAPTURE_NULL (angles_deg);
  CAPTURE_NULL (volume);
  GFUNC_CHECK_INIT_STATUS (proj_img);
  GFUNC_CHECK_INIT_STATUS (volume);
  
  compute_rotated_basis (angles_deg, omega_x, omega_y, omega);

  /* The projections of the increments in x, y and z directions as well as xmin are precomputed for 
   * speed reasons.
   * The letter "P" stands for "projected"
   */
  vec3_copy (xm, volume->xmin);
  vec3_axpby (1, xm, -1, proj_img->x0);
  Pxmin[0] = vec3_dot (omega_x, xm);
  Pxmin[1] = vec3_dot (omega_y, xm);
  // TODO: include tilt axis shift

  Pdx[0] = omega_x[0] * volume->csize[0];
  Pdx[1] = omega_y[0] * volume->csize[0];

  Pdy[0] = omega_x[1] * volume->csize[1];
  Pdy[1] = omega_y[1] * volume->csize[1];

  Pdz[0] = omega_x[2] * volume->csize[2];
  Pdz[1] = omega_y[2] * volume->csize[2];

  Pv[0] = Pxmin[0];
  Pv[1] = Pxmin[1];
  Pv[2] = 0.0;

  Try
  {
    for (iz = 0; iz < volume->shape[2]; iz++)
      {
        Pvy0[0] = Pv[0];
        Pvy0[1] = Pv[1];

        for (iy = 0; iy < volume->shape[1]; iy++)
          {
            Pvx0[0] = Pv[0];
            Pvx0[1] = Pv[1];

            for (ix = 0; ix < volume->shape[0]; ix++, idx++)
              {
                volume->fvals[idx] += gfunc3_interp_linear_2d (proj_img, Pv);

                Pv[0] += Pdx[0];
                Pv[1] += Pdx[1];
              }
            Pv[0] = Pvx0[0] + Pdy[0];
            Pv[1] = Pvx0[1] + Pdy[1];
          }
        Pv[0] = Pvy0[0] + Pdz[0];
        Pv[1] = Pvy0[1] + Pdz[1];
      }
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_sax (gfunc3 const *proj_img, float const theta_deg, float const axis_shift_y_px, 
                         gfunc3 *volume)
{
  CEXCEPTION_T e = EXC_NONE;
  int ix, iy, iz;
  size_t idx = 0l, fiy, fi;
  float Pxmin_y, Pdy, Pdz, Pvy, Pvy0;
  float cos_theta, sin_theta;
  float x, *wlx = NULL, *wux = NULL, idxf_x, wly, wuy, idxf_y;
  int *idx_x = NULL, idx_y;
  float fv;
  
  CAPTURE_NULL (proj_img);
  CAPTURE_NULL (volume);
  GFUNC_CHECK_INIT_STATUS (proj_img);
  GFUNC_CHECK_INIT_STATUS (volume);
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  if (autocenter_vol_flag)
    {
      volume->x0[0] = proj_img->x0[0];
      volume->x0[1] = proj_img->x0[1] + axis_shift_y_px * proj_img->csize[1];
      gfunc3_compute_xmin_xmax (volume);
    }

  Pxmin_y = cos_theta * (volume->xmin[1] - axis_shift_y_px * proj_img->csize[1]) 
    + sin_theta * volume->xmin[2];

  Pdy = cos_theta * volume->csize[1];
  Pdz = sin_theta * volume->csize[2];

  Try
  {
    /* Precompute weights and indices for the 'x' component since they don't change across the loop */
    wlx = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    wux = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    idx_x = (int *) ali16_malloc (volume->shape[0] * sizeof (int));
  }  CATCH_RETURN_VOID (e);

  x = volume->xmin[0];
  for (ix = 0; ix < volume->shape[0]; ix++)
    {
      if ((x <= proj_img->xmin[0]) || (x >= proj_img->xmax[0]))
        {
          idx_x[ix] = -1;
          x += volume->csize[0];
          continue;
        }
        
      idxf_x    = (x - proj_img->xmin[0]) / proj_img->csize[0];
      idx_x[ix] = (int) idxf_x;
      wux[ix]   = idxf_x - idx_x[ix];
      wlx[ix]   = 1.0 - wux[ix];
      
      x += volume->csize[0];
    }

  
  Pvy = Pxmin_y;
  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvy0 = Pvy;
      
      for (iy = 0; iy < volume->shape[1]; iy++)
        {
          if ((Pvy <= proj_img->xmin[1]) || (Pvy >= proj_img->xmax[1]))
            {
              Pvy += Pdy;
              idx += volume->shape[0];
              continue;
            }
          
          idxf_y = (Pvy - proj_img->xmin[1]) / proj_img->csize[1];
          idx_y = (int) idxf_y;
          wuy = idxf_y - idx_y;
          wly = 1.0f - wuy;
          
          fiy = idx_y * proj_img->shape[0];

          for (ix = 0; ix < volume->shape[0]; ix++, idx++)
            {
              if (idx_x[ix] == -1)
                continue;

              fi = fiy + idx_x[ix];
              
              fv  = (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wly;
              fi += proj_img->shape[0];
              fv += (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wuy;

              volume->fvals[idx] += fv;
            }
          Pvy += Pdy;
        }
      Pvy = Pvy0 + Pdz;
    }

  free (wlx);
  free (wux);
  free (idx_x);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_convolution (gfunc3 *gf1, gfunc3 *gf2)
{
  CEXCEPTION_T e = EXC_NONE;
  int gf1_was_hc, gf2_was_hc;
  float factor;

  CAPTURE_NULL (gf1);
  CAPTURE_NULL (gf2);
  GFUNC_CHECK_INIT_STATUS (gf1);
  GFUNC_CHECK_INIT_STATUS (gf2);

  Try
  {
    if (GFUNC_IS_2D (gf1))
      factor = 1.0 / (2 * M_PI);
    else
      factor = 1.0 / (2 * M_PI * M_SQRT2PI);

    gf1_was_hc = gf1->is_halfcomplex;
    gf2_was_hc = gf2->is_halfcomplex;

    if (!gf1->is_halfcomplex)
      fft_forward (gf1);
      
    if (!gf2->is_halfcomplex)
      fft_forward (gf2);

    gfunc3_mul (gf1, gf2);

    if (!gf1_was_hc)
      fft_backward (gf1);

    if (!gf2_was_hc)
      fft_backward (gf2);

    gfunc3_scale (gf1, factor);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
image_rotation (gfunc3 *proj_img, float const psi_deg)
{
  CEXCEPTION_T e = EXC_NONE;
  int ix, iy;
  int Nx_new, Ny_new;
  size_t idx = 0, ntotal_new;
  float cos_psi, sin_psi;
  float dx_rot[2], dy_rot[2];
  float xlen, ylen, xlen_new, ylen_new, xmin_new, ymin_new;
  float vx0[2];
  vec3 v = {0.0, 0.0, 0.0}, img_x0;
  float *fvals_rot = NULL;
  
  CAPTURE_NULL (proj_img);
  GFUNC_CHECK_INIT_STATUS (proj_img);
  
  if (!GFUNC_IS_2D(proj_img))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 2d functions supported.");

  if (fabsf (psi_deg) > 180.0)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Rotation angle must be between -180 and +180 degrees");

  if (psi_deg == 0.0)
    return;

  PRINT_VERBOSE ("psi_deg = %f\n", psi_deg);

  /* Image is rotated counter-clockwise with angle psi_deg, so the original function is evaluated
   * (interpolated) at the inversely rotated points.
   * 
   * Use exact sine and cosine values in the four corner cases
   */
  
  if      (psi_deg ==  90.0)          { cos_psi =  0.0; sin_psi =  1.0; }
  else if (psi_deg == -90.0)          { cos_psi =  0.0; sin_psi = -1.0; }
  else if (fabsf (psi_deg == 180.0))  { cos_psi = -1.0; sin_psi =  0.0; }
  else
    {
      cos_psi = cosf (psi_deg * ONE_DEGREE);
      sin_psi = sinf (psi_deg * ONE_DEGREE);
    }
    
  /* Shift image to origin for intrinsic rotation; origin is stored for backup */
  vec3_copy (img_x0, proj_img->x0);
  vec3_set_all (proj_img->x0, 0.0);
  gfunc3_compute_xmin_xmax (proj_img);
  
  /* Enlarge the grid such that all rotated points are inside */
  xlen = proj_img->xmax[0] - proj_img->xmin[0];
  ylen = proj_img->xmax[1] - proj_img->xmin[1];

  xlen_new = xlen * fabsf (cos_psi) + ylen * fabsf (sin_psi);
  ylen_new = xlen * fabsf (sin_psi) + ylen * fabsf (cos_psi);
  Nx_new = (int) (ceilf (xlen_new / proj_img->csize[0])) + 1;
  Ny_new = (int) (ceilf (ylen_new / proj_img->csize[1])) + 1;
  ntotal_new = (size_t) Nx_new * Ny_new;
 
  xmin_new = proj_img->x0[0] -  Nx_new / 2 * proj_img->csize[0];
  ymin_new = proj_img->x0[1] -  Ny_new / 2 * proj_img->csize[1];
  
  v[0] = cos_psi * xmin_new - sin_psi * ymin_new;
  v[1] = sin_psi * xmin_new + cos_psi * ymin_new;
  
  /* Rotated increments */
  dx_rot[0] =  cos_psi * proj_img->csize[0];
  dx_rot[1] =  sin_psi * proj_img->csize[0];
  dy_rot[0] = -sin_psi * proj_img->csize[1];
  dy_rot[1] =  cos_psi * proj_img->csize[1];

  /* Initialize new array */
  Try { fvals_rot = (float *) ali16_malloc (ntotal_new * sizeof (float)); }  CATCH_RETURN_VOID (e);
  for (idx = 0; idx < ntotal_new; idx++) fvals_rot[idx] = 0.0;


  /* TODO: inline the interpolation */
  for (iy = 0, idx = 0; iy < Ny_new; iy++)
    {
      vx0[0] = v[0];
      vx0[1] = v[1];
      
      for (ix = 0; ix < Nx_new; ix++, idx++)
        {
          fvals_rot[idx] = gfunc3_interp_linear_2d (proj_img, v);
          v[0] += dx_rot[0];
          v[1] += dx_rot[1];
        }
        
      v[0] = vx0[0] + dy_rot[0];
      v[1] = vx0[1] + dy_rot[1];
    }
    
  free (proj_img->fvals);
  proj_img->fvals     = fvals_rot;
  proj_img->shape[0]  = Nx_new;
  proj_img->shape[1]  = Ny_new;
  proj_img->ntotal    = ntotal_new;
  vec3_copy (proj_img->x0, img_x0);
  gfunc3_compute_xmin_xmax (proj_img);
  
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
histogram_normalization (gfunc3 *proj_img, idx3 bg_ix0, idx3 const bg_shp)
{
  CEXCEPTION_T e = EXC_NONE;
  int i;
  int static count = 0;
  size_t j, *idcs;
  float bg_avg = 0.0, bg_var = 0.0;
  vec3 bg_x0;
  gfunc3 *bg_patch;

  CAPTURE_NULL (proj_img);
  CAPTURE_NULL (bg_ix0);
  CAPTURE_NULL (bg_shp);
  GFUNC_CHECK_INIT_STATUS (proj_img);
  
  count++;
  Try
  {
    bg_patch = new_gfunc3 ();
    
    /* TODO: Use 4 corner patches if bg_ix0[2] == -1 (see options) */
    if (bg_ix0[2] == -1)
      bg_ix0[2] = 0;
    
    /* Init bg_patch grid */
    for (i = 0; i < 3; i++)
      bg_x0[i] = proj_img->xmin[i] + (bg_ix0[i] + bg_shp[i] / 2 ) * proj_img->csize[i];

    gfunc3_init (bg_patch, bg_x0, proj_img->csize, bg_shp, REAL);

    /* Extract values from proj_img */
    idcs = gfunc3_subgrid_flatidcs (proj_img, bg_patch);
    
    for (j = 0; j < bg_patch->ntotal; j++)
      bg_patch->fvals[j] = proj_img->fvals[idcs[j]];
     
    free (idcs);
    
    if (DEBUGGING)
      temp_mrc_out (bg_patch, "bg_patch_", count); 
    
    /* Compute statistics */
    bg_avg = gfunc3_mean (bg_patch);
    bg_var = gfunc3_variance (bg_patch, &bg_avg);
    
    PRINT_VERBOSE ("Background mean    : %f\n", bg_avg);
    PRINT_VERBOSE ("Background variance: %f\n", bg_var);
        
    /* Normalize proj_img -> (proj_img - bg_avg) / sqrt(bg_var) */
    gfunc3_add_constant (proj_img, -bg_avg);
    
    gfunc3_scale (proj_img, 1.0 / sqrtf (bg_var));
    
    gfunc3_free (&bg_patch);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
probability_normalization (gfunc3 *gf)
{
  float integral;
  
  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);
 
  integral = lp_integral (gf, TRAPEZOIDAL);
 
  PRINT_VERBOSE ("Function integral: %f\n", integral);
 
  if (fabsf (integral) > EPS_DENOM)
    gfunc3_scale (gf, 1.0 / integral);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

float
lp_integral (gfunc3 const *gf, integration_rule rule)
{
  size_t i;
  float integral = 0.0;
  
  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);
 
  switch (rule)
    {
      case TRAPEZOIDAL:
        for (i = 0; i < gf->ntotal; i++)
          integral += gf->fvals[i];
        
        integral *= vec3_product (gf->csize);
    }
    
  return integral;
}

/*-------------------------------------------------------------------------------------------------*/

