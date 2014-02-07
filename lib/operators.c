/*
 * operators.c -- operations on grid functions
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
#include <float.h>

#include "CException.h"

#include "vec3.h"
#include "misc.h"

#include "tiltangles.h"

#include "fft.h"
#include "gfunc3.h"
#include "gfunc3_private.h"
#include "mrc.h"
#include "operators.h"

// TODO: descriptions

/*-------------------------------------------------------------------------------------------------*/

float *
perp_plane_freqs (gfunc3 const *ft_proj_img_grid, vec3 const normal_angles_deg)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  float *freqs = NULL;
  
  vec3 p;
  int ix, iy;
  float sin_psi, cos_psi, sin_theta, cos_theta, sin_phi, cos_phi;
  float *cur_freq;

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
          /* p[0] is the coefficient of the first unit vector in omega^\perp, 
           * p[1] the coefficient of the second one. These vectors correspond to the first two
           * columns of the transposed rotation matrix in 'x' convention, see
           * https://de.wikipedia.org/wiki/Eulersche_Winkel
           */
          cur_freq[0] =   p[0] * ( cos_phi * cos_psi - sin_phi * cos_theta * sin_psi)
                        + p[1] * (-cos_phi * sin_psi - sin_phi * cos_theta * cos_psi);
          cur_freq[1] =   p[0] * ( sin_phi * cos_psi + cos_phi * cos_theta * sin_psi)
                        + p[1] * (-sin_phi * sin_psi + cos_phi * cos_theta * cos_psi);
          cur_freq[2] =   p[0] * ( sin_theta * sin_psi)
                        + p[1] * ( sin_theta * cos_psi);
          
          p[0] += ft_proj_img_grid->csize[0];
        }
      p[0]  = ft_proj_img_grid->xmin[0];
      p[1] += ft_proj_img_grid->csize[1];
    }

  return freqs;
}

/*-------------------------------------------------------------------------------------------------*/

float *
perp_plane_stack_freqs (gfunc3 const *ft_proj_img_grid, tiltangles const *tilts)
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
        freqs = perp_plane_freqs (ft_proj_img_grid, angles_deg);
      } CATCH_RETURN (_e, NULL);
      
      for (j = 0; j < 3 * ft_proj_img_grid->ntotal; j++)
        *(cur_freq++) = freqs[j];
    }
  
  free (freqs);
  
  return allfreqs;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_projection (gfunc3 const *volume, vec3 const angles_deg, gfunc3 *proj_img)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  float *freqs = NULL;
  
  CAPTURE_NULL_VOID (volume);
  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (angles_deg);
  
  GFUNC_CAPTURE_UNINIT_VOID (volume);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);

  if (proj_img->type != COMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Projection image must be of COMPLEX type.");
      return;
    }

  Try { gfunc3_grid_fwd_reciprocal (proj_img); }  CATCH_RETURN_VOID (_e);

  Try { freqs = perp_plane_freqs (proj_img, angles_deg); }  CATCH_RETURN_VOID (_e);
  
  Try { 
    nfft_transform (volume, freqs, proj_img->ntotal, (float complex **) (&proj_img->fvals)); 
  } CATCH_RETURN_VOID (_e);
  
  Try {
    fft_backward (proj_img);
    gfunc3_scale (proj_img, M_SQRT2PI);
  } CATCH_RETURN_VOID (_e);
  
  free (freqs);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_r_r (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume, float weight)
{
  int ix, iy, iz;
  size_t idx = 0;
  vec3 omega_x, omega_y, omega;
  vec3 xm;
  float Pxmin[2], Pdx[2], Pdy[2], Pdz[2], Pv[3], Pvx0[2], Pvy0[2];
  float *cur_fval = volume->fvals;
  
  compute_rotated_basis (angles_deg, omega_x, omega_y, omega);

  /* The projections of the increments in x, y and z directions as well as xmin are precomputed for 
   * speed reasons.
   * The letter "P" stands for "projected"
   */
  vec3_copy (xm, volume->xmin);
  vec3_axpby (1, xm, -1, proj_img->x0);
  Pxmin[0] = vec3_dot (omega_x, xm);
  Pxmin[1] = vec3_dot (omega_y, xm);

  Pdx[0] = omega_x[0] * volume->csize[0];
  Pdx[1] = omega_y[0] * volume->csize[0];

  Pdy[0] = omega_x[1] * volume->csize[1];
  Pdy[1] = omega_y[1] * volume->csize[1];

  Pdz[0] = omega_x[2] * volume->csize[2];
  Pdz[1] = omega_y[2] * volume->csize[2];

  Pv[0] = Pxmin[0];
  Pv[1] = Pxmin[1];
  Pv[2] = 0.0;


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
              *(cur_fval++) += weight * gfunc3_interp_linear_2d_r (proj_img, Pv);

              Pv[0] += Pdx[0];
              Pv[1] += Pdx[1];
            }
          Pv[0] = Pvx0[0] + Pdy[0];
          Pv[1] = Pvx0[1] + Pdy[1];
        }
      Pv[0] = Pvy0[0] + Pdz[0];
      Pv[1] = Pvy0[1] + Pdz[1];
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_r_c (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume, float weight)
{
  int ix, iy, iz;
  vec3 omega_x, omega_y, omega;
  vec3 xm;
  float Pxmin[2], Pdx[2], Pdy[2], Pdz[2], Pv[3], Pvx0[2], Pvy0[2];
  float complex *cur_fval = (float complex *) volume->fvals;
  
  compute_rotated_basis (angles_deg, omega_x, omega_y, omega);

  /* The projections of the increments in x, y and z directions as well as xmin are precomputed for 
   * speed reasons.
   * The letter "P" stands for "projected"
   */
  vec3_copy (xm, volume->xmin);
  vec3_axpby (1, xm, -1, proj_img->x0);
  Pxmin[0] = vec3_dot (omega_x, xm);
  Pxmin[1] = vec3_dot (omega_y, xm);

  Pdx[0] = omega_x[0] * volume->csize[0];
  Pdx[1] = omega_y[0] * volume->csize[0];

  Pdy[0] = omega_x[1] * volume->csize[1];
  Pdy[1] = omega_y[1] * volume->csize[1];

  Pdz[0] = omega_x[2] * volume->csize[2];
  Pdz[1] = omega_y[2] * volume->csize[2];

  Pv[0] = Pxmin[0];
  Pv[1] = Pxmin[1];
  Pv[2] = 0.0;


  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvy0[0] = Pv[0];
      Pvy0[1] = Pv[1];

      for (iy = 0; iy < volume->shape[1]; iy++)
        {
          Pvx0[0] = Pv[0];
          Pvx0[1] = Pv[1];

          for (ix = 0; ix < volume->shape[0]; ix++)
            {
              *(cur_fval++) += weight * gfunc3_interp_linear_2d_r (proj_img, Pv);

              Pv[0] += Pdx[0];
              Pv[1] += Pdx[1];
            }
          Pv[0] = Pvx0[0] + Pdy[0];
          Pv[1] = Pvx0[1] + Pdy[1];
        }
      Pv[0] = Pvy0[0] + Pdz[0];
      Pv[1] = Pvy0[1] + Pdz[1];
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_c_c (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume, float weight)
{
  int ix, iy, iz;
  vec3 omega_x, omega_y, omega;
  vec3 xm;
  float Pxmin[2], Pdx[2], Pdy[2], Pdz[2], Pv[3], Pvx0[2], Pvy0[2];
  float complex *cur_fval = (float complex *) volume->fvals;
  
  compute_rotated_basis (angles_deg, omega_x, omega_y, omega);

  /* The projections of the increments in x, y and z directions as well as xmin are precomputed for 
   * speed reasons.
   * The letter "P" stands for "projected"
   */
  vec3_copy (xm, volume->xmin);
  vec3_axpby (1, xm, -1, proj_img->x0);
  Pxmin[0] = vec3_dot (omega_x, xm);
  Pxmin[1] = vec3_dot (omega_y, xm);

  Pdx[0] = omega_x[0] * volume->csize[0];
  Pdx[1] = omega_y[0] * volume->csize[0];

  Pdy[0] = omega_x[1] * volume->csize[1];
  Pdy[1] = omega_y[1] * volume->csize[1];

  Pdz[0] = omega_x[2] * volume->csize[2];
  Pdz[1] = omega_y[2] * volume->csize[2];

  Pv[0] = Pxmin[0];
  Pv[1] = Pxmin[1];
  Pv[2] = 0.0;


  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvy0[0] = Pv[0];
      Pvy0[1] = Pv[1];

      for (iy = 0; iy < volume->shape[1]; iy++)
        {
          Pvx0[0] = Pv[0];
          Pvx0[1] = Pv[1];

          for (ix = 0; ix < volume->shape[0]; ix++)
            {
              *(cur_fval++) += weight * gfunc3_interp_linear_2d_c (proj_img, Pv);

              Pv[0] += Pdx[0];
              Pv[1] += Pdx[1];
            }
          Pv[0] = Pvx0[0] + Pdy[0];
          Pv[1] = Pvx0[1] + Pdy[1];
        }
      Pv[0] = Pvy0[0] + Pdz[0];
      Pv[1] = Pvy0[1] + Pdz[1];
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume, float weight)
{
  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (angles_deg);
  CAPTURE_NULL_VOID (volume);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);
  GFUNC_CAPTURE_UNINIT_VOID (volume);
  
  if (!GFUNC_IS_2D(proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Projection image must be 2-dimensional.");
      return;
    }

  if ((proj_img->type == HALFCOMPLEX) || (volume->type == HALFCOMPLEX))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Functions must not be of HALFCOMPLEX type.");
      return;
    }
  
  if ((proj_img->type == COMPLEX) && (volume->type == REAL))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Volume must be COMPLEX if projection image is.");
      return;
    }
  
  if (volume->type == REAL) 
    xray_backprojection_r_r (proj_img, angles_deg, volume, weight);
  else if ((volume->type == COMPLEX) && (proj_img == REAL))
    xray_backprojection_r_c (proj_img, angles_deg, volume, weight);
  else 
    xray_backprojection_c_c (proj_img, angles_deg, volume, weight);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

// extern int autocenter_vol_flag;

void
xray_backprojection_sax_r_r (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_y_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t fiy, fi;
  float Pxmin_y, Pdy, Pdz, Pvy, Pvy0;
  float cos_theta, sin_theta;
  float x, *wlx = NULL, *wux = NULL, idxf_x, wly, wuy, idxf_y;
  int *idx_x = NULL, idx_y;
  float fv;
  float *cur_fval = volume->fvals;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  /* TODO: find a solution to handle this (undefined reference now) */
  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0];
      // volume->x0[1] = proj_img->x0[1] + axis_shift_y_px * proj_img->csize[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_y = cos_theta * (volume->xmin[1] - axis_shift_y_px * proj_img->csize[1]) 
    + sin_theta * volume->xmin[2];

  Pdy = cos_theta * volume->csize[1];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'x' component since they don't change across the loop */
    wlx = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    wux = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    idx_x = (int *) ali16_malloc (volume->shape[0] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

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
              cur_fval += volume->shape[0];
              continue;
            }
          
          idxf_y = (Pvy - proj_img->xmin[1]) / proj_img->csize[1];
          idx_y = (int) idxf_y;
          wuy = idxf_y - idx_y;
          wly = 1.0f - wuy;
          
          fiy = idx_y * proj_img->shape[0];

          for (ix = 0; ix < volume->shape[0]; ix++, cur_fval++)
            {
              if (idx_x[ix] == -1)
                continue;

              fi = fiy + idx_x[ix];
              
              fv  = (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wly;
              fi += proj_img->shape[0];
              fv += (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wuy;

              *cur_fval += weight * fv;
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

// extern int autocenter_vol_flag;

void
xray_backprojection_sax_r_c (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_y_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t fiy, fi;
  float Pxmin_y, Pdy, Pdz, Pvy, Pvy0;
  float cos_theta, sin_theta;
  float x, *wlx = NULL, *wux = NULL, idxf_x, wly, wuy, idxf_y;
  int *idx_x = NULL, idx_y;
  float fv;
  float complex *cur_fval = (float complex *) volume->fvals;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  /* TODO: find a solution to handle this (undefined reference now) */
  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0];
      // volume->x0[1] = proj_img->x0[1] + axis_shift_y_px * proj_img->csize[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_y = cos_theta * (volume->xmin[1] - axis_shift_y_px * proj_img->csize[1]) 
    + sin_theta * volume->xmin[2];

  Pdy = cos_theta * volume->csize[1];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'x' component since they don't change across the loop */
    wlx = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    wux = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    idx_x = (int *) ali16_malloc (volume->shape[0] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

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
              cur_fval += volume->shape[0];
              continue;
            }
          
          idxf_y = (Pvy - proj_img->xmin[1]) / proj_img->csize[1];
          idx_y = (int) idxf_y;
          wuy = idxf_y - idx_y;
          wly = 1.0f - wuy;
          
          fiy = idx_y * proj_img->shape[0];

          for (ix = 0; ix < volume->shape[0]; ix++, cur_fval++)
            {
              if (idx_x[ix] == -1)
                continue;

              fi = fiy + idx_x[ix];
              
              fv  = (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wly;
              fi += proj_img->shape[0];
              fv += (wlx[ix] * proj_img->fvals[fi] + wux[ix] * proj_img->fvals[fi + 1]) * wuy;

              *cur_fval += weight * fv;
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

// extern int autocenter_vol_flag;

void
xray_backprojection_sax_c_c (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_y_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t fiy, fi;
  float Pxmin_y, Pdy, Pdz, Pvy, Pvy0;
  float cos_theta, sin_theta;
  float x, *wlx = NULL, *wux = NULL, idxf_x, wly, wuy, idxf_y;
  int *idx_x = NULL, idx_y;
  float fv;
  float complex *cur_fval = (float complex *) volume->fvals;
  float complex *proj_fvals_c = (float complex *) proj_img->fvals;;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  /* TODO: find a solution to handle this (undefined reference now) */
  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0];
      // volume->x0[1] = proj_img->x0[1] + axis_shift_y_px * proj_img->csize[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_y = cos_theta * (volume->xmin[1] - axis_shift_y_px * proj_img->csize[1]) 
    + sin_theta * volume->xmin[2];

  Pdy = cos_theta * volume->csize[1];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'x' component since they don't change across the loop */
    wlx = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    wux = (float *) ali16_malloc (volume->shape[0] * sizeof (float));
    idx_x = (int *) ali16_malloc (volume->shape[0] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

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
              cur_fval += volume->shape[0];
              continue;
            }
          
          idxf_y = (Pvy - proj_img->xmin[1]) / proj_img->csize[1];
          idx_y = (int) idxf_y;
          wuy = idxf_y - idx_y;
          wly = 1.0f - wuy;
          
          fiy = idx_y * proj_img->shape[0];

          for (ix = 0; ix < volume->shape[0]; ix++, cur_fval++)
            {
              if (idx_x[ix] == -1)
                continue;

              fi = fiy + idx_x[ix];
              
              fv  = (wlx[ix] * proj_fvals_c[fi] + wux[ix] * proj_fvals_c[fi + 1]) * wly;
              fi += proj_img->shape[0];
              fv += (wlx[ix] * proj_fvals_c[fi] + wux[ix] * proj_fvals_c[fi + 1]) * wuy;

              *cur_fval += weight * fv;
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
xray_backprojection_say_r_r (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_x_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t idx0 = 0l, fi;
  float Pxmin_x, Pdx, Pdz, Pvx, Pvx0;
  float cos_theta, sin_theta;
  float y, *wly = NULL, *wuy = NULL, idxf_y, wlx, wux, idxf_x;
  int *idx_y = NULL, idx_x;
  float fv;
  float *cur_fval = volume->fvals;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0] - axis_shift_x_px * proj_img->csize[0];
      // volume->x0[1] = proj_img->x0[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_x = cos_theta * (volume->xmin[0] - axis_shift_x_px * proj_img->csize[0]) 
    + sin_theta * volume->xmin[2];

  Pdx = cos_theta * volume->csize[0];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'y' component since they don't change across the loop */
    wly = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    wuy = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    idx_y = (int *) ali16_malloc (volume->shape[1] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

  y = volume->xmin[1];
  for (iy = 0; iy < volume->shape[1]; iy++)
    {
      if ((y <= proj_img->xmin[1]) || (y >= proj_img->xmax[1]))
        {
          idx_y[iy] = -1;
          y += volume->csize[1];
          continue;
        }
        
      idxf_y    = (y - proj_img->xmin[1]) / proj_img->csize[1];
      idx_y[iy] = (int) idxf_y;
      wuy[iy]   = idxf_y - idx_y[iy];
      wly[iy]   = 1.0 - wuy[iy];
      
      y += volume->csize[1];
    }

  
  Pvx = Pxmin_x;
  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvx0 = Pvx;
      idx0 = iz * volume->shape[1] * volume->shape[0];
      
      for (ix = 0; ix < volume->shape[0]; ix++)
        {
          if ((Pvx <= proj_img->xmin[0]) || (Pvx >= proj_img->xmax[0]))
            {
              Pvx += Pdx;
              idx0++;
              continue;
            }
          
          idxf_x = (Pvx - proj_img->xmin[0]) / proj_img->csize[0];
          idx_x = (int) idxf_x;
          wux = idxf_x - idx_x;
          wlx = 1.0f - wux;
          
          for (iy = 0; iy < volume->shape[1]; iy++, cur_fval += volume->shape[0])
            {
              if (idx_y[iy] == -1)
                continue;

              fi = idx_x + idx_y[iy] * proj_img->shape[0];
              
              fv  = (wlx * proj_img->fvals[fi] + wux * proj_img->fvals[fi + 1]) * wly[iy];
              fi += proj_img->shape[0];
              fv += (wlx * proj_img->fvals[fi] + wux * proj_img->fvals[fi + 1]) * wuy[iy];

              *cur_fval += weight * fv;
            }
          cur_fval = &volume->fvals[idx0++];
          Pvx += Pdx;
        }
      Pvx = Pvx0 + Pdz;
    }

  free (wly);
  free (wuy);
  free (idx_y);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_say_r_c (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_x_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t idx0 = 0l, fi;
  float Pxmin_x, Pdx, Pdz, Pvx, Pvx0;
  float cos_theta, sin_theta;
  float y, *wly = NULL, *wuy = NULL, idxf_y, wlx, wux, idxf_x;
  int *idx_y = NULL, idx_x;
  float fv;
  float complex *cur_fval = (float complex *) volume->fvals;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0] - axis_shift_x_px * proj_img->csize[0];
      // volume->x0[1] = proj_img->x0[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_x = cos_theta * (volume->xmin[0] - axis_shift_x_px * proj_img->csize[0]) 
    + sin_theta * volume->xmin[2];

  Pdx = cos_theta * volume->csize[0];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'y' component since they don't change across the loop */
    wly = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    wuy = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    idx_y = (int *) ali16_malloc (volume->shape[1] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

  y = volume->xmin[1];
  for (iy = 0; iy < volume->shape[1]; iy++)
    {
      if ((y <= proj_img->xmin[1]) || (y >= proj_img->xmax[1]))
        {
          idx_y[iy] = -1;
          y += volume->csize[1];
          continue;
        }
        
      idxf_y    = (y - proj_img->xmin[1]) / proj_img->csize[1];
      idx_y[iy] = (int) idxf_y;
      wuy[iy]   = idxf_y - idx_y[iy];
      wly[iy]   = 1.0 - wuy[iy];
      
      y += volume->csize[1];
    }

  
  Pvx = Pxmin_x;
  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvx0 = Pvx;
      idx0 = iz * volume->shape[1] * volume->shape[0];
      
      for (ix = 0; ix < volume->shape[0]; ix++)
        {
          if ((Pvx <= proj_img->xmin[0]) || (Pvx >= proj_img->xmax[0]))
            {
              Pvx += Pdx;
              idx0++;
              continue;
            }
          
          idxf_x = (Pvx - proj_img->xmin[0]) / proj_img->csize[0];
          idx_x = (int) idxf_x;
          wux = idxf_x - idx_x;
          wlx = 1.0f - wux;
          
          for (iy = 0; iy < volume->shape[1]; iy++, cur_fval += volume->shape[0])
            {
              if (idx_y[iy] == -1)
                continue;

              fi = idx_x + idx_y[iy] * proj_img->shape[0];
              
              fv  = (wlx * proj_img->fvals[fi] + wux * proj_img->fvals[fi + 1]) * wly[iy];
              fi += proj_img->shape[0];
              fv += (wlx * proj_img->fvals[fi] + wux * proj_img->fvals[fi + 1]) * wuy[iy];

              *cur_fval += weight * fv;
            }
          cur_fval = (float complex *) &volume->fvals[2 * idx0++];
          Pvx += Pdx;
        }
      Pvx = Pvx0 + Pdz;
    }

  free (wly);
  free (wuy);
  free (idx_y);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_say_c_c (gfunc3 const *proj_img, float const theta_deg, 
                             float const axis_shift_x_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t idx0 = 0l, fi;
  float Pxmin_x, Pdx, Pdz, Pvx, Pvx0;
  float cos_theta, sin_theta;
  float y, *wly = NULL, *wuy = NULL, idxf_y, wlx, wux, idxf_x;
  int *idx_y = NULL, idx_x;
  float fv;
  float complex *cur_fval = (float complex *) volume->fvals;
  float complex *proj_fvals_c = (float complex *) proj_img->fvals;
  
  /* The projections of the increments in y and z directions as well as xmin are precomputed */
  cos_theta = cosf (theta_deg * ONE_DEGREE);
  sin_theta = sinf (theta_deg * ONE_DEGREE);

  // if (autocenter_vol_flag)
    // {
      // volume->x0[0] = proj_img->x0[0] - axis_shift_x_px * proj_img->csize[0];
      // volume->x0[1] = proj_img->x0[1];
      // gfunc3_compute_xmin_xmax (volume);
    // }

  Pxmin_x = cos_theta * (volume->xmin[0] - axis_shift_x_px * proj_img->csize[0]) 
    + sin_theta * volume->xmin[2];

  Pdx = cos_theta * volume->csize[0];
  Pdz = sin_theta * volume->csize[2];

  Try {
    /* Precompute weights and indices for the 'y' component since they don't change across the loop */
    wly = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    wuy = (float *) ali16_malloc (volume->shape[1] * sizeof (float));
    idx_y = (int *) ali16_malloc (volume->shape[1] * sizeof (int));
  }  CATCH_RETURN_VOID (_e);

  y = volume->xmin[1];
  for (iy = 0; iy < volume->shape[1]; iy++)
    {
      if ((y <= proj_img->xmin[1]) || (y >= proj_img->xmax[1]))
        {
          idx_y[iy] = -1;
          y += volume->csize[1];
          continue;
        }
        
      idxf_y    = (y - proj_img->xmin[1]) / proj_img->csize[1];
      idx_y[iy] = (int) idxf_y;
      wuy[iy]   = idxf_y - idx_y[iy];
      wly[iy]   = 1.0 - wuy[iy];
      
      y += volume->csize[1];
    }

  
  Pvx = Pxmin_x;
  for (iz = 0; iz < volume->shape[2]; iz++)
    {
      Pvx0 = Pvx;
      idx0 = iz * volume->shape[1] * volume->shape[0];
      
      for (ix = 0; ix < volume->shape[0]; ix++)
        {
          if ((Pvx <= proj_img->xmin[0]) || (Pvx >= proj_img->xmax[0]))
            {
              Pvx += Pdx;
              idx0++;
              continue;
            }
          
          idxf_x = (Pvx - proj_img->xmin[0]) / proj_img->csize[0];
          idx_x = (int) idxf_x;
          wux = idxf_x - idx_x;
          wlx = 1.0f - wux;
          
          for (iy = 0; iy < volume->shape[1]; iy++, cur_fval += volume->shape[0])
            {
              if (idx_y[iy] == -1)
                continue;

              fi = idx_x + idx_y[iy] * proj_img->shape[0];
              
              fv  = (wlx * proj_fvals_c[fi] + wux * proj_fvals_c[fi + 1]) * wly[iy];
              fi += proj_img->shape[0];
              fv += (wlx * proj_fvals_c[fi] + wux * proj_fvals_c[fi + 1]) * wuy[iy];

              *cur_fval += weight * fv;
            }
          cur_fval = (float complex *) &volume->fvals[2 * idx0++];
          Pvx += Pdx;
        }
      Pvx = Pvx0 + Pdz;
    }

  free (wly);
  free (wuy);
  free (idx_y);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
xray_backprojection_single_axis (gfunc3 const *proj_img, float theta_deg, int axis,
                                 float tilt_axis_par_shift_px, gfunc3 *volume, float weight)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (volume);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);
  GFUNC_CAPTURE_UNINIT_VOID (volume);
  
  if (!GFUNC_IS_2D(proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Projection image must be 2-dimensional.");
      return;
    }

  if ((proj_img->type == HALFCOMPLEX) || (volume->type == HALFCOMPLEX))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Functions must not be of HALFCOMPLEX type.");
      return;
    }
  
  if ((proj_img->type == COMPLEX) && (volume->type == REAL))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Volume must be COMPLEX if projection image is.");
      return;
    }
  
  Try {  
    if (axis == 0)
      {
        if ((volume->type == REAL) && (proj_img->type == REAL))
          xray_backprojection_sax_r_r (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
        else if ((volume->type == COMPLEX) && (proj_img->type == REAL))
          xray_backprojection_sax_r_c (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
        else /* Both COMPLEX */
          xray_backprojection_sax_c_c (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
      }
    else /* axis == 1 */
      {
        if ((volume->type == REAL) && (proj_img->type == REAL))
          xray_backprojection_say_r_r (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
        else if ((volume->type == COMPLEX) && (proj_img->type == REAL))
          xray_backprojection_say_r_c (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
        else /* Both COMPLEX */
          xray_backprojection_say_c_c (proj_img, theta_deg, tilt_axis_par_shift_px, volume, weight);
      }
  } CATCH_RETURN_VOID (_e);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_convolution (gfunc3 *gf1, int trafo_gf1, gfunc3 *gf2, int trafo_gf2)
{
  CEXCEPTION_T _e = EXC_NONE;
  float factor;

  CAPTURE_NULL_VOID (gf1);
  CAPTURE_NULL_VOID (gf2);
  GFUNC_CAPTURE_UNINIT_VOID (gf1);
  GFUNC_CAPTURE_UNINIT_VOID (gf2);

  if (GFUNC_IS_2D (gf1))
    {
      if (!GFUNC_IS_2D (gf2))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Grid functions must agree in dimension.");
          return;
        }
      factor = 1.0 / (2 * M_PI);
    }
  else
    {
      if (GFUNC_IS_2D (gf2))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Grid functions must agree in dimension.");
          return;
        }
      factor = 1.0 / (2 * M_PI * M_SQRT2PI);
    }

  if (trafo_gf1 && (gf1->type == HALFCOMPLEX))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Cannot transform HALFCOMPLEX grid function 1.");
      return;
    }
  if (trafo_gf2 && (gf2->type == HALFCOMPLEX))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Cannot transform HALFCOMPLEX grid function 2.");
      return;
    }

  if (trafo_gf1)
    {
      Try { fft_forward (gf1); }  CATCH_RETURN_VOID (_e);
    }
    
  if (trafo_gf2)
    {
      Try { fft_forward (gf2);  }  CATCH_RETURN_VOID (_e);
    }

  Try { gfunc3_mul (gf1, gf2);  }  CATCH_RETURN_VOID (_e);

  if (trafo_gf1)
    {
      Try { fft_backward (gf1); }  CATCH_RETURN_VOID (_e);
    }

  if (trafo_gf2)
    {
      Try { fft_backward (gf2); }  CATCH_RETURN_VOID (_e);
    }

  gfunc3_scale (gf1, factor);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
image_rotation_r (gfunc3 *proj_img, float const psi_deg)
{
  CEXCEPTION_T _e = EXC_NONE;
  int ix, iy;
  int Nx_new, Ny_new;
  size_t idx = 0, ntotal_new;
  float cos_psi, sin_psi;
  float dx_rot[2], dy_rot[2];
  float xlen, ylen, xlen_new, ylen_new, xmin_new, ymin_new;
  float vx0[2];
  vec3 v = {0.0, 0.0, 0.0}, img_x0;
  float *fvals_rot = NULL;
  
  CAPTURE_NULL_VOID (proj_img);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);

  if (!GFUNC_IS_2D(proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 2d functions supported.");
      return;
    }
  if (fabsf (psi_deg) > 180.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Rotation angle must be between -180 and +180 degrees");
      return;
    }

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
  Try { fvals_rot = (float *) ali16_malloc (ntotal_new * sizeof (float)); }  CATCH_RETURN_VOID (_e);
  for (idx = 0; idx < ntotal_new; idx++) fvals_rot[idx] = 0.0;


  /* TODO: inline the interpolation */
  for (iy = 0, idx = 0; iy < Ny_new; iy++)
    {
      vx0[0] = v[0];
      vx0[1] = v[1];
      
      for (ix = 0; ix < Nx_new; ix++, idx++)
        {
          fvals_rot[idx] = gfunc3_interp_linear_2d_r (proj_img, v);
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
image_rotation_c (gfunc3 *proj_img, float const psi_deg)
{
  CEXCEPTION_T _e = EXC_NONE;
  int ix, iy;
  int Nx_new, Ny_new;
  size_t idx = 0, ntotal_new;
  float cos_psi, sin_psi;
  float dx_rot[2], dy_rot[2];
  float xlen, ylen, xlen_new, ylen_new, xmin_new, ymin_new;
  float vx0[2];
  vec3 v = {0.0, 0.0, 0.0}, img_x0;
  float complex *fvals_rot = NULL;
  

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
  Try { 
    fvals_rot = (float complex *) ali16_malloc (ntotal_new * sizeof (float complex)); 
  }  CATCH_RETURN_VOID (_e);
  for (idx = 0; idx < ntotal_new; idx++) fvals_rot[idx] = 0.0;


  /* TODO: inline the interpolation */
  for (iy = 0, idx = 0; iy < Ny_new; iy++)
    {
      vx0[0] = v[0];
      vx0[1] = v[1];
      
      for (ix = 0; ix < Nx_new; ix++, idx++)
        {
          fvals_rot[idx] = gfunc3_interp_linear_2d_c (proj_img, v);
          v[0] += dx_rot[0];
          v[1] += dx_rot[1];
        }
        
      v[0] = vx0[0] + dy_rot[0];
      v[1] = vx0[1] + dy_rot[1];
    }
    
  free (proj_img->fvals);
  proj_img->fvals     = (float *) fvals_rot;
  proj_img->shape[0]  = Nx_new;
  proj_img->shape[1]  = Ny_new;
  proj_img->ntotal    = ntotal_new;
  vec3_copy (proj_img->x0, img_x0);
  gfunc3_compute_xmin_xmax (proj_img);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
image_rotation (gfunc3 *proj_img, float const psi_deg)
{
  CAPTURE_NULL_VOID (proj_img);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);

  if (!GFUNC_IS_2D(proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 2d functions supported.");
      return;
    }
  if (proj_img->type == HALFCOMPLEX)  
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "HALFCOMPLEX functions not supported.");
      return;
    }
  if (fabsf (psi_deg) > 180.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Rotation angle must be between -180 and +180 degrees");
      return;
    }

  if (psi_deg == 0.0)
    return;

  PRINT_VERBOSE ("psi_deg = %f\n", psi_deg);

  if (proj_img->type == REAL)
    image_rotation_r (proj_img, psi_deg);
  else 
    image_rotation_c (proj_img, psi_deg);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
histogram_normalization (gfunc3 *proj_img, idx3 bg_ix0, idx3 const bg_shp)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i;
  int static count = 0;
  size_t j, *idcs = NULL;
  float bg_avg = 0.0, bg_var = 0.0;
  vec3 bg_x0;
  gfunc3 *bg_patch;

  CAPTURE_NULL_VOID (proj_img);
  CAPTURE_NULL_VOID (bg_ix0);
  CAPTURE_NULL_VOID (bg_shp);
  GFUNC_CAPTURE_UNINIT_VOID (proj_img);
  
  /* TODO: implement complex version */
  if (!GFUNC_IS_REAL (proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
      return;
    }
  if (!GFUNC_IS_2D(proj_img))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 2d functions supported.");
      return;
    }

  count++;
  
  Try { bg_patch = new_gfunc3 (); }  CATCH_RETURN_VOID (_e);
    
  /* TODO: Use 4 corner patches if bg_ix0[2] == -1 (see options) */
  if (bg_ix0[2] == -1)
    bg_ix0[2] = 0;
  
  /* Init bg_patch grid */
  for (i = 0; i < 3; i++)
    bg_x0[i] = proj_img->xmin[i] + (bg_ix0[i] + bg_shp[i] / 2 ) * proj_img->csize[i];

  gfunc3_init (bg_patch, bg_x0, proj_img->csize, bg_shp, REAL);

  /* Extract values from proj_img */
  Try { idcs = gfunc3_subgrid_flatidcs (proj_img, bg_patch); }  CATCH_RETURN_VOID (_e);
  
  for (j = 0; j < bg_patch->ntotal; j++)
    bg_patch->fvals[j] = proj_img->fvals[idcs[j]];
   
  free (idcs);
  
  if (DEBUGGING)
    {
      Try { temp_mrc_out (bg_patch, "bg_patch_", count); }  CATCH_RETURN_VOID (_e);
    }
  
  /* Compute statistics */
  bg_avg = gfunc3_mean (bg_patch);
  bg_var = gfunc3_variance (bg_patch, &bg_avg);
  
  PRINT_VERBOSE ("Background mean    : %f\n", bg_avg);
  PRINT_VERBOSE ("Background variance: %f\n", bg_var);
      
  /* Normalize proj_img -> (proj_img - bg_avg) / sqrt(bg_var) */
  gfunc3_add_constant (proj_img, -bg_avg);
  gfunc3_scale (proj_img, 1.0 / sqrtf (bg_var));
  
  gfunc3_free (&bg_patch);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
probability_normalization (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;

  float integral = 0.0;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
 
  /* TODO: implement complex version */
  if (!GFUNC_IS_REAL (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
      return;
    }

  Try { integral = lp_integral (gf, TRAPEZOIDAL);  }  CATCH_RETURN_VOID (_e);
 
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
  
  CAPTURE_NULL (gf, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);
 
  /* TODO: implement complex version */
  if (!GFUNC_IS_REAL (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
      return FLT_MAX;
    }

  switch (rule)
    {
      case TRAPEZOIDAL:
        for (i = 0; i < gf->ntotal; i++)
          integral += gf->fvals[i];
        
        integral *= vec3_product (gf->csize);
        
      default:
        EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Unknown integration rule.");
        return FLT_MAX;
    }
    
  return integral;
}

/*-------------------------------------------------------------------------------------------------*/

