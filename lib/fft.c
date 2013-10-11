/*
 * fft.c -- fast Fourier transform related routines
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
#include <fftw3.h>
// #include <nfft3.h>  // This is for the future

#include "CException.h"

#include "matvec3.h"
#include "misc.h"

#include "fft.h"
#include "gfunc3.h"
#include "vfunc.h"


#define FFTW_DATA_THRESHOLD 16777216
#define TAYLOR_THRESHOLD 5e-02


extern int fft_padding;

// TODO: use interpolation functions
// TODO: write NFFT functions

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf)
{
  int i;

  CAPTURE_NULL_VOID (gf);

  /* Store old x0 internally and set new one to zero */
  vec3_copy (gf->_fbuf, gf->x0);
  vec3_set_all (gf->x0, 0.0);

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf->csize[i] = 2 * M_PI / (gf->shape[i] * gf->csize[i]);

  /* Store old shape[0]; New shape[0]: only half of the frequencies are stored (half-complex) */
  gf->_ntmp = gf->shape[0];
  gf->shape[0] = gf->shape[0] / 2 + 1;

  gfunc3_compute_xmin_xmax (gf);

  gf->ntotal = idx3_product (gf->shape);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_grid_bwd_reciprocal (gfunc3 *gf_hc)
{
  int i;
  
  CAPTURE_NULL_VOID (gf_hc);

  /* Restore old x0 and shape[0] */
  vec3_copy (gf_hc->x0, gf_hc->_fbuf);
  gf_hc->shape[0] = gf_hc->_ntmp;
  gf_hc->_ntmp = 0;

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf_hc->csize[i] = 2 * M_PI / (gf_hc->shape[i] * gf_hc->csize[i]);
 
  /* Restore xmin and xmax */
  gfunc3_compute_xmin_xmax (gf_hc);

  gf_hc->ntotal = idx3_product (gf_hc->shape);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_fwd_premod (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  size_t idx;
  float fac;

  fac = (GFUNC_IS_2D (gf)) ? 1.0 : gf->csize[2] / M_SQRT2PI;
  fac *= gf->csize[1] * gf->csize[0] / (2 * M_PI);

  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, idx++)
            {
              gf->fvals[idx] *= fac;

              /* Multiply by (-1)^(ix+ix+iz) */
              if (isum % 2)
                gf->fvals[idx] = -gf->fvals[idx];
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_fwd_postmod (gfunc3 *gf_hc)
{
  int ix, iy, iz, isum, Nsum, Nfac_re = 0, Nfac_im = 0;
  float spx, spy, spz, sp, co, si, re, im;
  size_t idx;

  /* Compute (-i)^(Nx + Ny[ + Nz]) */
  if (GFUNC_IS_2D (gf_hc))
    Nsum = gf_hc->_ntmp + gf_hc->shape[1];

  else
    Nsum = gf_hc->_ntmp + gf_hc->shape[1] + gf_hc->shape[2];

  switch (Nsum % 4)
    {
    case 0:
      Nfac_re = 1;
      Nfac_im = 0;
      break;
    case 1:
      Nfac_re = 0;
      Nfac_im = -1;
      break;
    case 2:
      Nfac_re = -1;
      Nfac_im = 0;
      break;
    case 3:
      Nfac_re = 0;
      Nfac_im = 1;
      break;
    }

  idx = 0;
  for (iz = 0; iz < gf_hc->shape[2]; iz++)
    {
      spz = gf_hc->_fbuf[2] * (gf_hc->xmin[2] + iz * gf_hc->csize[2]);

      for (iy = 0; iy < gf_hc->shape[1]; iy++)
        {
          spy = gf_hc->_fbuf[1] * (gf_hc->xmin[1] + iy * gf_hc->csize[1]);
          isum = iz + iy;

          for (ix = 0; ix < gf_hc->shape[0]; ix++, idx += 2, isum++)
            {
              /* Multiply by (complex) exp(-i<x0,xi_j) */
              spx = gf_hc->_fbuf[0] * (gf_hc->xmin[0] + ix * gf_hc->csize[0]);

              sp = -(spx + spy + spz);
              co = cosf (sp);
              si = sinf (sp);

              re = gf_hc->fvals[idx] * co - gf_hc->fvals[idx + 1] * si;
              im = gf_hc->fvals[idx] * si + gf_hc->fvals[idx + 1] * co;

              /* Multiply by (-1)^(ix+iy+iz) */
              if (isum % 2)
                {
                  re = -re;
                  im = -im;
                }

              gf_hc->fvals[idx]     = Nfac_re * re - Nfac_im * im;
              gf_hc->fvals[idx + 1] = Nfac_re * im + Nfac_im * re;
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_forward (gfunc3 *gf)
{
  CEXCEPTION_T e = EXC_NONE;
  int nx;
  size_t n_ft;
  idx3 padding = {0, 0, 0};
  fftwf_complex *d_ft;
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Expected REAL type function.");

  /* Zero-padding first */

  Try
  {
    if (fft_padding > 0)
      {
        idx3_set_all (padding, fft_padding);
        
        if (GFUNC_IS_2D (gf))
          padding[2] = 0;
        
        gfunc3_zeropad (gf, padding);
      }

    nx = gf->shape[0] / 2 + 1;
    n_ft = nx * gf->shape[1] * gf->shape[2];
    d_ft = (fftwf_complex *) ali16_malloc (n_ft * sizeof(fftw_complex));

    fft_fwd_premod (gf);

    if (GFUNC_IS_2D (gf))
      p = fftwf_plan_dft_r2c_2d (gf->shape[1], gf->shape[0], gf->fvals, d_ft, FFTW_ESTIMATE);

    else
      p = fftwf_plan_dft_r2c_3d (gf->shape[2], gf->shape[1], gf->shape[0], 
        gf->fvals, d_ft, FFTW_ESTIMATE);
      
    fftwf_execute (p);
    fftwf_destroy_plan (p);

    fftwf_free (gf->fvals);
    gf->fvals = (float *) d_ft;

    gf->is_halfcomplex = 1;
    gfunc3_grid_fwd_reciprocal (gf);

    fft_fwd_postmod (gf);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_bwd_premod (gfunc3 *gf_hc)
{
  int ix, iy, iz, isum, Nfac_re = 0, Nfac_im = 0, Nsum;
  float fac, spx, spy, spz, sp, co, si, re, im;
  size_t idx;

  fac = (GFUNC_IS_2D (gf_hc)) ? 1.0 : gf_hc->csize[2] / M_SQRT2PI;
  fac *= gf_hc->csize[1] * gf_hc->csize[0] / (2 * M_PI);

  /* Compute (i)^(Nx + Ny[ + Nz]) */
  if (GFUNC_IS_2D (gf_hc))
    Nsum = gf_hc->_ntmp + gf_hc->shape[1];

  else
    Nsum = gf_hc->_ntmp + gf_hc->shape[1] + gf_hc->shape[2];

  switch (Nsum % 4)
    {
    case 0:
      Nfac_re = 1;
      Nfac_im = 0;
      break;
    case 1:
      Nfac_re = 0;
      Nfac_im = 1;
      break;
    case 2:
      Nfac_re = -1;
      Nfac_im = 0;
      break;
    case 3:
      Nfac_re = 0;
      Nfac_im = -1;
      break;
    }

  idx = 0;
  for (iz = 0; iz < gf_hc->shape[2]; iz++)
    {
      spz = gf_hc->_fbuf[2] * (gf_hc->xmin[2] + iz * gf_hc->csize[2]);

      for (iy = 0; iy < gf_hc->shape[1]; iy++)
        {
          spy = gf_hc->_fbuf[1] * (gf_hc->xmin[1] + iy * gf_hc->csize[1]);
          isum = iz + iy;

          for (ix = 0; ix < gf_hc->shape[0]; ix++, idx += 2, isum++)
            {
              /* Multiply by (complex) exp(i*<x0*xi_j>) */
              spx = gf_hc->_fbuf[0] * (gf_hc->xmin[0] + ix * gf_hc->csize[0]);
              sp = (spx + spy + spz);
              co = cosf (sp);
              si = sinf (sp);

              re = fac * (gf_hc->fvals[idx] * co - gf_hc->fvals[idx + 1] * si);
              im = fac * (gf_hc->fvals[idx] * si + gf_hc->fvals[idx + 1] * co);

              /* Multiply by (-1)^(ix+iy+iz) */
              if (isum % 2)
                {
                  re = -re;
                  im = -im;
                }

              gf_hc->fvals[idx] = Nfac_re * re - Nfac_im * im;
              gf_hc->fvals[idx + 1] = Nfac_re * im + Nfac_im * re;
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_bwd_postmod (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  size_t idx;

  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, idx++)
            {
              if (isum % 2)
                gf->fvals[idx] = -gf->fvals[idx];
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_backward (gfunc3 *gf)
{
  size_t n_ift;
  float *d_ift;
  idx3 padding = {0, 0, 0};
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (!gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Expected HALFCOMPLEX type function.");

  n_ift = gf->_ntmp * gf->shape[1] * gf->shape[2];
  d_ift = (float *) ali16_malloc (n_ift * sizeof (float));

  fft_bwd_premod (gf);

  if (GFUNC_IS_2D (gf))
    p = fftwf_plan_dft_c2r_2d (gf->shape[1], gf->_ntmp,
      (fftwf_complex *) gf->fvals, d_ift, FFTW_ESTIMATE);

  else
    p = fftwf_plan_dft_c2r_3d (gf->shape[2], gf->shape[1],
      gf->_ntmp, (fftwf_complex *) gf->fvals, d_ift, FFTW_ESTIMATE);
    
  fftwf_execute (p);
  fftwf_destroy_plan (p);

  fftwf_free (gf->fvals);
  gf->fvals = d_ift;

  gf->is_halfcomplex = 0;
  gfunc3_grid_bwd_reciprocal (gf);

  if (fft_padding > 0)
    {
      idx3_set_all (padding, fft_padding);
      
      if (GFUNC_IS_2D (gf))
        padding[2] = 0;
        
      gfunc3_unpad (gf, padding);
    }

  fft_bwd_postmod (gf);

  return;
}

/*-------------------------------------------------------------------------------------------------*/
