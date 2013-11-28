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
#include <nfft3.h>

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

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf)
{
  int i;

  CAPTURE_NULL_VOID (gf);
  if (gf->is_initialized)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFINIT, "Grid function must not be initialized!");

  /* Store old x0 internally and set new one to zero */
  vec3_copy (gf->_fbuf, gf->x0);
  vec3_set_all (gf->x0, 0.0);

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf->csize[i] = 2 * M_PI / (gf->shape[i] * gf->csize[i]);

  /* Store old shape[0]; New shape[0]: only half of the frequencies are stored (half-complex) */
  gf->_ntmp = gf->shape[0];
  gf->shape[0] = gf->shape[0] / 2 + 1;

  gf->is_halfcomplex = TRUE;

  gfunc3_compute_xmin_xmax (gf);

  gf->ntotal = idx3_product (gf->shape);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_hc_grid_bwd_reciprocal (gfunc3 *gf_hc)
{
  int i;
  
  CAPTURE_NULL_VOID (gf_hc);
  if (gf_hc->is_initialized)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFINIT, "Grid function must not be initialized!");

  /* Restore old x0 and shape[0] */
  vec3_copy (gf_hc->x0, gf_hc->_fbuf);
  gf_hc->shape[0] = gf_hc->_ntmp;
  gf_hc->_ntmp = 0;

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf_hc->csize[i] = 2 * M_PI / (gf_hc->shape[i] * gf_hc->csize[i]);

  gf_hc->is_halfcomplex = FALSE;
 
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
  CEXCEPTION_T _e = EXC_NONE;
  int nx;
  size_t n_ft;
  idx3 padding = {0, 0, 0};
  fftwf_complex *d_ft = NULL;
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (gf->is_halfcomplex)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Expected REAL type function.");
      return;
    }

  /* Zero-padding first */

  if (fft_padding > 0)
    {
      idx3_set_all (padding, fft_padding);
      
      if (GFUNC_IS_2D (gf))
        padding[2] = 0;
      
      Try { gfunc3_zeropad (gf, padding); } CATCH_RETURN_VOID (_e);
    }

  nx = gf->shape[0] / 2 + 1;
  n_ft = nx * gf->shape[1] * gf->shape[2];
  Try { d_ft = (fftwf_complex *) ali16_malloc (n_ft * sizeof(fftw_complex)); } CATCH_RETURN_VOID (_e);

  fft_fwd_premod (gf);

  if (GFUNC_IS_2D (gf))
    p = fftwf_plan_dft_r2c_2d (gf->shape[1], gf->shape[0], gf->fvals, d_ft, FFTW_ESTIMATE);

  else
    p = fftwf_plan_dft_r2c_3d (gf->shape[2], gf->shape[1], gf->shape[0], 
      gf->fvals, d_ft, FFTW_ESTIMATE);
    
  fftwf_execute (p);
  fftwf_destroy_plan (p);

  fftwf_free (gf->fvals);
  gf->is_initialized = FALSE;
  gfunc3_grid_fwd_reciprocal (gf);

  gf->fvals = (float *) d_ft;
  gf->is_initialized = TRUE;

  fft_fwd_postmod (gf);
  
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
fft_backward (gfunc3 *gf_hc)
{
  CEXCEPTION_T _e = EXC_NONE;
  size_t n_ift;
  float *d_ift = NULL;
  idx3 padding = {0, 0, 0};
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf_hc);
  GFUNC_CAPTURE_UNINIT_VOID (gf_hc);
  if (!gf_hc->is_halfcomplex)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Expected HALFCOMPLEX type function.");
      return;
    }

  n_ift = gf_hc->_ntmp * gf_hc->shape[1] * gf_hc->shape[2];
  Try { d_ift = (float *) ali16_malloc (n_ift * sizeof (float)); } CATCH_RETURN_VOID (_e);

  fft_bwd_premod (gf_hc);

  if (GFUNC_IS_2D (gf_hc))
    p = fftwf_plan_dft_c2r_2d (gf_hc->shape[1], gf_hc->_ntmp,
      (fftwf_complex *) gf_hc->fvals, d_ift, FFTW_ESTIMATE);

  else
    p = fftwf_plan_dft_c2r_3d (gf_hc->shape[2], gf_hc->shape[1],
      gf_hc->_ntmp, (fftwf_complex *) gf_hc->fvals, d_ift, FFTW_ESTIMATE);
    
  fftwf_execute (p);
  fftwf_destroy_plan (p);

  fftwf_free (gf_hc->fvals);
  gf_hc->is_initialized = FALSE;
  gfunc3_hc_grid_bwd_reciprocal (gf_hc);
  
  gf_hc->fvals = d_ift;
  gf_hc->is_initialized = TRUE;

  if (fft_padding > 0)
    {
      idx3_set_all (padding, fft_padding);
      
      if (GFUNC_IS_2D (gf_hc))
        padding[2] = 0;
        
      Try { gfunc3_unpad (gf_hc, padding); } CATCH_RETURN_VOID (_e);
    }

  fft_bwd_postmod (gf_hc);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
nfft3_normalize_freqs (double *freqs, vec3 const csize, size_t nfreqs)
{
  int i;
  size_t j;
  double *cur_freq;
  
  for (j = 0, cur_freq = freqs; j < nfreqs; j++, cur_freq += 3)
    {
      for (i = 0; i < 3; i++)
        {
          cur_freq[i] *= csize[i] / (2*M_PI);
          if (fabsf (cur_freq[i]) > 0.5)
            {
              EXC_THROW_CUSTOMIZED_PRINT (EXC_COMPUTE, "Normalized frequencies must be between "
                "-0.5 and 0.5 (value %g for x[%lu][%d]).", cur_freq[i], j, i);
              return;
            }
        }
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
nfft3_postmod (float *ft_re, float *ft_im, vec3 const x0, vec3 const csize, float const *freqs, 
               size_t nfreqs)
{
  size_t j;
  float sp, co, si, re, im, cs_prod;
  float const *cur_freq = freqs;
  
  cs_prod = vec3_product (csize);
  
  // TODO: interpolation function!
  for (j = 0; j < nfreqs; j++, cur_freq += 3)
    {
      /* Multiply complex ft with exp(-<x0, freq_j>) * [product of csize] */
      sp = x0[0] * cur_freq[0] + x0[1] * cur_freq[1] + x0[2] * cur_freq[2];
      co = cosf (-sp);
      si = sinf (-sp);
      re = ft_re[j] * co - ft_im[j] * si;
      im = ft_re[j] * si + ft_im[j] * co;
      ft_re[j] = re * cs_prod;
      ft_im[j] = im;
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
nfft3_transform (gfunc3 const *f_re, gfunc3 const *f_im, float const *freqs, size_t nfreqs, 
                 float *ft_re, float *ft_im)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int i;
  size_t j;
  float *real = NULL, *imag = NULL;
  vec3 x0 = {0.0, 0.0, 0.0}, cs = {1.0, 1.0, 1.0};
  idx3 shp = {1, 1, 1};
  size_t ntotal = 0; 

  nfft_plan p;

  CAPTURE_NULL_VOID (freqs);
  
  if ( (f_re == NULL) && (f_im == NULL) )
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "F_RE and F_IM may not be NULL at the same time.");
      return;
    }
  if (f_im != NULL)
    {
      GFUNC_CAPTURE_UNINIT_VOID (f_im);
      if (GFUNC_IS_2D (f_im))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 3D functions are supported.");
          return;
        }
      vec3_copy (x0, f_im->x0);
      vec3_copy (cs, f_im->csize);
      idx3_copy (shp, f_im->shape);
      ntotal = f_im->ntotal;
    }
  if (f_re != NULL)
    {
      GFUNC_CAPTURE_UNINIT_VOID (f_re);
      if (GFUNC_IS_2D (f_re))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 3D functions are supported.");
          return;
        }
      vec3_copy (x0, f_re->x0);
      vec3_copy (cs, f_re->csize);
      idx3_copy (shp, f_re->shape);
      ntotal = f_re->ntotal;
    }  
  if ( (f_re != NULL) && (f_im != NULL) )
    {
      if (f_re->ntotal != f_im->ntotal)
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Real and imaginary parts must agree in size.");
          return;
        }
    }

  /* Use real and imag as dummies if the corresponding output array is NULL */
  if (ft_re == NULL)
    {
      Try { real = (float *) ali16_malloc (nfreqs * sizeof (float)); } CATCH_RETURN_VOID (_e);
    }
  else
    real = ft_re;

  if (ft_im == NULL)
    {
      Try { imag = (float *) ali16_malloc (nfreqs * sizeof (float)); } CATCH_RETURN_VOID (_e);
    }
  else
    imag = ft_im;

  nfft_init_3d (&p, shp[2], shp[1], shp[0], nfreqs);
  
  /* TODO: Swap axes 0 and 2 ?? */
  for (j = 0; j < ntotal; j++)
  {
    ((double *) p.f_hat)[2 * j]     = (f_re != NULL) ? (double) f_re->fvals[j] : 0.0;
    ((double *) p.f_hat)[2 * j + 1] = (f_im != NULL) ? (double) f_im->fvals[j] : 0.0;
  }
  
  /* Copy over the points due to only double precision support */
  /* TODO: Swap axes 0 and 2 ?? */
  for (j = 0; j < nfreqs; j++)
    {
      for (i = 0; i < 3; i++)
        p.x[3 * j + i] = (double) freqs[3 * j + i];
    }
  
  /* Prepare and execute the transform */
  Try { nfft3_normalize_freqs (p.x, cs, nfreqs); } CATCH_RETURN_VOID (_e);
  nfft_precompute_one_psi (&p);
  nfft_trafo (&p);
  
  for (j = 0; j < nfreqs; j++)
    {
      real[j] = (float) ((double *) p.f)[2 * j];
      imag[j] = (float) ((double *) p.f)[2 * j + 1];
    }

  /* The plan is no longer needed, so free the space */
  nfft_finalize (&p);
  
  /* Apply the modification due to shift and scaling */
  nfft3_postmod (real, imag, x0, cs, freqs, nfreqs);
  
  if (ft_re == NULL)
    free (real);

  if (ft_im == NULL)
    free (imag);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_backward_fullcomplex (gfunc3 *f_re, gfunc3 *f_im, vec3 const shift)
{
  
}

