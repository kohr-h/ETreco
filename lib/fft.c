/*
 * fft.c -- fast Fourier transform related routines
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
#include <complex.h>
#include <fftw3.h>
#include <nfft3.h>

#include "CException.h"

#include "vec3.h"
#include "misc.h"

#include "fft.h"
#include "gfunc3.h"
#include "mrc.h"
#include "vfunc.h"


#define TAYLOR_THRESHOLD 5e-02


extern int fft_padding;

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf)
{
  int i;

  CAPTURE_NULL_VOID (gf);
  if ((gf->is_initialized) && (gf->type == REAL))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFINIT, "REAL grid function must not be initialized!");
      return;
    }
  if (gf->type == HALFCOMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "HALFCOMPLEX type not applicable!");
      return;
    }

  /* Store old x0 internally and set new one to zero */
  vec3_copy (gf->_fbuf, gf->x0);
  vec3_set_all (gf->x0, 0.0);

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf->csize[i] = 2 * M_PI / (gf->shape[i] * gf->csize[i]);

  if (GFUNC_IS_REAL (gf))
    {
      /* Store old shape[0]; New shape[0]: only half of the frequencies are stored (half-complex) */
      gf->_ntmp = gf->shape[0];
      gf->shape[0] = gf->shape[0] / 2 + 1;
      gf->type = HALFCOMPLEX;
    }

  gfunc3_compute_xmin_xmax (gf);

  gf->ntotal = idx3_product (gf->shape);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_grid_bwd_reciprocal (gfunc3 *gf)
{
  int i;
  
  CAPTURE_NULL_VOID (gf);
  if (gf->type == REAL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "REAL type not applicable!");
      return;
    }

  /* Restore old x0 */
  vec3_copy (gf->x0, gf->_fbuf);
  
  if (gf->type == HALFCOMPLEX)
    {
      /* Function is the transform of a real-valued function. First dimension (shape[0]) was 
       * halved in forward transform. Restore it. 
       */
      gf->shape[0] = gf->_ntmp;
      gf->_ntmp = 0;
      gf->type = REAL;
    }

  /* New csize = 2*PI / (old shape * old csize) */
  for (i = 0; i < 3; i++)
    gf->csize[i] = 2 * M_PI / (gf->shape[i] * gf->csize[i]);

 
  /* Restore xmin and xmax */
  gfunc3_compute_xmin_xmax (gf);

  gf->ntotal = idx3_product (gf->shape);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

float fft_fwd_interp_kernel (float t)
{
  if (fabsf (t) > TAYLOR_THRESHOLD)
    return 2.0 * (1.0 - cosf (t)) / (t * t);
  
  else
    {
      float tmp, ker_val = 1.0;
      
      tmp      = t * t / 12.0;
      ker_val -= tmp;
      tmp     *= t * t / 30.0;
      ker_val += tmp;
      tmp     *= t * t / 56.0;
      ker_val -= tmp;
      
      return ker_val;
    }
}

/*-------------------------------------------------------------------------------------------------*/

float fft_bwd_interp_kernel (float t)
{
  if (fabsf (t) > TAYLOR_THRESHOLD)
    return t * t / (2.0 * (1.0 - cosf (t)));
  
  else
    {
      float tmp, ker_val = 1.0;
      
      tmp      = t * t / 12.0;
      ker_val += tmp;
      tmp     *= t * t / 20.0;
      ker_val += tmp;
      tmp     *= t * t / 25.2;
      ker_val += tmp;
      
      return ker_val;
    }
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_fwd_premod_r (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  float cs_fac, *pfval = gf->fvals;

  cs_fac = (GFUNC_IS_2D (gf)) ? 1.0 : gf->csize[2];
  cs_fac *= gf->csize[1] * gf->csize[0];

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              *pfval *= cs_fac;

              /* Multiply by (-1)^(ix+ix+iz) */
              if (isum % 2)
                *pfval = -(*pfval);
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_fwd_premod_c (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  float cs_fac;
  float complex *pfval = (float complex *) gf->fvals;

  cs_fac = (GFUNC_IS_2D (gf)) ? 1.0 : gf->csize[2];
  cs_fac *= gf->csize[1] * gf->csize[0];

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              *pfval *= cs_fac;

              /* Multiply by (-1)^(ix+ix+iz) */
              if (isum % 2)
                *pfval = -(*pfval);
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
fft_fwd_premod (gfunc3 *gf)
{
  if (gf->type == REAL)
    fft_fwd_premod_r (gf);
  else if (gf->type == COMPLEX)
    fft_fwd_premod_c (gf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_fwd_postmod (gfunc3 *gf)
{
  int ix, iy, iz, isum, Nx, Nsum;
  float spx, spy, spz, sp, ker_fac_z, ker_fac_y, ker_fac;
  float complex Nfac = 1.0, *pfval = (float complex *) gf->fvals;

  Nx = (gf->type == HALFCOMPLEX) ? gf->_ntmp : gf->shape[0];

  /* Compute (-i)^(Nx + Ny[ + Nz]) */
  if (GFUNC_IS_2D (gf))  Nsum = Nx + gf->shape[1];
  else                   Nsum = Nx + gf->shape[1] + gf->shape[2];

  switch (Nsum % 4)
    {
    case 0: Nfac =  1.0;     break;
    case 1: Nfac = -1.0 * I; break;
    case 2: Nfac = -1.0;     break;
    case 3: Nfac =  1.0 * I; break;
    }

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      ker_fac_z = (GFUNC_IS_2D (gf)) ? 
        1.0 : 
        fft_fwd_interp_kernel (-M_PI + 2 * M_PI * iz / gf->shape[2]) / M_SQRT2PI;
      spz = gf->_fbuf[2] * (gf->xmin[2] + iz * gf->csize[2]);

      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          ker_fac_y = fft_fwd_interp_kernel (-M_PI + 2.0 * M_PI * iy / gf->shape[1]) / M_SQRT2PI;
          spy = gf->_fbuf[1] * (gf->xmin[1] + iy * gf->csize[1]);
          isum = iz + iy;

          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              /* Multiply by (-1)^(Nx+Ny[+Nz]) * (-1)^(ix+iy+iz) * exp(-i<x0,xi_j>) * kernel_fac */
              ker_fac = fft_fwd_interp_kernel (-M_PI + 2.0 * M_PI * ix / Nx) / M_SQRT2PI 
                * ker_fac_y * ker_fac_z;
              spx = gf->_fbuf[0] * (gf->xmin[0] + ix * gf->csize[0]);
              sp = spx + spy + spz;
              
              *pfval *=  cexpf (-sp * I);
              if (isum % 2) *pfval = -(*pfval);
              *pfval *= ker_fac * Nfac;
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
  size_t n_ft;
  idx3 padding = {0, 0, 0};
  fftwf_complex *d_ft = NULL;
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (gf->type == HALFCOMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "HALFCOMPLEX type not applicable!");
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

  if (DEBUGGING)
    temp_mrc_out (gf, "zero-padded", 0);

  fft_fwd_premod (gf);

  if (DEBUGGING)
    temp_mrc_out (gf, "fwd_premod", 0);

  if (GFUNC_IS_REAL (gf))
    {
      /* Need to allocate a new array which may be slightly larger than the original one. Store 
       * the transform there using half-complex FFT.
       */
      n_ft = (gf->shape[0] / 2 + 1) * gf->shape[1] * gf->shape[2];
      Try { 
        d_ft = (fftwf_complex *) ali16_malloc (n_ft * sizeof(fftw_complex));
      } CATCH_RETURN_VOID (_e);

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
      if (DEBUGGING)
        temp_mrc_out (gf, "fwd_transformed", 0);
    }
  else if (GFUNC_IS_COMPLEX (gf))
    {
      /* We can use in-place transforms */
      if (GFUNC_IS_2D (gf))
        p = fftwf_plan_dft_2d (gf->shape[1], gf->shape[0], 
          (fftwf_complex *) gf->fvals, (fftwf_complex *) gf->fvals, 
          FFTW_FORWARD, FFTW_ESTIMATE);

      else
        p = fftwf_plan_dft_3d (gf->shape[2], gf->shape[1], gf->shape[0], 
          (fftwf_complex *) gf->fvals, (fftwf_complex *) gf->fvals, 
          FFTW_FORWARD, FFTW_ESTIMATE);
    
      fftwf_execute (p);
      fftwf_destroy_plan (p);

      gfunc3_grid_fwd_reciprocal (gf);
    }

  fft_fwd_postmod (gf);
  if (DEBUGGING)
    temp_mrc_out (gf, "fwd_postmod", 0);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_bwd_premod (gfunc3 *gf)
{
  int ix, iy, iz, isum, Nx, Nsum;
  float spx, spy, spz, sp, ker_fac_z, ker_fac_y, ker_fac;
  float complex Nfac = 1.0, *pfval = (float complex *) gf->fvals;

  Nx = (gf->type == HALFCOMPLEX) ? gf->_ntmp : gf->shape[0];

  /* Compute (i)^(Nx + Ny[ + Nz]) */
  if (GFUNC_IS_2D (gf))    Nsum = Nx + gf->shape[1];
  else                     Nsum = Nx + gf->shape[1] + gf->shape[2];

  switch (Nsum % 4)
    {
    case 0: Nfac =  1.0;     break;
    case 1: Nfac =  1.0 * I; break;
    case 2: Nfac = -1.0;     break;
    case 3: Nfac = -1.0 * I; break;
    }

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      ker_fac_z = (GFUNC_IS_2D (gf)) ? 
        1.0 : 
        fft_bwd_interp_kernel (-M_PI + 2 * M_PI * iz / gf->shape[2]) * M_SQRT2PI;
      spz = gf->_fbuf[2] * (gf->xmin[2] + iz * gf->csize[2]);

      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          ker_fac_y = fft_bwd_interp_kernel (-M_PI + 2 * M_PI * iy / gf->shape[1]) * M_SQRT2PI;
          spy = gf->_fbuf[1] * (gf->xmin[1] + iy * gf->csize[1]);
          isum = iz + iy;

          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              /* Multiply by (-1)^(Nx+Ny[+Nz]) * (-1)^(ix+iy+iz) * exp(i<x0,xi_j>) * interp_kernel */
              ker_fac = ker_fac_z * ker_fac_y * 
                fft_bwd_interp_kernel (-M_PI + 2 * M_PI * ix / Nx) * M_SQRT2PI;
              spx = gf->_fbuf[0] * (gf->xmin[0] + ix * gf->csize[0]);
              sp = spx + spy + spz;
              *pfval *= cexpf (sp * I);
              if (isum % 2) *pfval = -(*pfval);
              *pfval *= ker_fac * Nfac;
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_bwd_postmod_r (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  float *pfval = gf->fvals, cs_fac;

  cs_fac = (GFUNC_IS_2D (gf)) ? 1.0 : gf->csize[2] * gf->shape[2];
  cs_fac *= gf->csize[1] * gf->shape[1] * gf->csize[0] * gf->shape[0];

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              *pfval /= cs_fac;
              
              if (isum % 2)
                *pfval = -(*pfval);
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_bwd_postmod_c (gfunc3 *gf)
{
  int isum, ix, iy, iz;
  float cs_fac;
  float complex *pfval = (float complex *) gf->fvals;

  cs_fac = (GFUNC_IS_2D (gf)) ? 1.0 : gf->csize[2] * gf->shape[2];
  cs_fac *= gf->csize[1] * gf->shape[1] * gf->csize[0] * gf->shape[0];

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          isum = iz + iy;
          for (ix = 0; ix < gf->shape[0]; ix++, isum++, pfval++)
            {
              *pfval /= cs_fac;
              
              if (isum % 2)
                *pfval = -(*pfval);
            }
        }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
fft_bwd_postmod (gfunc3 *gf)
{
  if (gf->type == REAL)
    fft_bwd_postmod_r (gf);
  else if (gf->type == COMPLEX)
    fft_bwd_postmod_c (gf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
fft_backward (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  size_t n_ift;
  float *d_ift = NULL;
  idx3 padding = {0, 0, 0};
  fftwf_plan p;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (gf->type == REAL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "REAL type not applicable!");
      return;
    }

  fft_bwd_premod (gf);

  if (DEBUGGING)
    temp_mrc_out (gf, "bwd_premod", 0);

  if (gf->type == HALFCOMPLEX)
    {
      /* Allocate an array with the size of the (original) real-valued function. Store the 
       * inverse half-complex FFT there.
       */
      n_ift = gf->_ntmp * gf->shape[1] * gf->shape[2];
      Try { d_ift = (float *) ali16_malloc (n_ift * sizeof (float)); } CATCH_RETURN_VOID (_e);

      if (GFUNC_IS_2D (gf))
        p = fftwf_plan_dft_c2r_2d (gf->shape[1], gf->_ntmp,
          (fftwf_complex *) gf->fvals, d_ift, FFTW_ESTIMATE);

      else
        p = fftwf_plan_dft_c2r_3d (gf->shape[2], gf->shape[1],
          gf->_ntmp, (fftwf_complex *) gf->fvals, d_ift, FFTW_ESTIMATE);
        
      fftwf_execute (p);
      fftwf_destroy_plan (p);

      fftwf_free (gf->fvals);
      gf->is_initialized = FALSE;
      gfunc3_grid_bwd_reciprocal (gf);
      
      gf->fvals = d_ift;
      gf->is_initialized = TRUE;
    }
  else if (gf->type == COMPLEX)
    {
      /* Perform the transform in-place */
      if (GFUNC_IS_2D (gf))
        p = fftwf_plan_dft_2d (gf->shape[1], gf->shape[0],
          (fftwf_complex *) gf->fvals, (fftwf_complex *) gf->fvals, FFTW_BACKWARD, FFTW_ESTIMATE);

      else
        p = fftwf_plan_dft_3d (gf->shape[2], gf->shape[1], gf->shape[0],
          (fftwf_complex *) gf->fvals, (fftwf_complex *) gf->fvals, FFTW_BACKWARD, FFTW_ESTIMATE);
        
      fftwf_execute (p);
      fftwf_destroy_plan (p);

      gf->is_initialized = FALSE;
      gfunc3_grid_bwd_reciprocal (gf);
      gf->is_initialized = TRUE;
    }

  if (fft_padding > 0)
    {
      idx3_set_all (padding, fft_padding);
      
      if (GFUNC_IS_2D (gf))
        padding[2] = 0;
        
      Try { gfunc3_unpad (gf, padding); } CATCH_RETURN_VOID (_e);
    }

  fft_bwd_postmod (gf);

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
nfft3_postmod (float complex *ftvals, vec3 const x0, vec3 const csize, 
               float const *freqs, size_t nfreqs)
{
  size_t j;
  float sp, cs_fac, ker_fac;
  float const *cur_freq = freqs;
  float complex *pfval = ftvals;
  
  cs_fac = vec3_product (csize);
  
  for (j = 0; j < nfreqs; j++, cur_freq += 3, pfval++)
    {
      /* Multiply complex ft with exp(-i*<x0, freq_j>) * [product of csize] * interp_kernel */
      ker_fac = fft_fwd_interp_kernel (csize[0] * cur_freq[0]) *
                fft_fwd_interp_kernel (csize[1] * cur_freq[1]) *
                fft_fwd_interp_kernel (csize[2] * cur_freq[2]) /
                (2 * M_PI * M_SQRT2PI);
      
      sp = x0[0] * cur_freq[0] + x0[1] * cur_freq[1] + x0[2] * cur_freq[2];
      *pfval *= cs_fac * ker_fac * cexpf (- I * sp);
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
nfft3_transform (gfunc3 const *gf, float const *freqs, size_t nfreqs, float complex *ftvals)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  int ix, iy, iz;
  size_t j, idx, incx;
  float const *cur_freq;
  double *cur_x;
  float complex *pfval;
  double complex *pfhat_val;

  nfft_plan p;

  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (freqs);
  CAPTURE_NULL_VOID (ftvals);
  
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  if (GFUNC_IS_2D (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only 3D grid functions are supported.");
      return;
    }
  if (gf->type != COMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Only COMPLEX grid functions are supported.");
      return;
    }

  nfft_init_3d (&p, gf->shape[0], gf->shape[1], gf->shape[2], nfreqs);
  
  /* Copy and cast the values to the double complex array in the plan. Swap x and z axes. */
  pfval = (float complex *) gf->fvals;
  pfhat_val = p.f_hat;
  idx = 0;
  incx = gf->shape[1] * gf->shape[0];
  for (iz = 0; iz < gf->shape[0]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          idx = iy * gf->shape[0] + iz;
          for (ix = 0; ix < gf->shape[2]; ix++, idx += incx)
            *(pfhat_val++) = (double complex) pfval[idx];
        }
    }

  /* Copy and cast the normalized frequencies to the double array in the plan (double precision 
   * support only).
   */
  cur_freq = freqs;
  cur_x = p.x; 
  for (j = 0; j < 3 * nfreqs; j++)
    *(cur_x++) = (double) (*(cur_freq++));

  Try { nfft3_normalize_freqs (p.x, gf->csize, nfreqs); } CATCH_RETURN_VOID (_e);
  
  /* Prepare and execute the transform */
  nfft_precompute_one_psi (&p);
  nfft_trafo (&p);
  
  /* Copy and cast the transformed values to the output array */
  pfval = ftvals;
  pfhat_val = p.f;
  for (j = 0; j < nfreqs; j++, pfval++, pfhat_val++)
    *pfval = (float complex) (*pfhat_val);

  /* The plan is no longer needed, so free the space */
  nfft_finalize (&p);
  
  /* Apply the modification due to shift and scaling (use unnormalized frequencies) */
  nfft3_postmod (ftvals, gf->x0, gf->csize, freqs, nfreqs);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
