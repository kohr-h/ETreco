/*
 * gfunc3.c -- 3-dimensional grid functions
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
#include <limits.h>

#if HAVE_CBLAS == 1
#include <gsl/gsl_blas.h>
#endif

#include "CException.h"

#include "matvec3.h"
#include "misc.h"
#include "simd.h"

#include "gfunc3.h"
#include "vfunc.h"


/*-------------------------------------------------------------------------------------------------*
 * Allocation
 *-------------------------------------------------------------------------------------------------*/

gfunc3 *
new_gfunc3 (void)
{
  gfunc3 *gf;
  CEXCEPTION_T e = EXC_NONE;
  
  Try
  {
    gf = (gfunc3 *) ali16_malloc (sizeof (gfunc3));
    
    gf->is_initialized = 0;
    gf->is_halfcomplex = 0;
    gf->fvals          = NULL;
    
    return gf;
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return NULL;
}

/*-------------------------------------------------------------------------------------------------*/
// Structure initialization
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type)
{
  CEXCEPTION_T e = EXC_NONE;
  float nul[2] = {0.0f, 0.0f};
  
  CAPTURE_NULL (gf);
  CAPTURE_NULL (cs);
  CAPTURE_NULL (shp);

  /* Initialize grid */
  idx3_copy (gf->shape, shp);
  
  if (shp[2] == 0)
    gf->shape[2] = 1;
  gf->ntotal = idx3_product (gf->shape);

  if (x0 == NULL)
    vec3_set_all (gf->x0, 0.0);
  else
    vec3_copy (gf->x0, x0);

  vec3_copy (gf->csize, cs);

  if (GFUNC_IS_2D (gf))
    {
      gf->x0[2] = 0.0;
      gf->csize[2] = 1.0;
    }

  if (gf_type == HALFCOMPLEX)
    gf->is_halfcomplex = 1;

  gfunc3_compute_xmin_xmax (gf);

  Try
  {
    /* Allocate memory for values and zero array */
    gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));

    gf->is_initialized = 1;
    gfunc3_set_all (gf, nul);
  }    
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init_from_foreign_grid (gfunc3 *gf, gfunc3 const *gf_template)
{
  CEXCEPTION_T e = EXC_NONE;
  float nul[2] = {0.0f, 0.0f};
  
  CAPTURE_NULL (gf);
  CAPTURE_NULL (gf_template);

  if (gf->is_initialized)
    free (gf->fvals);

  vec3_copy (gf->x0,    gf_template->x0);
  vec3_copy (gf->csize, gf_template->csize);
  idx3_copy (gf->shape, gf_template->shape);

  vec3_copy (gf->xmin, gf_template->xmin);
  vec3_copy (gf->xmax, gf_template->xmax);
  gf->ntotal = gf_template->ntotal;
  
  vec3_copy (gf->_fbuf, gf_template->_fbuf);
  gf->_ntmp = gf_template->_ntmp;
  
  Try
  {
    if (gf_template->is_halfcomplex)
      {
        gf->fvals = (float *) ali16_malloc (2 * gf->ntotal * sizeof (float));
        gf->is_halfcomplex = 1;
      }
    else
      {
        gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));
        gf->is_halfcomplex = 0;
      }

    gf->is_initialized = 1;
    gfunc3_set_all (gf, nul);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_csize (gfunc3 *gf, vec3 const cs)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (cs);

  vec3_copy (gf->csize, cs);
  gfunc3_compute_xmin_xmax (gf);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_x0 (gfunc3 *gf, vec3 const x0)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (x0);

  vec3_copy (gf->x0, x0);
  gfunc3_compute_xmin_xmax (gf);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_compute_xmin_xmax (gfunc3 *gf)
{
  CAPTURE_NULL (gf);

  int i;
  
  if (gf->is_halfcomplex)
    {
      gf->xmin[0] = gf->x0[0] - (gf->shape[0] - 1) * gf->csize[0];
      gf->xmax[0] = gf->x0[0];
    }
  else
    {
      gf->xmin[0] = gf->x0[0] -  gf->shape[0]      / 2 * gf->csize[0];
      gf->xmax[0] = gf->x0[0] + (gf->shape[0] - 1) / 2 * gf->csize[0];
    }
            
  for (i = 1; i < 3; i++)
    {
      gf->xmin[i] = gf->x0[i] -  gf->shape[i]      / 2 * gf->csize[i];
      gf->xmax[i] = gf->x0[i] + (gf->shape[i] - 1) / 2 * gf->csize[i];
    }
  
  if (GFUNC_IS_2D (gf))
    {
      gf->xmin[2]  = 0.0;
      gf->xmax[2]  = 0.0;
      gf->csize[2] = 1.0;
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_free (gfunc3 **pgf)
{
  if (pgf == NULL)
    return;
  
  if ((*pgf) == NULL)
    return;
  
  if ((*pgf)->fvals != NULL)
    free ((*pgf)->fvals);
  
  (*pgf)->is_initialized = 0;
  (*pgf)->is_halfcomplex = 0;

  free (*pgf);

  return;
}

/*-------------------------------------------------------------------------------------------------*/
// Function value initialization
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_assign_fvals_from_vfunc (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx, incr;
  float nul[2] = {0.0f, 0.0f};
  vec3 p;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (vf);
  GFUNC_CHECK_INIT_STATUS (gf);

  incr = (gf->is_halfcomplex) ? 2 : 1;

  gfunc3_set_all (gf, nul);

  vec3_copy (p, gf->xmin);

  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx +=incr)
            {
              VFUNC_EVAL (vf, &gf->fvals[idx], p);
              p[0] += gf->csize[0];
            }
          p[0]  = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1]  = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/
// Screen output
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_print_grid (gfunc3 const *gf, char const *intro_text)
{
  CAPTURE_NULL (gf);

  printf ("\n");
  if (intro_text != NULL)
    printf ("%s", intro_text);
  if (GFUNC_IS_2D (gf))
    printf ("    [ 2D ]");
  else
    printf ("    [ 3D ]");
  
  if (gf->is_halfcomplex)
    puts ("    [half-complex]\n");
  else
    puts ("\n");

  puts ("Shape:");
  idx3_print (gf->shape);
  printf ("\n");

  puts ("Cell size:");
  vec3_print (gf->csize);
  printf ("\n");
  
  puts ("Origin:");
  vec3_print (gf->x0);
  printf ("\n");

  puts ("Minimum:");
  vec3_print (gf->xmin);
  printf ("\n");
  
  puts ("Maximum:");
  vec3_print (gf->xmax);
  printf ("\n\n");

  return;
}

/*-------------------------------------------------------------------------------------------------*/
// Attributes
/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_min (gfunc3 const *gf)
{
  size_t i;
  float min;

  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Minimum not defined for complex functions.");
  
  min = gf->fvals[0];

  if (gf->ntotal == 1)
    return min;

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 mmin = _mm_set1_ps (FLT_MAX);
  __m128 const *p = (__m128 const *) gf->fvals;
  float *pmin = (float *) &mmin;
  
  for (i = 0; i < N4; i++, p++)
    mmin = _mm_min_ps (*p, mmin);

  min = fminf (pmin[0], pmin[1]);
  min = fminf (min, pmin[2]);
  min = fminf (min, pmin[3]);
  
  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    min = fminf (min, gf->fvals[i]);

  #else  /* !HAVE_SSE */
  for (i = 1; i < gf->ntotal; i++)
    {
      if (gf->fvals[i] < min)
        min = gf->fvals[i];
    }
  
  #endif
  return min;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_max (gfunc3 const *gf)
{
  size_t i;
  float max;
  
  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Minimum not defined for complex functions.");
  
  max = gf->fvals[0];

  if (gf->ntotal == 1)
    return max;

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 mmax = _mm_set1_ps (-FLT_MAX);
  __m128 const *p = (__m128 const *) gf->fvals;
  float *pmax = (float *) &mmax;
  
  for (i = 0; i < N4; i++, p++)
    mmax = _mm_max_ps (*p, mmax);

  max = fmaxf (pmax[0], pmax[1]);
  max = fmaxf (max, pmax[2]);
  max = fmaxf (max, pmax[3]);
  
  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    max = fmaxf (max, gf->fvals[i]);

  #else  /* !HAVE_SSE */
  for (i = 1; i < gf->ntotal; i++)
    {
      if (gf->fvals[i] > max)
        max = gf->fvals[i];
    }

  #endif
  return max;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_mean (gfunc3 const *gf)
{
  size_t i;
  float sum = 0.0;

  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Minimum not defined for complex functions.");
  
  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 msum = _mm_setzero_ps ();
  __m128 const *p = (__m128 const *) gf->fvals;
  float *psum = (float *) &msum;
  
  for (i = 0; i < N4; i++, p++)
    msum = _mm_add_ps (*p, msum);

  sum = psum[0] + psum[1] + psum[2] + psum[3];
  
  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    sum += gf->fvals[i];

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    sum += gf->fvals[i];

  #endif
  return sum / gf->ntotal;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_variance (gfunc3 const *gf, float const *pmean)
{
  size_t i;
  float mean = 0.0, tmp;
  float sum_sq = 0.0;

  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Variance not defined for complex functions.");
  
  if (pmean == NULL)
    mean = gfunc3_mean (gf);
  else
    mean = *pmean;
  
  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 msum_sq = _mm_setzero_ps (), sq, mmean = _mm_set1_ps (-mean);
  __m128 const *p = (__m128 const *) gf->fvals;
  float *psum_sq = (float *) &msum_sq;
  
  for (i = 0; i < N4; i++, p++)
    {
      sq      = _mm_add_ps (*p, mmean);
      sq      = _mm_mul_ps (sq, sq);
      msum_sq = _mm_add_ps (sq, msum_sq);
    }
    
  sum_sq = psum_sq[0] + psum_sq[1] + psum_sq[2] + psum_sq[3];
  
  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    {
      tmp     = gf->fvals[i] - mean;
      sum_sq += tmp * tmp;
    }

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    {
      tmp     = gf->fvals[i] - mean;
      sum_sq += tmp * tmp;
    }

  #endif
  return sum_sq / (gf->ntotal - 1);
}

/*-------------------------------------------------------------------------------------------------*/

int
gfunc3_grids_are_equal (gfunc3 const *gf1, gfunc3 const *gf2)
{
  CAPTURE_NULL (gf1);
  CAPTURE_NULL (gf2);
  GFUNC_CHECK_INIT_STATUS (gf1);
  GFUNC_CHECK_INIT_STATUS (gf2);
  
  return (idx3_eq (gf1->shape, gf2->shape)                 &&
          vec3_about_eq (gf1->csize, gf2->csize, EPS_GRID) &&
          vec3_about_eq (gf1->x0, gf2->x0, EPS_GRID));
}

/*-------------------------------------------------------------------------------------------------*/

int
gfunc3_grid_is_subgrid (gfunc3 const *gf, gfunc3 const *gf_sub)
{
  int i;
  float tmp;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (gf_sub);
  
  /* Subgrid shape must not exceed grid shape */
  if (!idx3_le (gf_sub->shape, gf->shape))
    return FALSE;
  
  /* Minimum must not be smaller, maximum not larger */
  if (!vec3_ge (gf_sub->xmin, gf->xmin) || !vec3_le (gf_sub->xmax, gf->xmax))
    return FALSE;

  /* Cell size must not be smaller */
  if (!vec3_ge (gf_sub->csize, gf->csize))
    return FALSE;

  for (i = 0; i < 3; i++)
    {
      /* Test if shifts between minima/maxima are integer multiples of cell size */
      tmp = (gf_sub->xmin[i] - gf->xmin[i]) / gf->csize[i];
      if (fabsf (tmp - roundf (tmp)) > EPS_GRID)
        return FALSE;

      tmp = (gf->xmax[i] - gf_sub->xmax[i]) / gf->csize[i];
      if (fabsf (tmp - roundf (tmp)) > EPS_GRID)
        return FALSE;

      /* Test if cell sizes of subgrid are integer multiples of grid cell size */
      tmp = (gf_sub->csize[i] - gf->csize[i]) / gf->csize[i];
      if (fabsf (tmp - roundf (tmp)) > EPS_GRID)
        return FALSE;
    }

  return TRUE;
}

/*-------------------------------------------------------------------------------------------------*/
// Operations
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_copy (gfunc3 *dest, gfunc3 const *src)
{
  CEXCEPTION_T e = EXC_NONE;
  size_t ntotal_flt;
  
  CAPTURE_NULL (dest);
  CAPTURE_NULL (src);
  GFUNC_CHECK_INIT_STATUS (src);

  if (dest->is_initialized)
    free (dest->fvals);
  
  vec3_copy (dest->x0,    src->x0);
  vec3_copy (dest->csize, src->csize);
  idx3_copy (dest->shape, src->shape);

  vec3_copy (dest->xmin, src->xmin);
  vec3_copy (dest->xmax, src->xmax);
  dest->ntotal = src->ntotal;

  if (src->is_halfcomplex)
    ntotal_flt = 2 * src->ntotal;
  else
    ntotal_flt = src->ntotal;
  
  Try 
  {
    dest->fvals = (float *) ali16_malloc (ntotal_flt * sizeof (float));
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  dest->is_initialized = 1;

  #if HAVE_CBLAS
  cblas_scopy (ntotal_flt, src->fvals, 1, dest->fvals, 1);

  #elif HAVE_SSE
  size_t i, N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
  __m128 *p1 = (__m128 *) dest->fvals;
  __m128 const *p2 = (__m128 const *) src->fvals;

  for (i = 0; i < N4; i++, p1++, p2++)
    *p1 = *p2;
  
  for (i = ntotal_flt - N4rem; i < ntotal_flt; i++)
    dest->fvals[i] = src->fvals[i];

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t i;
  for (i = 0; i < ntotal_flt; i++)
    dest->fvals[i] = src->fvals[i];

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy (float a, gfunc3 *gf1, gfunc3 const *gf2)
{
  size_t ntotal_flt;

  CAPTURE_NULL (gf1);
  CAPTURE_NULL (gf2);
  GFUNC_CHECK_INIT_STATUS (gf1);
  GFUNC_CHECK_INIT_STATUS (gf2);

  /* TODO: implement mixed axpy */
  if (gf1->is_halfcomplex != gf2->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Not implemented for mixed function types.");

  if (gf2->is_halfcomplex)
    ntotal_flt = 2 * gf2->ntotal;
  else
    ntotal_flt = gf2->ntotal;

  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_CBLAS
      cblas_sscal (ntotal_flt, a, gf1->fvals, 1);
      cblas_saxpy (ntotal_flt, 1.0, gf2->fvals, 1, gf1->fvals, 1);

      #elif HAVE_SSE
      size_t i, N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;
      __m128 ma = _mm_set1_ps (a);

      for (i = 0; i < N4; i++, p1++, p2++)
        {
          *p1 = _mm_mul_ps (*p1, ma);
          *p1 = _mm_add_ps (*p1, *p2);
        }
      
      for (i = ntotal_flt - N4rem; i < ntotal_flt; i++)
        gf1->fvals[i] += a * gf2->fvals[i];

      #else  /* !(HAVE_CBLAS || HAVE_SSE) */
      size_t i;
      for (i = 0; i < ntotal_flt; i++)
        gf1->fvals[i] = a * gf1->fvals[i] + gf2->fvals[i];
        
      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs;
      CEXCEPTION_T e = EXC_NONE;
      
      Try 
      {
        idcs = gfunc3_subgrid_flatidcs (gf1, gf2);

        #if HAVE_SSE
        size_t i, l = 0, N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
        __m128 const *p2 = (__m128 const *) gf2->fvals;
        __m128 ma = _mm_set1_ps (a), v1;
        float *pv1 = (float *) &v1;
  
        for (i = 0; i < N4; i++, l += 4, p2++)
          {
            v1 = _mm_set_ps (gf1->fvals[idcs[l + 3]], gf1->fvals[idcs[l + 2]], 
                  gf1->fvals[idcs[l + 1]], gf1->fvals[idcs[l]]);
            v1 = _mm_mul_ps (v1, ma);
            v1 = _mm_add_ps (v1, *p2);
  
            gf1->fvals[idcs[l]]     = pv1[0];
            gf1->fvals[idcs[l + 1]] = pv1[1];
            gf1->fvals[idcs[l + 2]] = pv1[2];
            gf1->fvals[idcs[l + 3]] = pv1[3];
          }
        
        for (i = ntotal_flt - N4rem; i < ntotal_flt; i++)
          gf1->fvals[idcs[i]] += a * gf2->fvals[i];
  
        #else  /* !HAVE_SSE */
        size_t i;
        for (i = 0; i < ntotal_flt; i++)
          gf1->fvals[idcs[i]] = a * gf1->fvals[idcs[i]] + gf2->fvals[i];
                    
        #endif
        
        free (idcs);
      }
      Catch (e)
      {
        EXC_RETHROW_REPRINT (e);
      }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy_vfunc_re (float a, gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx;
  float vfval;
  vec3 p;

  vec3_copy (p, gf->xmin);
  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              gf->fvals[idx] = a * gf->fvals[idx] + vfval;
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1] = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy_vfunc_hc (float a, gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx;
  cplx vfval;
  vec3 p;

  vec3_copy (p, gf->xmin);
  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx += 2)
            {
              VFUNC_EVAL (vf, vfval, p);
              gf->fvals[idx]     = a * gf->fvals[idx]     + vfval[0];
              gf->fvals[idx + 1] = a * gf->fvals[idx + 1] + vfval[1];
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1] = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy_vfunc (float a, gfunc3 *gf, const vfunc *vf)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (vf);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (gf->is_halfcomplex)
    gfunc3_axpy_vfunc_hc (a, gf, vf);
  else
    gfunc3_axpy_vfunc_re (a, gf, vf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_re (gfunc3 *gf1, gfunc3 const *gf2)
{
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_SSE
      size_t i, N4 = gf2->ntotal / 4, N4rem = gf2->ntotal % 4;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;

      for (i = 0; i < N4; i++, p1++, p2++)
        *p1 = _mm_mul_ps (*p1, *p2);

      for (i = gf2->ntotal - N4rem; i < gf2->ntotal; i++)
        gf1->fvals[i] *= gf2->fvals[i];
      
      #else  /* !HAVE_SSE */
      size_t i;
      
      for (i = 0; i < gf1->ntotal; i++)
        gf1->fvals[i] *= gf2->fvals[i];

      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs;
      CEXCEPTION_T e = EXC_NONE;
      
      Try 
      {
        idcs = gfunc3_subgrid_flatidcs (gf1, gf2);

        #if HAVE_SSE
        size_t i, l = 0, N4 = gf2->ntotal / 4, N4rem = gf2->ntotal % 4;
        __m128 const *p2 = (__m128 const *) gf2->fvals;
        __m128 v1;
        float *pv1 = (float *) &v1;
  
        for (i = 0; i < N4; i++, l += 4, p2++)
          {
            v1 = _mm_set_ps (gf1->fvals[idcs[l + 3]], gf1->fvals[idcs[l + 2]], 
                  gf1->fvals[idcs[l + 1]], gf1->fvals[idcs[l]]);
            v1 = _mm_mul_ps (v1, *p2);
  
            gf1->fvals[idcs[l]]     = pv1[0];
            gf1->fvals[idcs[l + 1]] = pv1[1];
            gf1->fvals[idcs[l + 2]] = pv1[2];
            gf1->fvals[idcs[l + 3]] = pv1[3];
          }
        
        for (i = gf2->ntotal - N4rem; i < gf2->ntotal; i++)
          gf1->fvals[idcs[i]] *= gf2->fvals[i];
  
        #else  /* !HAVE_SSE */
        size_t i;
        for (i = 0; i < gf2->ntotal; i++)
          gf1->fvals[idcs[i]] *= gf2->fvals[i];
  
        #endif
                
        free (idcs);
      }
      Catch (e)
      {
        EXC_RETHROW_REPRINT (e);
      }
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_hc (gfunc3 *gf1, gfunc3 const *gf2)
{
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_SSE
      size_t i, iarr;
      size_t N2 = gf2->ntotal / 2, N2rem = gf2->ntotal % 2;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;

      __m128 prod_eq, prod_across;
      float *p1f, *p2f, *pprod_eq = (float *) &prod_eq, *pprod_across = (float *) &prod_across;
      float re, im;
      
      for (i = 0; i < N2; i++, p1++, p2++)
        {
          prod_eq = _mm_mul_ps (*p1, *p2);
          p2f = (float *) p2;
          /* Switch re and im; index order would be 3,2,1,0 normally */
          prod_across = _mm_set_ps (p2f[2], p2f[3], p2f[0], p2f[1]);
          prod_across = _mm_mul_ps (prod_across, *p1);

          p1f = (float *) p1;
          p1f[0] = pprod_eq[0] - pprod_eq[1];
          p1f[1] = pprod_across[0] + pprod_across[1];
          p1f[2] = pprod_eq[2] - pprod_eq[3];
          p1f[3] = pprod_across[2] + pprod_across[3];
        }
      if (N2rem != 0)
        {
          iarr = 2 * (gf2->ntotal - N2rem);
          
          re = gf1->fvals[iarr] * gf2->fvals[iarr]     - gf1->fvals[iarr + 1] * gf2->fvals[iarr + 1];
          im = gf1->fvals[iarr] * gf2->fvals[iarr + 1] + gf1->fvals[iarr + 1] * gf2->fvals[iarr];
          gf1->fvals[iarr]     = re;
          gf1->fvals[iarr + 1] = im;
        }
        
      #else  /* !HAVE_SSE */
      size_t iarr;
      float re, im;

      for (iarr = 0; iarr < 2 * gf2->ntotal; iarr += 2)
        {
          re = gf1->fvals[iarr] * gf2->fvals[iarr]     - gf1->fvals[iarr + 1] * gf2->fvals[iarr + 1];
          im = gf1->fvals[iarr] * gf2->fvals[iarr + 1] + gf1->fvals[iarr + 1] * gf2->fvals[iarr];
          gf1->fvals[iarr]     = re;
          gf1->fvals[iarr + 1] = im;
        }

      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs;
      CEXCEPTION_T e = EXC_NONE;
      
      Try 
      {
        idcs = gfunc3_subgrid_flatidcs (gf1, gf2);

        #if HAVE_SSE
        size_t i, i2, i2arr, i1, i1next, ntotal_flt = 2 * gf2->ntotal;
        size_t N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
        __m128 const *p2 = (__m128 const *) gf2->fvals;
        __m128 v1, prod_eq, prod_across;
        float *p2f, *pprod_eq = (float *) &prod_eq, *pprod_across = (float *) &prod_across;
        float re, im;
        
        for (i = 0, i2 = 0; i < N4; i2 += 2, p2++)
          {
            i1     = 2 * idcs[i2];
            i1next = 2 * idcs[i2 + 1];
            v1 = _mm_set_ps (gf1->fvals[i1next + 1], gf1->fvals[i1next], 
              gf1->fvals[i1 + 1], gf1->fvals[i1]);
            prod_eq = _mm_mul_ps (v1, *p2);
            p2f = (float *) p2;
            /* Switch re and im; index order would be 3,2,1,0 normally */
            prod_across = _mm_set_ps (p2f[2], p2f[3], p2f[0], p2f[1]);
            prod_across = _mm_mul_ps (prod_across, v1);
  
            gf1->fvals[i1]         = pprod_eq[0] - pprod_eq[1];
            gf1->fvals[i1 + 1]     = pprod_across[0] + pprod_across[1];
            gf1->fvals[i1next]     = pprod_eq[2] - pprod_eq[3];
            gf1->fvals[i1next + 1] = pprod_across[2] + pprod_across[3];
          }
        if (N4rem != 0)
          {
            i2    = gf2->ntotal - 1;
            i1    = 2 * idcs[i2];
            i2arr = 2 * i2;
            
            re = gf1->fvals[i1] * gf2->fvals[i2arr]     - gf1->fvals[i1 + 1] * gf2->fvals[i2arr + 1];
            im = gf1->fvals[i1] * gf2->fvals[i2arr + 1] + gf1->fvals[i1 + 1] * gf2->fvals[i2arr];
            gf1->fvals[i1]     = re;
            gf1->fvals[i1 + 1] = im;
          }
  
        #else  /* !HAVE_SSE */
        size_t i1, i2, i2arr;
        
        for (i2 = 0; i2 < gf2->ntotal; i2++, i2arr += 2)
          {
            i1 = 2 * idcs[i2];
            
            re = gf1->fvals[i1] * gf2->fvals[i2arr]     - gf1->fvals[i1 + 1] * gf2->fvals[i2arr + 1];
            im = gf1->fvals[i1] * gf2->fvals[i2arr + 1] + gf1->fvals[i1 + 1] * gf2->fvals[i2arr];
            gf1->fvals[i1]     = re;
            gf1->fvals[i1 + 1] = im;
          }
          
        #endif
                
        free (idcs);
      }
      Catch (e)
      {
        EXC_RETHROW_REPRINT (e);
      }
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_hc_re (gfunc3 *gf1, gfunc3 const *gf2)
{
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_SSE
      size_t i, i2, iarr;
      size_t N2 = gf2->ntotal / 2, N2rem = gf2->ntotal % 2;

      __m128 *p1 = (__m128 *) gf1->fvals;
      float const *p2 = gf2->fvals;
      __m128 v2;

      for (i = 0; i < N2; i++, i2 += 2, p1++, p2 += 2)
        {
          v2 = _mm_set_ps (p2[1], p2[1], p2[0], p2[0]);
          *p1 = _mm_mul_ps (*p1, v2);
        }
      if (N2rem != 0)
        {
          iarr = 2 * (gf2->ntotal - N2rem);
          p2--;
          
          gf1->fvals[iarr]     *= p2[0];
          gf1->fvals[iarr + 1] *= p2[1];
        }
        
      #else  /* !HAVE_SSE */
      size_t iarr;
      float *p2 = gf2->fvals;

      for (iarr = 0; iarr < 2 * gf2->ntotal; iarr += 2, p2++)
        {
          gf1->fvals[iarr]     *= *p2;
          gf1->fvals[iarr + 1] *= *p2;
        }

      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs;
      CEXCEPTION_T e = EXC_NONE;
      
      Try 
      {
        idcs = gfunc3_subgrid_flatidcs (gf1, gf2);

        #if HAVE_SSE
        size_t i, i2, i1, i1next;
        size_t N2 = gf2->ntotal / 2, N2rem = gf2->ntotal % 2;
        float const *p2 = gf2->fvals;
        __m128 v1, v2, prod;
        float *pprod = (float *) &prod;
        
        for (i = 0, i2 = 0; i < N2; i2 += 2, p2 += 2)
          {
            i1     = 2 * idcs[i2];
            i1next = 2 * idcs[i2 + 1];
            v1 = _mm_set_ps (gf1->fvals[i1next + 1], gf1->fvals[i1next], 
              gf1->fvals[i1 + 1], gf1->fvals[i1]);
            v2 = _mm_set_ps (p2[1], p2[1], p2[0], p2[0]);
            prod = _mm_mul_ps (v1, v2);

            gf1->fvals[i1]         = pprod[0];
            gf1->fvals[i1 + 1]     = pprod[1];
            gf1->fvals[i1next]     = pprod[2];
            gf1->fvals[i1next + 1] = pprod[3];
          }
        if (N2rem != 0)
          {
            i2    = gf2->ntotal - 1;
            i1    = 2 * idcs[i2];
            p2--;
            
            gf1->fvals[i1]     *= p2[0];
            gf1->fvals[i1 + 1] *= p2[0];
          }
  
        #else  /* !HAVE_SSE */
        size_t i1, i2, i2arr;
        
        for (i2 = 0; i2 < gf2->ntotal; i2++)
          {
            i1 = 2 * idcs[i2];
            
            gf1->fvals[i1]     *= gf2->fvals[i2];
            gf1->fvals[i1 + 1] *= gf2->fvals[i2];
          }
          
        #endif
                
        free (idcs);
      }
      Catch (e)
      {
        EXC_RETHROW_REPRINT (e);
      }
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul (gfunc3 *gf1, gfunc3 const *gf2)
{
  CEXCEPTION_T e = EXC_NONE;
  
  CAPTURE_NULL (gf1);
  CAPTURE_NULL (gf2);
  GFUNC_CHECK_INIT_STATUS (gf1);
  GFUNC_CHECK_INIT_STATUS (gf2);

  Try
  {
    if ((gf1->is_halfcomplex) && (gf2->is_halfcomplex))
      gfunc3_mul_hc (gf1, gf2);
    else if ((!gf1->is_halfcomplex) && (!gf2->is_halfcomplex))
      gfunc3_mul_re (gf1, gf2);
    else if ((gf1->is_halfcomplex) && (!gf2->is_halfcomplex))
      gfunc3_mul_hc_re (gf1, gf2);
    else
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Not implemented for half-complex * real.");
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_vfunc_re (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx;
  float vfval;
  vec3 p;

  vec3_copy (p, gf->xmin);
  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              gf->fvals[idx] *= vfval;
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1] = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_vfunc_hc (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx;
  float re, im;
  cplx vfval;
  vec3 p;

  vec3_copy (p, gf->xmin);
  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx += 2)
            {
              VFUNC_EVAL (vf, vfval, p);
              re = gf->fvals[idx] * vfval[0] - gf->fvals[idx + 1] * vfval[1];
              im = gf->fvals[idx] * vfval[1] + gf->fvals[idx + 1] * vfval[0];
              gf->fvals[idx]     = re;
              gf->fvals[idx + 1] = im;
              
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1] = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_vfunc (gfunc3 *gf, const vfunc *vf)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (vf);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (gf->is_halfcomplex)
    gfunc3_mul_vfunc_hc (gf, vf);
  else
    gfunc3_mul_vfunc_re (gf, vf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div_re (gfunc3 *gf1, gfunc3 const *gf2)
{
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_SSE
      size_t i, N4 = gf2->ntotal / 4, N4rem = gf2->ntotal % 4;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;

      for (i = 0; i < N4; i++, p1++, p2++)
        *p1 = _mm_div_ps (*p1, *p2);
      
      for (i = gf2->ntotal - N4rem; i < gf2->ntotal; i++)
        gf1->fvals[i] /= gf2->fvals[i];

      #else  /* !HAVE_SSE */
      size_t i;
      for (i = 0; i < gf1->ntotal; i++)
        gf1->fvals[i] /= gf2->fvals[i];

      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs;
      CEXCEPTION_T e = EXC_NONE;
      
      Try 
      {
        idcs = gfunc3_subgrid_flatidcs (gf1, gf2);

        #if HAVE_SSE
        size_t i, l = 0, N4 = gf2->ntotal / 4, N4rem = gf2->ntotal % 4;
        __m128 const *p2 = (__m128 const *) gf2->fvals;
        __m128 v1;
        float *pv1 = (float *) &v1;
  
        for (i = 0; i < N4; i++, l += 4, p2++)
          {
            v1 = _mm_set_ps (gf1->fvals[idcs[l + 3]], gf1->fvals[idcs[l + 2]], 
                  gf1->fvals[idcs[l + 1]], gf1->fvals[idcs[l]]);
            v1 = _mm_div_ps (v1, *p2);
  
            gf1->fvals[idcs[l]]     = pv1[0];
            gf1->fvals[idcs[l + 1]] = pv1[1];
            gf1->fvals[idcs[l + 2]] = pv1[2];
            gf1->fvals[idcs[l + 3]] = pv1[3];
          }
        
        for (i = gf2->ntotal - N4rem; i < gf2->ntotal; i++)
          gf1->fvals[idcs[i]] /= gf2->fvals[i];
  
        #else  /* !HAVE_SSE */
        size_t i;
        for (i = 0; i < gf2->ntotal; i++)
          gf1->fvals[idcs[i]] /= gf2->fvals[i];
  
        #endif
        
        free (idcs);
      }
      Catch (e)
      {
        EXC_RETHROW_REPRINT (e);
      }
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div (gfunc3 *gf1, gfunc3 const *gf2)
{
  CEXCEPTION_T e = EXC_NONE;

  CAPTURE_NULL (gf1);
  CAPTURE_NULL (gf2);
  GFUNC_CHECK_INIT_STATUS (gf1);
  GFUNC_CHECK_INIT_STATUS (gf2);
  
  /* TODO: implement mixed div */
  if (gf1->is_halfcomplex != gf2->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Not implemented for mixed function types.");

  Try
  {
    if (gf1->is_halfcomplex)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
    else
      gfunc3_mul_re (gf1, gf2);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div_vfunc_re (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx;
  vec3 p;
  float vfval;

  vec3_copy (p, gf->xmin);
  idx = 0;
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, idx++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              gf->fvals[idx] /= vfval;
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          p[1] += gf->csize[1];
        }
      p[1] = gf->xmin[1];
      p[2] += gf->csize[2];
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div_vfunc (gfunc3 *gf, const vfunc *vf)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (vf);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
  else
    gfunc3_div_vfunc_re (gf, vf);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_scale (gfunc3 *gf, float a)
{
  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  #if HAVE_CBLAS
  cblas_sscal (gf->ntotal, a, gf->fvals, 1);
  
  #elif HAVE_SSE
  size_t i, N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 ma = _mm_set1_ps (a);
  
  for (i = 0; i < N4; i++, p++)
    *p = _mm_mul_ps (*p, ma);

  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    gf->fvals[i] *= a;

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t i;
  for (i = 0; i < gf->ntotal; i++)
    gf->fvals[i] *= a;

  #endif  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_add_constant (gfunc3 *gf, float c)
{
  size_t i;

  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 mc = _mm_set1_ps (c);
  
  for (i = 0; i < N4; i++, p++)
    *p = _mm_add_ps (*p, mc);

  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    gf->fvals[i] += c;

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    gf->fvals[i] += c;

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_translate (gfunc3 *gf, vec3 const s)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (s);
  
  vec3_axpby (1, gf->x0,   1, s);
  vec3_axpby (1, gf->xmin, 1, s);
  vec3_axpby (1, gf->xmax, 1, s);

  if (GFUNC_IS_2D (gf) && (s[2] != 0.0))
    {
      fprintf (stderr, "Warning: z component of shift ignored!\n");

      gf->x0[2] = 0.0;
      gf->xmin[2] = 0.0;
      gf->xmax[2] = 0.0;
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_scale_grid (gfunc3 *gf, float a)
{
  CAPTURE_NULL (gf);

  if (a <= 0)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "a must be positive.");
    
  vec3_scale (gf->x0, a);
  vec3_scale (gf->csize, a);
  vec3_scale (gf->xmin, a);
  vec3_scale (gf->xmax, a);

  if (GFUNC_IS_2D (gf))
    gf->csize[2] = 1.0;

  return;
  
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_dilate (gfunc3 *gf, float a)
{
  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (a <= 0)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "a must be positive.");

  gfunc3_scale_grid (gf, a);

  if (GFUNC_IS_2D (gf))
    gfunc3_scale (gf, 1.0 / a);
  else
    gfunc3_scale (gf, 1.0 / (a * sqrtf (a)));

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set (gfunc3 *gf, idx3 const idx, float const *pval)
{
  size_t fi;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (idx);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  if (!idx3_inside_range (idx, gf->shape))
    EXC_THROW_PRINT (EXC_INDEX);
  
  fi = idx3_flat (idx, gf->shape);

  if (gf->is_halfcomplex)
    {
      fi *= 2;
      gf->fvals[fi]     = pval[0];
      gf->fvals[fi + 1] = pval[1];
    }
  else
    gf->fvals[fi] = *pval;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all_re (gfunc3 *gf, float const *pval)
{
  size_t i;

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 mval = _mm_set1_ps (*pval);
  
  for (i = 0; i < N4; i++, p++)
    *p = mval;

  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    gf->fvals[i] = *pval;

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    gf->fvals[i] = *pval;

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all_hc (gfunc3 *gf, float const *pval)
{
  size_t i, ntotal_flt = 2 * gf->ntotal;

  #if HAVE_SSE
  size_t N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 mval = _mm_set_ps (pval[1], pval[0], pval[1], pval[0]);
  
  for (i = 0; i < N4; i++, p++)
    *p = mval;

  if (N4rem != 0)
    {
      i = ntotal_flt - N4rem;
      gf->fvals[i]     = pval[0];
      gf->fvals[i + 1] = pval[1];
    }
    
  #else  /* !HAVE_SSE */
  for (i = 0; i < ntotal_flt; i += 2)
    {
      gf->fvals[i]     = pval[0];
      gf->fvals[i + 1] = pval[1];
    }

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all (gfunc3 *gf, float const *pval)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (pval);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  if (gf->is_halfcomplex)
    gfunc3_set_all_hc (gf, pval);
  else
    gfunc3_set_all_re (gf, pval);
    
  return;
}
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_make_nonneg (gfunc3 *gf)
{
  size_t i;

  CAPTURE_NULL (gf);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Cannot compare complex values against zero.");

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 m0 = _mm_setzero_ps ();
  
  for (i = 0; i < N4; i++, p++)
    *p = _mm_max_ps (*p, m0);

  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    gf->fvals[i] = fmaxf (gf->fvals[i], 0.0);

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    gf->fvals[i] = fmaxf (gf->fvals[i], 0.0);

  #endif
  return;
}


/*-------------------------------------------------------------------------------------------------*/
// Evaluation
/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_eval (gfunc3 *gf, idx3 const idx)
{
  size_t fi;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (idx);
  GFUNC_CHECK_INIT_STATUS (gf);

  if (!idx3_inside_range (idx, gf->shape))
    EXC_THROW_PRINT (EXC_INDEX);
  
  fi = idx3_flat (idx, gf->shape);

  return gf->fvals[fi];
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_nearest (gfunc3 const *gf, vec3 const pt)
{
  idx3 idx;
  size_t fi;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (pt);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  if (!vec3_between (pt, gf->xmin, gf->xmax))
    return 0.0;

  #if HAVE_SSE
  __m128 const *mpt = (__m128 const *) pt, *mxmin = (__m128 const *) gf->xmin;
  __m128 const *mcs = (__m128 const *) gf->csize;
  __m128 mx = _mm_sub_ps (*mpt, *mxmin);
  mx = _mm_div_ps (mx, *mcs);
  float *idxf = (float *) &mx;
  
  #if HAVE_SSE41
  mx = _mm_round_ps (mx, _MM_FROUND_TO_NEAREST_INT);

  #else
  __m128 m05 = _mm_set1_ps (0.5);  /* use floor (x + 0.5) */
  mx = _mm_add_ps (mx, m05);
  idxf[0] = floorf (idxf[0]);
  idxf[1] = floorf (idxf[1]);
  idxf[2] = floorf (idxf[2]);

  #endif  /* HAVE_SSE41 */
  
  idx[0] = (int) idxf[0];
  idx[1] = (int) idxf[1];
  idx[2] = (int) idxf[2];

  #else  /* !HAVE_SSE */
  int i;
  vec3 x;
  for (i = 0; i < 3; i++)
    {
      x[i] = (pt[i] - gf->xmin[i]) / gf->csize[i];
      idx[i] = (int) (roundf (x[i]));
    }

  #endif  /* HAVE_SSE */

  fi = idx3_flat (idx, gf->shape);

  return gf->fvals[fi];
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_nearest_2d (gfunc3 const *gf, vec3 const pt)
{
  idx3 idx;
  size_t fi;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (pt);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  if ((pt[0] <= gf->xmin[0]) || (pt[0] >= gf->xmax[0]) || 
      (pt[1] <= gf->xmin[1]) || (pt[1] >= gf->xmax[1]))
    return 0.0;

  if (!GFUNC_IS_2D(gf))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Function must be 2-dimensional.");

  #if HAVE_SSE
  __m128 const *mpt = (__m128 const *) pt, *mxmin = (__m128 const *) gf->xmin;
  __m128 const *mcs = (__m128 const *) gf->csize;
  __m128 mx = _mm_sub_ps (*mpt, *mxmin);
  mx = _mm_div_ps (mx, *mcs);
  float *idxf = (float *) &mx;
      
  #if HAVE_SSE41
  mx = _mm_round_ps (mx, _MM_FROUND_TO_NEAREST_INT);

  #else
  __m128 m05 = _mm_set1_ps (0.5);  /* use floor (x + 0.5) */
  mx = _mm_add_ps (mx, m05);
  idxf[0] = floorf (idxf[0]);
  idxf[1] = floorf (idxf[1]);
  idxf[2] = floorf (idxf[2]);

  #endif  /* HAVE_SSE41 */
      
  idx[0] = (int) idxf[0];
  idx[1] = (int) idxf[1];
  idx[2] = (int) idxf[2];

  #else  /* !HAVE_SSE */
  int i;
  vec3 x;
  for (i = 0; i < 3; i++)
    {
      x[i] = (pt[i] - gf->xmin[i]) / gf->csize[i];
      idx[i] = (int) (roundf (x[i]));
    }

  #endif  /*HAVE_SSE */

  fi = idx3_flat (idx, gf->shape);

  return gf->fvals[fi];
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_linear (gfunc3 const *gf, vec3 const pt)
{
  int inc1, inc2;
  idx3 idx;
  size_t fi;
  float fv1, fv2, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (pt);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  if (!vec3_between (pt, gf->xmin, gf->xmax))
    return 0.0;

  /* Compute the cell index and the weights for the function values */

  #if HAVE_SSE
  __m128 const *mpt = (__m128 const *) pt, *mxmin = (__m128 const *) gf->xmin;
  __m128 const *mcs = (__m128 const *) gf->csize;
  __m128 mx = _mm_sub_ps (*mpt, *mxmin);
  mx = _mm_div_ps (mx, *mcs);

  #if HAVE_SSE41
  midxf = _mm_floor_ps (mx);
  idxf = (float *) &midxf;

  #else  
  midxf = mx;
  float *px = (float *) &mx;
  idxf = (float *) &midxf;
  idxf[0] = floorf (px[0]);
  idxf[1] = floorf (px[1]);
  idxf[2] = floorf (px[2]);

  #endif  /* HAVE_SSE41 */

  mwu = _mm_sub_ps (mx, midxf);
  mwl = _mm_set1_ps (1.0f);
  mwl = _mm_sub_ps (mwl, mwu);

  idx[0] = (int) idxf[0];
  idx[1] = (int) idxf[1];
  idx[2] = (int) idxf[2];

  pwl = (float *) &mwl;
  pwu = (float *) &mwu;

  #else  /* !HAVE_SSE */
  int i;
  vec3 x, wl, wu;
  
  for (i = 0; i < 3; i++)
    {
      x[i] = (pt[i] - gf->xmin[i]) / gf->csize[i];
      idx[i] = (int) x[i];
      wu[i] = x[i] - idx[i];
      wl[i] = 1 - wu[i];
    }
    
  pwl = wl;
  pwu = wu;

  #endif  /* HAVE_SSE */
  
  /* Now compute the interpolated value */

  /* TODO: possible to use SIMD here, too? */
  inc2 = gf->shape[1] * gf->shape[0];
  inc1 = gf->shape[0];

  fi = idx[2] * inc2 + idx[1] * inc1 + idx[0];

  fv1  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 1])        * pwl[1];
  fv1 += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 1]) * pwu[1];
  fv1 *= pwl[2];
  fi += inc2;
  fv2  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 1])        * pwl[1];
  fv2 += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 1]) * pwu[1];
  fv2 *= pwu[2];

  return fv1 + fv2;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_linear_2d (gfunc3 const *gf, vec3 const pt)
{
  int inc1;
  idx3 idx;
  size_t fi;
  float fv, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;


  CAPTURE_NULL (gf);
  CAPTURE_NULL (pt);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  if (!GFUNC_IS_2D(gf))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Function must be 2-dimensional.");
  

  if ((pt[0] <= gf->xmin[0]) || (pt[0] >= gf->xmax[0]) || 
      (pt[1] <= gf->xmin[1]) || (pt[1] >= gf->xmax[1]))
    return 0.0;

  /* Compute the cell index and the weights for the function values */

  #if HAVE_SSE
  __m128 const *mpt = (__m128 const *) pt, *mxmin = (__m128 const *) gf->xmin;
  __m128 const *mcs = (__m128 const *) gf->csize;
  __m128 mx = _mm_sub_ps (*mpt, *mxmin);
  mx = _mm_div_ps (mx, *mcs);

  #if HAVE_SSE41
  midxf = _mm_floor_ps (mx);
  idxf = (float *) &midxf;

  #else  
  midxf = mx;
  float *px = (float *) &mx;
  idxf = (float *) &midxf;
  idxf[0] = floorf (px[0]);
  idxf[1] = floorf (px[1]);
  idxf[2] = floorf (px[2]);

  #endif  /* HAVE_SSE41 */

  mwu = _mm_sub_ps (mx, midxf);
  mwl = _mm_set1_ps (1.0f);
  mwl = _mm_sub_ps (mwl, mwu);

  idx[0] = (int) idxf[0];
  idx[1] = (int) idxf[1];
  idx[2] = (int) idxf[2];

  pwl = (float *) &mwl;
  pwu = (float *) &mwu;

  #else  /* !HAVE_SSE */
  int i;
  vec3 x, wl, wu;

  for (i = 0; i < 2; i++)
    {
      x[i] = (pt[i] - gf->xmin[i]) / gf->csize[i];
      idx[i] = (int) x[i];
      wu[i] = x[i] - idx[i];
      wl[i] = 1 - wu[i];
    }
    
  pwl = wl;
  pwu = wu;

  #endif  /* HAVE_SSE */

  /* Now compute the interpolated value */

  inc1 = gf->shape[0];
  fi = idx[1] * inc1 + idx[0];
  
  if (fi + inc1 + 1 >= gf->ntotal)
    {
      vec3_print (pt);
      printf ("idx: %d, %d\n", idx[0], idx[1]);
      printf ("xmax-reldiff: %e %e\n", (gf->xmax[0] - pt[0]) / gf->csize[0], 
      (gf->xmax[1] - pt[1]) / gf->csize[1]);
    }
  
  fv  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 1])        * pwl[1];
  fv += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 1]) * pwu[1];
  
  return fv;
}

/*-------------------------------------------------------------------------------------------------
 * Domain change
 *-------------------------------------------------------------------------------------------------*/

size_t *
gfunc3_subgrid_flatidcs (gfunc3 const *gf, gfunc3 const *gf_sub)
{
  CEXCEPTION_T e = EXC_NONE;
  int i, iz, iy, ix;
  size_t idx, z_offset, y_offset, *active_idcs;
  idx3 off, sfac;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (gf_sub);

  if (!gfunc3_grid_is_subgrid (gf, gf_sub))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_SUBGRID, "2nd argument grid not contained in 1st argument grid.");

  /* Calculation of offset and cell size factors of the subgrid */
  for (i = 0; i < 3; i++)
    {
      off[i]  = (int) (roundf ((gf_sub->xmin[i] - gf->xmin[i]) / gf->csize[i]));
      sfac[i] = (int) (roundf (gf_sub->csize[i] / gf->csize[i]));
    }

  Try
  {
    active_idcs = (size_t *) ali16_malloc (gf_sub->ntotal * sizeof (size_t));
  
    z_offset = off[2] * gf->shape[1] * gf->shape[0];
    for (iz = 0, idx = 0; iz < gf_sub->shape[2]; iz++)
      {
        y_offset = off[1] * gf->shape[0];
        for (iy = 0; iy < gf_sub->shape[1]; iy++)
          {
            for (ix = 0; ix < gf_sub->shape[0]; ix++)
              active_idcs[idx++] = z_offset + y_offset + off[0] + ix * sfac[0];
  
            y_offset += sfac[1] * gf->shape[0];
          }
        z_offset += sfac[2] * gf->shape[1] * gf->shape[0];
      }

    return active_idcs;
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }

  return NULL;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_zeropad (gfunc3 *gf, idx3 const padding)
{
  CEXCEPTION_T e = EXC_NONE;
  int i;
  size_t idx, ntotal_old, *idcs;

  float *fvals_old, nul[2] = {0.0f, 0.0f};
  gfunc3 *gf_tmp;
  
  CAPTURE_NULL (gf);
  CAPTURE_NULL (padding);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  /* gf_tmp is used only to hold the old grid and as input for the subgrid_flatidcs function */
  Try
  {
    gf_tmp = new_gfunc3 ();
    idx3_copy (gf_tmp->shape, gf->shape);
    vec3_copy (gf_tmp->csize, gf->csize);
    vec3_copy (gf_tmp->x0, gf->x0);
    gfunc3_compute_xmin_xmax (gf_tmp);
    gf_tmp->ntotal = gf->ntotal;

    for (i = 0; i < 3; i++)
      gf->shape[i] += 2 * padding[i];

    ntotal_old = gf->ntotal;
    gf->ntotal = idx3_product (gf->shape);
    gfunc3_compute_xmin_xmax (gf);

    idcs = gfunc3_subgrid_flatidcs (gf, gf_tmp);
    
    gfunc3_free (&gf_tmp);

    fvals_old = gf->fvals;
    
    gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));
    gfunc3_set_all (gf, nul);
    for (idx = 0; idx < ntotal_old; idx++)
      gf->fvals[idcs[idx]] = fvals_old[idx];
    
    free (idcs);
    free (fvals_old);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_unpad (gfunc3 *gf, idx3 const padding)
{
  CEXCEPTION_T e = EXC_NONE;
  int i;
  size_t idx, *idcs;

  idx3 ptmp;
  float *fvals_old;
  gfunc3 *gf_tmp;
  
  CAPTURE_NULL (gf);
  CAPTURE_NULL (padding);
  GFUNC_CHECK_INIT_STATUS (gf);

  /* TODO: implement half-complex version */
  if (gf->is_halfcomplex)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");

  /* Check if it is possible to remove 2*padding values from the function */
  idx3_copy (ptmp, padding);
  idx3_scale (ptmp, 2);
  if (!idx3_lt (ptmp, gf->shape))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "padding too large.");

  Try
  {
    /* gf_tmp is used only to hold the old grid and as input for the subgrid_flatidcs function */
    gf_tmp = new_gfunc3 ();
    idx3_copy (gf_tmp->shape, gf->shape);
    vec3_copy (gf_tmp->csize, gf->csize);
    vec3_copy (gf_tmp->x0, gf->x0);
    gfunc3_compute_xmin_xmax (gf_tmp);
    gf_tmp->ntotal = gf->ntotal;

    for (i = 0; i < 3; i++)
      gf->shape[i] -= 2 * padding[i];

    gf->ntotal = idx3_product (gf->shape);
    gfunc3_compute_xmin_xmax (gf);

    idcs = gfunc3_subgrid_flatidcs (gf_tmp, gf);

    gfunc3_free (&gf_tmp);

    fvals_old = gf->fvals;
    gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));

    for (idx = 0; idx < gf->ntotal; idx++)
      gf->fvals[idx] = fvals_old[idcs[idx]];
    
    free (idcs);
    free (fvals_old);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------
 * Old backup code
 *-------------------------------------------------------------------------------------------------*/

// Store the corner coordinates as rows of the matrix m_out. The corners
// are numbered according to the bit pattern, where "0" stands for min 
// and "1" for max. Example: 110 <-> (max,max,min) 
/*
void grid3_vertices(grid3 g, float m_out[8][3], int coord)
{
  int i, j, bitmask; 
  
  VECTOR3_COPY(m_out[0], g.xmin);
  VECTOR3_COPY(m_out[7], g.xmax);
  
  for (j=1; j<7; j++)
  {
    bitmask = 4;
    for (i=0; i<3; i++)
    {
      if ( j & (bitmask >> i) )
        m_out[j][i] = g.xmax[i];
      else
        m_out[j][i] = g.xmin[i];
    }
  }
  
  if (coord == COORD_SPACE)
  {
    for (i=0; i<8; i++)
      MATRIX33_T_VECTOR3(g.rotation, m_out[i]);
  }
  
  return;
}
*/


/*-------------------------------------------------------------------------------------------------*/


/* 

void
gfunc3_restrict (gfunc3 * gf, grid3 g_sub)
{
  // FIXME: flatidcs function is broken
  size_t i, *flatidcs;
  float *dnew;

  if (grid3_equal (gf->grid, g_sub))
    return;

  dnew = (float *) fftwf_alloc_real (g_sub.ntotal);
  MALLOC_CHECK(dnew, g_sub.ntotal * sizeof (float));

  flatidcs = grid3_subgrid_flatidcs (gf->grid, g_sub);

  for (i = 0; i < g_sub.ntotal; i++)
    dnew[i] = gf->fvals[flatidcs[i]];

  free (flatidcs);

  gf->grid = g_sub;
  free (gf->fvals);
  gf->fvals = dnew;

  return;
}
*/

/*-------------------------------------------------------------------------------------------------*/

/*
void
gfunc3_continue (gfunc3 * gf, grid3 g_sup)
{
  // FIXME: flatidcs function is broken
  size_t i, *flatidcs;
  float *dnew;

  if (grid3_equal (gf->grid, g_sup))
    return;

  dnew = (float *) fftwf_alloc_real (g_sup.ntotal);
  MALLOC_CHECK(dnew, g_sup.ntotal * sizeof (float));

  flatidcs = grid3_subgrid_flatidcs (g_sup, gf->grid);

  for (i = 0; i < g_sup.ntotal; i++)
    dnew[flatidcs[i]] = gf->fvals[i];

  free (flatidcs);

  gf->grid = g_sup;
  free (gf->fvals);
  gf->fvals = dnew;

  return;
}
*/

/*-------------------------------------------------------------------------------------------------*/

/*
void
gfunc3_equate (gfunc3 gf1, gfunc3 gf2)
{
  // FIXME: flatidcs function is broken
  size_t i, *flatidcs;

  if (grid3_equal (gf1.grid, gf2.grid))
    {
      cblas_scopy (gf1.grid.ntotal, gf2.data, 1, gf1.data, 1);
      return;
    }

  flatidcs = grid3_subgrid_flatidcs (gf1.grid, gf2.grid);

  for (i = 0; i < gf2.grid.ntotal; i++)
    gf1.data[flatidcs[i]] = gf2.data[i];

  free (flatidcs);

  return;
}
*/

/*-------------------------------------------------------------------------------------------------*/

/*
void
gfunc3_zeropad (gfunc3 * gf, int const *padding)
{
  // FIXME: flatidcs function is broken
  size_t i, *flatidcs;
  float *dnew;
  grid3 g_padded;

  flatidcs = grid3_zeropad_flatidcs (gf->grid, padding, &g_padded);

  dnew = (float *) fftwf_alloc_real (g_padded.ntotal);
  MALLOC_CHECK(dnew, g_padded.ntotal * sizeof (float));

  for (i = 0; i < gf->ntotal; i++)
    dnew[flatidcs[i]] = gf->fvals[i];

  free (flatidcs);

  gf->grid = g_padded;
  free (gf->fvals);
  gf->fvals = dnew;

  return;
}
*/
/*-------------------------------------------------------------------------------------------------*/

/*
void
gfunc3_unpad (gfunc3 * gf, int const *padding)
{
  // FIXME: flatidcs function is broken
  size_t i, *flatidcs;
  float *dnew;
  grid3 g_unpadded;

  flatidcs = grid3_unpad_flatidcs (gf->grid, padding, &g_unpadded);

  dnew = (float *) fftwf_alloc_real (g_unpadded.ntotal);
  MALLOC_CHECK(dnew, g_unpadded.ntotal * sizeof (float));

  for (i = 0; i < g_unpadded.ntotal; i++)
    dnew[i] = gf->fvals[flatidcs[i]];

  free (flatidcs);

  gf->grid = g_unpadded;
  free (gf->fvals);
  gf->fvals = dnew;

  return;
}
*/

/*-------------------------------------------------------------------------------------------------*/

/*
  // Rotated evaluation -- old macro implementation


  float incx[3] = {1.0, 0.0, 0.0};
  float incy[3] = {0.0, 1.0, 0.0};
  float incz[3] = {0.0, 0.0, 1.0};
  float pbufy[3], pbufz[3];
// 
  // Increments according to cell size with rotation
  VECTOR3_SCALE(incx, gf->csize[0]);
  VECTOR3_SCALE(incy, gf->csize[1]);
  VECTOR3_SCALE(incz, gf->csize[2]);

  MATRIX33_T_VECTOR3(gf->rotation, incx);
  MATRIX33_T_VECTOR3(gf->rotation, incy);
  MATRIX33_T_VECTOR3(gf->rotation, incz);

  
  vec3_copy(p, gf->xmin);
  MATRIX33_T_VECTOR3(gf->rotation, p);
  
  idx = 0;
  for (iz=0; iz<gf->shape[2]; iz++)
  {
    vec3_copy(pbufz, p);                           // See below
    for (iy=0; iy<gf->shape[1]; iy++)
    {
      vec3_copy(pbufy, p);                         // Store current 'y' state
      for (ix=0; ix<gf->shape[0]; ix++, idx++)
      {
        gf->fvals[idx] = VFUNC_EVAL(vf, p);
        VECTOR3_ADD(p, incx);
      }
      vec3_copy(p, pbufy);                         // Get 'y' state from before the loop
      VECTOR3_ADD(p, incy);                           // and update
    }
    vec3_copy(p, pbufz);                           // See above
    VECTOR3_ADD(p, incz);
  }
*/
