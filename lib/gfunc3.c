/*
 * gfunc3.c -- 3-dimensional grid functions
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
#include <limits.h>
#include <complex.h>

#if HAVE_CBLAS
#include <gsl/gsl_blas.h>
#endif

#include "CException.h"

#include "vec3.h"
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
  gfunc3 *gf = NULL;
  CEXCEPTION_T _e = EXC_NONE;
  
  Try { gf = (gfunc3 *) ali16_malloc (sizeof (gfunc3)); }  CATCH_RETURN (_e, NULL);
        
  gf->is_initialized = FALSE;
  gf->type           = REAL;
  gf->fvals          = NULL;
  return gf;
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
  
  (*pgf)->is_initialized = FALSE;
  (*pgf)->type           = REAL;

  free (*pgf);
  return;
}

/*-------------------------------------------------------------------------------------------------*
 * Structure initialization
 *-------------------------------------------------------------------------------------------------*/

void
gfunc3_init (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  size_t nfloats = 0;

  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (cs);
  CAPTURE_NULL_VOID (shp);

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

  gf->type = gf_type;
  if ((gf_type == HALFCOMPLEX) || (gf_type == COMPLEX))
    nfloats = 2 * gf->ntotal;
  else if (gf_type == REAL)
    nfloats = gf->ntotal;

  gfunc3_compute_xmin_xmax (gf);

  Try { gf->fvals = (float *) ali16_malloc (nfloats * sizeof (float)); }  CATCH_RETURN_VOID (_e);

  gf->is_initialized = TRUE;
  gfunc3_set_all (gf, 0.0);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init_gridonly (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type)
{
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (cs);
  CAPTURE_NULL_VOID (shp);

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

  gf->type = gf_type;
  gfunc3_compute_xmin_xmax (gf);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init_from_foreign_grid (gfunc3 *gf, gfunc3 const *gf_template)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  size_t nfloats = 0;

  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (gf_template);

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
  gf->type = gf_template->type;
  
  if (GFUNC_IS_REAL (gf))
    nfloats = gf_template->ntotal;
  else if (GFUNC_IS_COMPLEX (gf))
    nfloats = 2 * gf_template->ntotal;

  Try { gf->fvals = (float *) ali16_malloc (nfloats * sizeof (float)); }  CATCH_RETURN_VOID (_e);

  gf->is_initialized = TRUE;
  gfunc3_set_all (gf, 0.0);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_stack_pointer (gfunc3 *stack_pt, gfunc3 *stack, int pos)
{
  idx3 shp;
  
  CAPTURE_NULL_VOID (stack_pt);
  CAPTURE_NULL_VOID (stack);
  GFUNC_CAPTURE_UNINIT_VOID (stack);
  
  if ((pos < 0) || (pos >= stack->shape[2]))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Stack position must lie between 0 and # of "
        "images.");
      return;
    }

  idx3_copy (shp, stack->shape);
  shp[2] = 1;
  gfunc3_init_gridonly (stack_pt, stack->x0, stack->csize, shp, stack->type);
  
  gfunc3_print_grid (stack_pt, "stack pointer grid");
  printf ("\nntotal: %lu\n\n", stack_pt->ntotal);
  
  if (stack->type == REAL)
    stack_pt->fvals = &stack->fvals[pos * stack_pt->ntotal];
  else
    stack_pt->fvals = &stack->fvals[2 * pos * stack_pt->ntotal];
    
  stack_pt->is_initialized = TRUE;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_csize (gfunc3 *gf, vec3 const cs)
{
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (cs);

  vec3_copy (gf->csize, cs);
  gfunc3_compute_xmin_xmax (gf);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_x0 (gfunc3 *gf, vec3 const x0)
{
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (x0);

  vec3_copy (gf->x0, x0);
  gfunc3_compute_xmin_xmax (gf);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_compute_xmin_xmax (gfunc3 *gf)
{
  CAPTURE_NULL_VOID (gf);

  int i;
  
  if (gf->type == HALFCOMPLEX)
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

/*-------------------------------------------------------------------------------------------------*
 * Function value initialization
 *-------------------------------------------------------------------------------------------------*/

void
gfunc3_assign_fvals_from_vfunc (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  size_t idx, incr;
  vec3 p;

  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (vf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  incr = GFUNC_IS_COMPLEX (gf) ? 2 : 1;

  gfunc3_set_all (gf, 0.0);
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

/*-------------------------------------------------------------------------------------------------*
 * Screen output
 *-------------------------------------------------------------------------------------------------*/

void
gfunc3_print_grid (gfunc3 const *gf, char const *intro_text)
{
  CAPTURE_NULL_VOID (gf);

  printf ("\n");
  if (intro_text != NULL)
    printf ("%s", intro_text);
  if (GFUNC_IS_2D (gf))
    printf ("    [ 2D ]");
  else
    printf ("    [ 3D ]");
  
  if (gf->type == REAL)
    puts ("    [real]\n");
  else if (gf->type == HALFCOMPLEX)
    puts ("    [half-complex]\n");
  else if (gf->type == COMPLEX)
    puts ("    [complex]\n");

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

/*-------------------------------------------------------------------------------------------------*
 * Attributes
 *-------------------------------------------------------------------------------------------------*/

float
gfunc3_min (gfunc3 const *gf)
{
  size_t i;
  float min;

  CAPTURE_NULL (gf, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);
  
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Minimum not defined for complex functions.");
      return FLT_MAX;
    }
  
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
    min = fminf (min, gf->fvals[i]);
    
  #endif  /* HAVE_SSE */
  return min;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_max (gfunc3 const *gf)
{
  size_t i;
  float max;
  
  CAPTURE_NULL (gf, -FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, -FLT_MAX);

  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Maximum not defined for complex functions.");
      return -FLT_MAX;
    }
  
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
        max = fmaxf (max, gf->fvals[i]);
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

  CAPTURE_NULL (gf, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

  /* TODO: implement complex version; this changes the return value to 'float complex' */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, 
        "Mean value not implemented for complex functions.");
      return FLT_MAX;
    }
  
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

  CAPTURE_NULL (gf, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Variance not implemented for complex functions.");
      return FLT_MAX;
    }
  
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
  CAPTURE_NULL (gf1, 0);
  CAPTURE_NULL (gf2, 0);
  
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

  CAPTURE_NULL (gf, 0);
  CAPTURE_NULL (gf_sub, 0);
  
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

/*-------------------------------------------------------------------------------------------------*
 * Operations
 *-------------------------------------------------------------------------------------------------*/

void
gfunc3_copy_r (gfunc3 *dest, gfunc3 const *src)
{
  CEXCEPTION_T _e = EXC_NONE;

  /* Copy grid */
  vec3_copy (dest->x0,    src->x0);
  vec3_copy (dest->csize, src->csize);
  idx3_copy (dest->shape, src->shape);

  vec3_copy (dest->xmin, src->xmin);
  vec3_copy (dest->xmax, src->xmax);
  dest->ntotal = src->ntotal;
  
  
  /* Init values array */
  if (dest->is_initialized)
    free (dest->fvals);
  
  Try { 
    dest->fvals = (float *) ali16_malloc (src->ntotal * sizeof (float)); 
  } CATCH_RETURN_VOID (_e);
  
  dest->is_initialized = TRUE;
  dest->type = src->type;


  /* Copy values */
  #if HAVE_CBLAS
  cblas_scopy (src->ntotal, src->fvals, 1, dest->fvals, 1);

  #elif HAVE_SSE
  size_t nfloats = src->ntotal;
  size_t i, N4 = nfloats / 4, N4rem = nfloats % 4;
  __m128 *p1 = (__m128 *) dest->fvals;
  __m128 const *p2 = (__m128 const *) src->fvals;

  for (i = 0; i < N4; i++, p1++, p2++)
    *p1 = *p2;
  
  for (i = nfloats - N4rem; i < nfloats; i++)
    dest->fvals[i] = src->fvals[i];

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t nfloats = src->ntotal;
  size_t i;
  for (i = 0; i < nfloats; i++)
    dest->fvals[i] = src->fvals[i];

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_copy_c (gfunc3 *dest, gfunc3 const *src)
{
  CEXCEPTION_T _e = EXC_NONE;

  /* Copy grid */
  vec3_copy (dest->x0,    src->x0);
  vec3_copy (dest->csize, src->csize);
  idx3_copy (dest->shape, src->shape);

  vec3_copy (dest->xmin, src->xmin);
  vec3_copy (dest->xmax, src->xmax);
  dest->ntotal = src->ntotal;
  
  
  /* Init fvals array */
  if (dest->is_initialized)
    free (dest->fvals);
  
  Try { 
    dest->fvals = (float *) ali16_malloc (src->ntotal * sizeof (float complex)); 
  } CATCH_RETURN_VOID (_e);
  
  dest->is_initialized = TRUE;
  dest->type = src->type;


  /* Copy values */
  #if HAVE_CBLAS
  cblas_ccopy (src->ntotal, src->fvals, 1, dest->fvals, 1);

  #elif HAVE_SSE
  size_t nfloats = 2 * src->ntotal;
  size_t i, N4 = nfloats / 4, N4rem = nfloats % 4;
  __m128 *p1 = (__m128 *) dest->fvals;
  __m128 const *p2 = (__m128 const *) src->fvals;

  for (i = 0; i < N4; i++, p1++, p2++)
    *p1 = *p2;
  
  for (i = nfloats - N4rem; i < nfloats; i++)
    dest->fvals[i] = src->fvals[i];

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t nfloats = src->ntotal;
  size_t i;
  for (i = 0; i < nfloats; i++)
    dest->fvals[i] = src->fvals[i];

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_copy (gfunc3 *dest, gfunc3 const *src)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (dest);
  CAPTURE_NULL_VOID (src);
  GFUNC_CAPTURE_UNINIT_VOID (src);

  if (src->type == REAL)
    {
      Try { gfunc3_copy_r (dest, src); }  CATCH_RETURN_VOID (_e);
    }
  else
    {
      Try { gfunc3_copy_c (dest, src); }  CATCH_RETURN_VOID (_e);
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

gfunc3 *
gfunc3_realpart (gfunc3 const *gf, gfunc3 *re)
{
  CEXCEPTION_T  _e = EXC_NONE;
  size_t j;
  float *pre_val;
  float complex *pfval;
  gfunc3 *re_out = NULL;
  
  CAPTURE_NULL (gf, NULL);
  GFUNC_CAPTURE_UNINIT (gf, NULL);
  
  if (re != NULL)
    {
      GFUNC_CAPTURE_UNINIT (re, NULL);
      if (!gfunc3_grids_are_equal (gf, re))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_SUBGRID, "Grids of both arguments must be equal.");
          return NULL; 
        }
      if (GFUNC_IS_COMPLEX (re))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Real part function must be of REAL type.");
          return NULL; 
        }
      re_out = re;
    }
  else
    {
      Try { re_out = new_gfunc3 (); }  CATCH_RETURN (_e, NULL);
      Try { gfunc3_init (re_out, gf->x0, gf->csize, gf->shape, REAL); }  CATCH_RETURN (_e, NULL);
    }
    
  if (GFUNC_IS_REAL (gf))
    {
      Try { gfunc3_copy (re_out, gf); }  CATCH_RETURN (_e, NULL);
    }
  else
    {
      pfval = (float complex *) gf->fvals;
      pre_val = re_out->fvals;
      for (j = 0; j < gf->ntotal; j++)
        *(pre_val++) = crealf (*(pfval++));
    }
  return re_out;
}

/*-------------------------------------------------------------------------------------------------*/

gfunc3 *
gfunc3_imagpart (gfunc3 const *gf, gfunc3 *im)
{
  CEXCEPTION_T  _e = EXC_NONE;
  size_t j;
  float *pim_val;
  float complex *pfval;
  gfunc3 *im_out = NULL;
  
  CAPTURE_NULL (gf, NULL);
  GFUNC_CAPTURE_UNINIT (gf, NULL);
  
  if (im != NULL)
    {
      GFUNC_CAPTURE_UNINIT (im, NULL);
      if (!gfunc3_grids_are_equal (gf, im))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_SUBGRID, "Grids of both arguments must be equal.");
          return NULL; 
        }
      if (GFUNC_IS_COMPLEX (im))
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Imaginary part function must be of REAL type.");
          return NULL; 
        }
      im_out = im;
    }
  else
    {
      Try { im_out = new_gfunc3 (); }  CATCH_RETURN (_e, NULL);
      Try { gfunc3_init (im_out, gf->x0, gf->csize, gf->shape, REAL); }  CATCH_RETURN (_e, NULL);
    }
    
  if (GFUNC_IS_REAL (gf))
    {
      Try { gfunc3_set_all (im_out, 0.0); }  CATCH_RETURN (_e, NULL);
    }
  else
    {
      pfval = (float complex *) gf->fvals;
      pim_val = im_out->fvals;
      for (j = 0; j < gf->ntotal; j++)
        *(pim_val++) = cimagf (*(pfval++));
    }
  return im_out;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_real2complex (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  size_t j;
  float *pfold_val, *fvals_old;
  float complex *pfval;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (!GFUNC_IS_REAL (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Grid function must be of REAL type.");
      return;
    }
  
  fvals_old = gf->fvals;
  Try { 
    gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float complex)); 
  } CATCH_RETURN_VOID (_e);
  
  pfold_val = fvals_old;
  pfval = (float complex *) gf->fvals;
  for (j = 0; j < gf->ntotal; j++)
    *(pfval++) = *(pfold_val++) + 0.0 * I;
  
  gf->type = COMPLEX;
  free (fvals_old);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_imag2complex (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  size_t j;
  float *pfold_val, *fvals_old;
  float complex *pfval;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (!GFUNC_IS_REAL (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Grid function must be of REAL type.");
      return;
    }
  
  fvals_old = gf->fvals;
  Try { 
    gf->fvals = (float *) ali16_malloc (2 * gf->ntotal * sizeof (float)); 
  } CATCH_RETURN_VOID (_e);
  
  pfold_val = fvals_old;
  pfval = (float complex *) gf->fvals;
  for (j = 0; j < gf->ntotal; j++)
    *(pfval++) = 0.0 + *(pfold_val++) * I;
  
  gf->type = COMPLEX;
  free (fvals_old);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
gfunc3_swapxz_r (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  int ix, iy, iz;
  size_t idx, incx;
  float *pfval_sw, *fvals_swapped = NULL;
  
  Try { 
    fvals_swapped = (float *) ali16_malloc (gf->ntotal * sizeof (float)); 
  } CATCH_RETURN_VOID (_e);
  
  pfval_sw = fvals_swapped;

  ISWAP (gf->shape[0], gf->shape[2]);
  FSWAP (gf->csize[0], gf->csize[2]);
  FSWAP (gf->x0[0], gf->x0[2]);
  FSWAP (gf->xmin[0], gf->xmin[2]);
  FSWAP (gf->xmax[0], gf->xmax[2]);
  
  idx = 0;
  incx = gf->shape[1] * gf->shape[2];
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          idx = iy * gf->shape[2] + iz;
          for (ix = 0; ix < gf->shape[0]; ix++, idx += incx)
            *(pfval_sw++) = gf->fvals[idx];
        }
    }
  
  free (gf->fvals);
  gf->fvals = fvals_swapped;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
gfunc3_swapxz_c (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  int ix, iy, iz;
  size_t idx, incx;
  float *fvals_swapped = NULL;
  float complex *pfval, *pfval_sw; 
  
  Try { 
    fvals_swapped = (float *) ali16_malloc (gf->ntotal * sizeof (float complex)); 
  } CATCH_RETURN_VOID (_e);
  
  pfval    = (float complex *) gf->fvals;
  pfval_sw = (float complex *) fvals_swapped;

  ISWAP (gf->shape[0], gf->shape[2]);
  FSWAP (gf->csize[0], gf->csize[2]);
  FSWAP (gf->x0[0], gf->x0[2]);
  FSWAP (gf->xmin[0], gf->xmin[2]);
  FSWAP (gf->xmax[0], gf->xmax[2]);
  
  idx = 0;
  incx = gf->shape[1] * gf->shape[2];
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          idx = iy * gf->shape[2] + iz;
          for (ix = 0; ix < gf->shape[0]; ix++, idx += incx)
            *(pfval_sw++) = pfval[idx];
        }
    }
  
  free (gf->fvals);
  gf->fvals = fvals_swapped;

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
gfunc3_swapxz (gfunc3 *gf)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (GFUNC_IS_2D (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Grid function must be 3D.");
      return;
    }
  if (gf->type == HALFCOMPLEX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Swapping not supported for HALFCOMPLEX functions.");
      return;
    }
  else if (gf->type == REAL)
    {
      Try { gfunc3_swapxz_r (gf); } CATCH_RETURN_VOID (_e);
    }
  else if (gf->type == COMPLEX)
    {
      Try { gfunc3_swapxz_c (gf); } CATCH_RETURN_VOID (_e);
    }

  return;  
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy (float a, gfunc3 *gf1, gfunc3 const *gf2)
{
  size_t nfloats = 0;

  CAPTURE_NULL_VOID (gf1);
  CAPTURE_NULL_VOID (gf2);
  GFUNC_CAPTURE_UNINIT_VOID (gf1);
  GFUNC_CAPTURE_UNINIT_VOID (gf2);

  /* TODO: implement mixed axpy */
  if ( ((GFUNC_IS_REAL (gf1)) && (GFUNC_IS_COMPLEX (gf2))) ||
       ((GFUNC_IS_REAL (gf2)) && (GFUNC_IS_COMPLEX (gf1))) )
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Not implemented for mixed function types.");
      return;
    }

  if (GFUNC_IS_REAL (gf2))
    nfloats = gf2->ntotal;
  else if (GFUNC_IS_COMPLEX (gf2))
    nfloats = 2 * gf2->ntotal;
  
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_CBLAS
      cblas_sscal (nfloats, a, gf1->fvals, 1);
      cblas_saxpy (nfloats, 1.0, gf2->fvals, 1, gf1->fvals, 1);

      #elif HAVE_SSE
      size_t i, N4 = nfloats / 4, N4rem = nfloats % 4;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;
      __m128 ma = _mm_set1_ps (a);

      for (i = 0; i < N4; i++, p1++, p2++)
        {
          *p1 = _mm_mul_ps (*p1, ma);
          *p1 = _mm_add_ps (*p1, *p2);
        }
      
      for (i = nfloats - N4rem; i < nfloats; i++)
        gf1->fvals[i] += a * gf2->fvals[i];

      #else  /* !(HAVE_CBLAS || HAVE_SSE) */
      size_t i;
      for (i = 0; i < nfloats; i++)
        gf1->fvals[i] = a * gf1->fvals[i] + gf2->fvals[i];
        
      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs = NULL;
      CEXCEPTION_T _e = EXC_NONE;
      
      Try { idcs = gfunc3_subgrid_flatidcs (gf1, gf2); }  CATCH_RETURN_VOID (_e);

      #if HAVE_SSE
      size_t i, l = 0, N4 = nfloats / 4, N4rem = nfloats % 4;
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
      
      for (i = nfloats - N4rem; i < nfloats; i++)
        gf1->fvals[idcs[i]] += a * gf2->fvals[i];

      #else  /* !HAVE_SSE */
      size_t i;
      for (i = 0; i < nfloats; i++)
        gf1->fvals[idcs[i]] = a * gf1->fvals[idcs[i]] + gf2->fvals[i];
                  
      #endif
      
      free (idcs);
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_axpy_vfunc_r (float a, gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  float vfval = 0.0, *pfval = gf->fvals;
  vec3 p;

  vec3_copy (p, gf->xmin);
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, pfval++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              *pfval = a * (*pfval) + vfval;
              
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
gfunc3_axpy_vfunc_c (float a, gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  float complex vfval, *pfval = (float complex *) gf->fvals;
  vec3 p;

  vec3_copy (p, gf->xmin);
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, pfval++)
            {
              VFUNC_EVAL (vf, (float *) &vfval, p);
              *pfval = a * (*pfval) + vfval;
              
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
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (vf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  if (GFUNC_IS_REAL (gf))
    gfunc3_axpy_vfunc_r (a, gf, vf);
  else if (GFUNC_IS_COMPLEX (gf))
    gfunc3_axpy_vfunc_c (a, gf, vf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_r (gfunc3 *gf1, gfunc3 const *gf2)
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
      size_t *idcs = NULL;
      CEXCEPTION_T _e = EXC_NONE;
      
      Try { idcs = gfunc3_subgrid_flatidcs (gf1, gf2); }  CATCH_RETURN_VOID (_e);

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

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_c (gfunc3 *gf1, gfunc3 const *gf2)
{
  if (gfunc3_grids_are_equal (gf1, gf2))
    {
      #if HAVE_SSE
      size_t i;
      size_t N2 = gf2->ntotal / 2, N2rem = gf2->ntotal % 2;
      __m128 *p1 = (__m128 *) gf1->fvals;
      __m128 const *p2 = (__m128 const *) gf2->fvals;

      __m128 prod_eq, prod_across;
      float *p1f, *p2f, *pprod_eq = (float *) &prod_eq, *pprod_across = (float *) &prod_across;
      float complex *pf1val, *pf2val;
      
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
          pf1val = (float complex *) &gf1->fvals[2 * (gf2->ntotal - N2rem)];
          pf2val = (float complex *) &gf2->fvals[2 * (gf2->ntotal - N2rem)];
  
          *pf1val *= *pf2val;
        }
        
      #else  /* !HAVE_SSE */
      size_t i;
      float complex *pf1val = (float complex *) gf1->fvals, *pf2val = (float complex *) gf2->fvals;

      for (i = 0; i < gf2->ntotal; i++)
        *(pf1val++) *= *(pf2val++);

      #endif
    }
  else  /* !gfunc3_grids_are_equal (gf1, gf2) */
    {
      size_t *idcs = NULL;
      CEXCEPTION_T _e = EXC_NONE;
      
      Try { idcs = gfunc3_subgrid_flatidcs (gf1, gf2); }  CATCH_RETURN_VOID (_e);

      #if HAVE_SSE
      size_t i, i2, i2arr, i1, i1next, nfloats = 2 * gf2->ntotal;
      size_t N4 = nfloats / 4, N4rem = nfloats % 4;
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
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_c_r (gfunc3 *gf1, gfunc3 const *gf2)
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
      size_t *idcs = NULL;
      CEXCEPTION_T _e = EXC_NONE;
      
      Try { idcs = gfunc3_subgrid_flatidcs (gf1, gf2); }  CATCH_RETURN_VOID (_e);

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
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul (gfunc3 *gf1, gfunc3 const *gf2)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (gf1);
  CAPTURE_NULL_VOID (gf2);
  GFUNC_CAPTURE_UNINIT_VOID (gf1);
  GFUNC_CAPTURE_UNINIT_VOID (gf2);

  Try {
    if ((GFUNC_IS_REAL (gf1)) && (GFUNC_IS_REAL (gf2)))
      gfunc3_mul_r (gf1, gf2);
    else if ((GFUNC_IS_COMPLEX (gf1)) && (GFUNC_IS_COMPLEX (gf2)))
      gfunc3_mul_c (gf1, gf2);
    else if ((GFUNC_IS_COMPLEX (gf1)) && (GFUNC_IS_REAL (gf2)))
      gfunc3_mul_c_r (gf1, gf2);
    else if ((GFUNC_IS_REAL (gf1)) && (GFUNC_IS_COMPLEX (gf2)))
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Not implemented for real x complex.");
  }  CATCH_RETURN_VOID (_e);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_mul_vfunc_r (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  float vfval = 0.0, *pfval = gf->fvals;
  vec3 p;

  vec3_copy (p, gf->xmin);
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              *(pfval++) *= vfval;
              
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
gfunc3_mul_vfunc_c (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  float complex vfval = 0.0, *pfval = (float complex *) gf->fvals;
  vec3 p;

  vec3_copy (p, gf->xmin);
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++)
            {
              VFUNC_EVAL (vf, (float *) &vfval, p);
              *(pfval++) *= vfval;
              
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
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (vf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  if (GFUNC_IS_REAL (gf))
    gfunc3_mul_vfunc_r (gf, vf);
  else if (GFUNC_IS_COMPLEX (gf))
    gfunc3_mul_vfunc_c (gf, vf);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div_r (gfunc3 *gf1, gfunc3 const *gf2)
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
      size_t *idcs = NULL;
      CEXCEPTION_T _e = EXC_NONE;
      
      Try { idcs = gfunc3_subgrid_flatidcs (gf1, gf2); }  CATCH_RETURN_VOID (_e);

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
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div (gfunc3 *gf1, gfunc3 const *gf2)
{
  CAPTURE_NULL_VOID (gf1);
  CAPTURE_NULL_VOID (gf2);
  GFUNC_CAPTURE_UNINIT_VOID (gf1);
  GFUNC_CAPTURE_UNINIT_VOID (gf2);
  
  /* TODO: implement mixed and complex divisions */
  if ((GFUNC_IS_REAL (gf1)) && (GFUNC_IS_REAL (gf2)))
    gfunc3_div_r (gf1, gf2);
  else if ((GFUNC_IS_COMPLEX (gf1)) || (GFUNC_IS_COMPLEX (gf2)))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Not implemented mixed or complex types.");
      return;
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_div_vfunc_r (gfunc3 *gf, const vfunc *vf)
{
  int ix, iy, iz;
  vec3 p;
  float vfval, *pfval = gf->fvals;

  vec3_copy (p, gf->xmin);
  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++)
            {
              VFUNC_EVAL (vf, &vfval, p);
              *(pfval++) /= vfval;
              
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
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (vf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  /* TODO: implement complex version */
  if (GFUNC_IS_REAL (gf))
   gfunc3_div_vfunc_r (gf, vf);
  else if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not implemented.");
      return;
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_scale_realfac (gfunc3 *gf, float a)
{
  size_t ntotal_flt = 0;
  
  if (GFUNC_IS_REAL (gf))
    ntotal_flt = gf->ntotal;
  else if (GFUNC_IS_COMPLEX (gf))
    ntotal_flt = 2 * gf->ntotal;
  
  #if HAVE_CBLAS
  cblas_sscal (ntotal_flt, a, gf->fvals, 1);
  
  #elif HAVE_SSE
  size_t i, N4 = ntotal_flt / 4, N4rem = ntotal_flt % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 ma = _mm_set1_ps (a);
  
  for (i = 0; i < N4; i++, p++)
    *p = _mm_mul_ps (*p, ma);

  for (i = ntotal_flt - N4rem; i < ntotal_flt; i++)
    gf->fvals[i] *= a;

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t i;
  for (i = 0; i < ntotal_flt; i++)
    gf->fvals[i] *= a;

  #endif  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_scale_complfac (gfunc3 *gf, float complex a)
{
  #if HAVE_CBLAS
  cblas_cscal (gf->ntotal, &a, gf->fvals, 1);
  
  #elif HAVE_SSE
  size_t i;
  size_t N2 = gf->ntotal / 2, N2rem = gf->ntotal % 2;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 ma = _mm_set_ps (cimagf (a), crealf (a), cimagf (a), crealf (a));
  __m128 ma_across = _mm_set_ps (crealf (a), cimagf (a), crealf (a), cimagf (a));
  __m128 prod_eq, prod_across;

  float *pf, *pprod_eq = (float *) &prod_eq, *pprod_across = (float *) &prod_across;
  float complex *pfval;
  
  for (i = 0; i < N2; i++, p++)
    {
      prod_eq = _mm_mul_ps (*p, ma);
      prod_across = _mm_mul_ps (*p, ma_across);

      pf = (float *) p;
      pf[0] = pprod_eq[0] - pprod_eq[1];
      pf[1] = pprod_across[0] + pprod_across[1];
      pf[2] = pprod_eq[2] - pprod_eq[3];
      pf[3] = pprod_across[2] + pprod_across[3];
    }
  if (N2rem != 0)
    {
      pfval = (float complex *) &gf->fvals[2 * (gf->ntotal - N2rem)];
      *pfval *= a;
    }

  #else  /* !(HAVE_CBLAS || HAVE_SSE) */
  size_t i;
  float complex *pfval = (float complex *) gf->fvals;
  
  for (i = 0; i < gf->ntotal; i++)
    *(pfval++) *= a;

  #endif  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_scale (gfunc3 *gf, float complex a)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (cimagf (a) == 0.0)
    gfunc3_scale_realfac (gf, crealf (a));
  else
    {
      if (gf->type == REAL)
        {
          Try { gfunc3_real2complex (gf); } CATCH_RETURN_VOID (_e);
        }
      
      gfunc3_scale_complfac (gf, a);
    }
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_add_constant (gfunc3 *gf, float c)
{
  size_t i;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
      return;
    }
  
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
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (s);
  
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
  CAPTURE_NULL_VOID (gf);

  if (a <= 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "a must be positive.");
      return;
    }
    
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
  CEXCEPTION_T _e = EXC_NONE;
  
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  if (a <= 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "a must be positive.");
      return;
    }

  gfunc3_scale_grid (gf, a);

  Try {
    if (GFUNC_IS_2D (gf))
      gfunc3_scale (gf, 1.0 / a);
    else
      gfunc3_scale (gf, 1.0 / (a * sqrtf (a)));
  }  CATCH_RETURN_VOID (_e);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set (gfunc3 *gf, idx3 const idx, float complex val)
{
  size_t fi;

  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (idx);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (!idx3_inside_range (idx, gf->shape))
    {
      EXC_THROW_PRINT (EXC_INDEX);
      return;
    }
  
  fi = idx3_flat (idx, gf->shape);

  if (GFUNC_IS_REAL (gf))
    gf->fvals[fi] = crealf (val);
  else if (GFUNC_IS_COMPLEX (gf))
    {
      fi *= 2;
      gf->fvals[fi]     = crealf (val);
      gf->fvals[fi + 1] = cimagf (val);
    }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all_r (gfunc3 *gf, float val)
{
  size_t i;

  #if HAVE_SSE
  size_t N4 = gf->ntotal / 4, N4rem = gf->ntotal % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 mval = _mm_set1_ps (val);
  
  for (i = 0; i < N4; i++, p++)
    *p = mval;

  for (i = gf->ntotal - N4rem; i < gf->ntotal; i++)
    gf->fvals[i] = val;

  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    gf->fvals[i] = val;

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all_c (gfunc3 *gf, float complex val)
{
  size_t i, nfloats = 2 * gf->ntotal;
  float complex *pfval = (float complex *) gf->fvals;

  #if HAVE_SSE
  size_t N4 = nfloats / 4, N4rem = nfloats % 4;
  __m128 *p = (__m128 *) gf->fvals;
  __m128 mval = _mm_set_ps (cimagf (val), crealf (val), cimagf (val), crealf (val));
  
  for (i = 0; i < N4; i++, p++)
    *p = mval;

  if (N4rem != 0)
    {
      pfval = (float complex *) &gf->fvals [nfloats - N4rem];
      *pfval = val;
    }
    
  #else  /* !HAVE_SSE */
  for (i = 0; i < gf->ntotal; i++)
    *(pfval++) = val;

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_set_all (gfunc3 *gf, float complex val)
{
  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (GFUNC_IS_REAL (gf))
    {
      if (cimagf (val) != 0.0)
        EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Imaginary part is ignored.");
      gfunc3_set_all_r (gf, crealf (val));
    }
  else if (GFUNC_IS_COMPLEX (gf))
    gfunc3_set_all_c (gf, val);
    
  return;
}
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_make_nonneg (gfunc3 *gf)
{
  size_t i;

  CAPTURE_NULL_VOID (gf);
  GFUNC_CAPTURE_UNINIT_VOID (gf);
  
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Cannot compare complex values against zero.");
      return;
    }

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


/*-------------------------------------------------------------------------------------------------*
 * Evaluation
 *-------------------------------------------------------------------------------------------------*/

float
gfunc3_eval (gfunc3 *gf, idx3 const idx)
{
  size_t fi;

  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (idx, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not yet implemented.");
      return FLT_MAX;
    }

  if (!idx3_inside_range (idx, gf->shape))
    {
      EXC_THROW_PRINT (EXC_INDEX);
      return FLT_MAX;
    }
  
  fi = idx3_flat (idx, gf->shape);

  return gf->fvals[fi];
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_nearest (gfunc3 const *gf, vec3 const pt)
{
  idx3 idx;
  size_t fi;

  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not implemented.");
      return FLT_MAX;
    }

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

  #else  /* !HAVE_SSE41 */
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

  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not implemented.");
      return FLT_MAX;
    }

  if ((pt[0] <= gf->xmin[0]) || (pt[0] >= gf->xmax[0]) || 
      (pt[1] <= gf->xmin[1]) || (pt[1] >= gf->xmax[1]))
    return 0.0;

  if (!GFUNC_IS_2D(gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Function must be 2-dimensional.");
      return FLT_MAX;
    }

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
gfunc3_interp_linear_r (gfunc3 const *gf, vec3 const pt)
{
  int inc1, inc2;
  idx3 idx;
  size_t fi;
  float fv1, fv2, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;

  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

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
gfunc3_interp_linear_c (gfunc3 const *gf, vec3 const pt)
{
  int inc1, inc2;
  idx3 idx;
  size_t fi;
  float fv1, fv2, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;

  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);

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
  inc2 = 2 * gf->shape[1] * gf->shape[0];
  inc1 = 2 * gf->shape[0];

  fi = 2 * (idx[2] * inc2 + idx[1] * inc1 + idx[0]);

  fv1  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 2])        * pwl[1];
  fv1 += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 2]) * pwu[1];
  fv1 *= pwl[2];
  fi += inc2;
  fv2  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 2])        * pwl[1];
  fv2 += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 2]) * pwu[1];
  fv2 *= pwu[2];

  return fv1 + fv2;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_linear_2d_r (gfunc3 const *gf, vec3 const pt)
{
  int inc1;
  idx3 idx;
  size_t fi;
  float fv, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;


  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);


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
    
  fv  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 1])        * pwl[1];
  fv += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 1]) * pwu[1];
  
  return fv;
}

/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_interp_linear_2d_c (gfunc3 const *gf, vec3 const pt)
{
  int inc1;
  idx3 idx;
  size_t fi;
  float fv, *pwl, *pwu;

  float *idxf;
  __m128 midxf, mwl, mwu;


  CAPTURE_NULL (gf, FLT_MAX);
  CAPTURE_NULL (pt, FLT_MAX);
  GFUNC_CAPTURE_UNINIT (gf, FLT_MAX);


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

  inc1 = 2 * gf->shape[0];
  fi = 2 * (idx[1] * inc1 + idx[0]);
    
  fv  = (pwl[0] * gf->fvals[fi]        + pwu[0] * gf->fvals[fi + 2])        * pwl[1];
  fv += (pwl[0] * gf->fvals[fi + inc1] + pwu[0] * gf->fvals[fi + inc1 + 2]) * pwu[1];
  
  return fv;
}

/*-------------------------------------------------------------------------------------------------*/

void
store_grid_points_euclidean (gfunc3 const *gf, float *points)
{
  vec3 p;
  int ix, iy, iz;
  float *cur_pt = points;

  vec3_copy (p, gf->xmin);

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, cur_pt += 3)
            {
              /* Just copy over the point */
              cur_pt[0] = p[0];
              cur_pt[1] = p[1];
              cur_pt[2] = p[2];
              
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

void
store_grid_points_polar (gfunc3 const *gf, float *points)
{
  vec3 p;
  int ix, iy, iz;
  float sin_theta, cos_theta, sin_phi, cos_phi;
  float *cur_pt = points;

  vec3_copy (p, gf->xmin);
  cos_theta = cosf (p[1]);
  sin_theta = sinf (p[1]);
  cos_phi = cosf (p[2]);
  sin_phi = sinf (p[2]);
  

  for (iz = 0; iz < gf->shape[2]; iz++)
    {
      for (iy = 0; iy < gf->shape[1]; iy++)
        {
          for (ix = 0; ix < gf->shape[0]; ix++, cur_pt += 3)
            {
              /* p[0] = radius, p[1] = polar angle, p[2] = azimuth angle 
               * The resulting point corresponds to the rotated version of (0,0,radius) 
               * following the 'x' convention in https://de.wikipedia.org/wiki/Eulersche_Winkel */
              cur_pt[0] =   p[0] * sin_theta * sin_phi;
              cur_pt[1] = - p[0] * sin_theta * cos_phi;
              cur_pt[2] =   p[0] * cos_theta;
              
              p[0] += gf->csize[0];
            }
          p[0] = gf->xmin[0];
          
          p[1] += gf->csize[1];
          cos_theta = cosf (p[1]);
          sin_theta = sinf (p[1]);
        }
      p[1] = gf->xmin[1];
      cos_theta = cosf (p[1]);
      sin_theta = sinf (p[1]);

      p[2] += gf->csize[2];
      cos_phi = cosf (p[2]);
      sin_phi = sinf (p[2]);
    }
  
  return;
}

float *
gfunc3_grid_points (gfunc3 const *gf, grid_type gridtype)
{
  CEXCEPTION_T _e = EXC_NONE;
  
  float *points = NULL;

  CAPTURE_NULL (gf, NULL);
  
  Try { 
    points = (float *) ali16_malloc ((gf->ntotal * 3) * sizeof (float)); 
  } CATCH_RETURN (_e, NULL);
  
  if (gridtype == EUCLIDEAN)
    store_grid_points_euclidean (gf, points);
  else if (gridtype == POLAR)
    store_grid_points_polar (gf, points);
    
  return points;
}

/*-------------------------------------------------------------------------------------------------
 * Domain change
 *-------------------------------------------------------------------------------------------------*/

size_t *
gfunc3_subgrid_flatidcs (gfunc3 const *gf, gfunc3 const *gf_sub)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i, iz, iy, ix;
  size_t idx, z_offset, y_offset;
  size_t *active_idcs = NULL;
  idx3 off, sfac;

  CAPTURE_NULL (gf, NULL);
  CAPTURE_NULL (gf_sub, NULL);

  if (!gfunc3_grid_is_subgrid (gf, gf_sub))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_SUBGRID, 
        "2nd argument grid not contained in 1st argument grid.");
      return NULL;
    }

  /* Calculation of offset and cell size factors of the subgrid */
  for (i = 0; i < 3; i++)
    {
      off[i]  = (int) (roundf ((gf_sub->xmin[i] - gf->xmin[i]) / gf->csize[i]));
      sfac[i] = (int) (roundf (gf_sub->csize[i] / gf->csize[i]));
    }

  Try { active_idcs = (size_t *) ali16_malloc (gf_sub->ntotal * sizeof (size_t)); }
  CATCH_RETURN (_e, NULL);
  
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

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_zeropad (gfunc3 *gf, idx3 const padding)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i;
  size_t idx, ntotal_old, *idcs = NULL;

  float *fvals_old, *pfval_old;
  float complex *pfval_c, *pfval_old_c;
  gfunc3 *gf_tmp;
  
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (padding);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_GFTYPE, "Not applicable to HALFCOMPLEX type.");
      return;
    }

  /* gf_tmp is used only to hold the old grid and as input for the subgrid_flatidcs function */
  Try { gf_tmp = new_gfunc3 (); }  CATCH_RETURN_VOID (_e);
  
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

  Try { idcs = gfunc3_subgrid_flatidcs (gf, gf_tmp); }  CATCH_RETURN_VOID (_e);
    
  gfunc3_free (&gf_tmp);

  fvals_old = gf->fvals;
    
  if (gf->type == REAL)
    {
      Try { gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float)); 
      } CATCH_RETURN_VOID (_e);
  
      gfunc3_set_all (gf, 0.0);
      pfval_old = fvals_old;
      for (idx = 0; idx < ntotal_old; idx++, pfval_old++)
        gf->fvals[idcs[idx]] = *pfval_old;
    }
  else if (gf->type == COMPLEX)
    {
      Try { gf->fvals = (float *) ali16_malloc (2 * gf->ntotal * sizeof (float)); 
      } CATCH_RETURN_VOID (_e);
  
      gfunc3_set_all (gf, 0.0);
      
      pfval_old_c = (float complex *) fvals_old;
      for (idx = 0; idx < ntotal_old; idx++, pfval_old_c++)
        {
          pfval_c = (float complex *) &gf->fvals[2 * idcs[idx]];
          *pfval_c = *pfval_old_c;
        }
    }

  free (idcs);
  free (fvals_old);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_unpad (gfunc3 *gf, idx3 const padding)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i;
  size_t idx, *idcs = NULL;

  idx3 ptmp;
  float *fvals_old;
  gfunc3 *gf_tmp;
  
  CAPTURE_NULL_VOID (gf);
  CAPTURE_NULL_VOID (padding);
  GFUNC_CAPTURE_UNINIT_VOID (gf);

  /* TODO: implement complex version */
  if (GFUNC_IS_COMPLEX (gf))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_UNIMPL, "Complex version not implemented.");
      return;
    }

  /* Check if it is possible to remove 2*padding values from the function */
  idx3_copy (ptmp, padding);
  idx3_scale (ptmp, 2);
  if (!idx3_lt (ptmp, gf->shape))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "padding too large.");
      return;
    }

    /* gf_tmp is used only to hold the old grid and as input for the subgrid_flatidcs function */
  Try { gf_tmp = new_gfunc3 (); } CATCH_RETURN_VOID (_e);
  
  idx3_copy (gf_tmp->shape, gf->shape);
  vec3_copy (gf_tmp->csize, gf->csize);
  vec3_copy (gf_tmp->x0, gf->x0);
  gfunc3_compute_xmin_xmax (gf_tmp);
  gf_tmp->ntotal = gf->ntotal;

  for (i = 0; i < 3; i++)
    gf->shape[i] -= 2 * padding[i];

  gf->ntotal = idx3_product (gf->shape);
  gfunc3_compute_xmin_xmax (gf);

  Try { idcs = gfunc3_subgrid_flatidcs (gf_tmp, gf); }  CATCH_RETURN_VOID (_e);

  gfunc3_free (&gf_tmp);

  fvals_old = gf->fvals;
  Try { gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float)); }  CATCH_RETURN_VOID (_e);

  for (idx = 0; idx < gf->ntotal; idx++)
    gf->fvals[idx] = fvals_old[idcs[idx]];
  
  free (idcs);
  free (fvals_old);
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
