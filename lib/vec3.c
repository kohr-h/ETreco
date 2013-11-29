/*
 * vec3.c -- manipulation of 3d vtors and matrices
 * 
 * Copyright 2013 Holger Kohr <kohr@num.uni-sb.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
#include <limits.h>
#include <gsl/gsl_math.h>

#include "vec3.h"
#include "misc.h"
#include "simd.h"

/*-------------------------------------------------------------------------------------------------*
 * Pretty printing
 *-------------------------------------------------------------------------------------------------*/

void
vec3_print (vec3 const v)
{
  printf("/% 10.2e \\\n"
         "|% 10.2e |\n"
         "\\% 10.2e /\n", 
         v[0], v[1], v[2]);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
idx3_print (idx3 const n)
{
  printf("/ %5d \\\n"
         "| %5d |\n"
         "\\ %5d /\n", 
         n[0], n[1], n[2]);

  return;
}

/*-------------------------------------------------------------------------------------------------*
 * Vector properties
 *-------------------------------------------------------------------------------------------------*/

int
idx3_ispos (idx3 const n)
{
  #if HAVE_SSE2
  __m128i const *p = (__m128i const *) n;
  __m128i nul = _mm_setzero_si128 ();
  __m128i pos = _mm_cmpgt_epi32 (*p, nul);
  int *ppos  = (int *) &pos;
  return ppos[0] && ppos[1] && ppos[2];
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (n[i] <= 0)
        return 0;
    }
  
  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_about_eq (vec3 const v1, vec3 const v2, float rel_tol)
{
  int i;
  
  for (i = 0; i < 3; i++)
    {
      if (gsl_fcmp (v1[i], v2[i], rel_tol))
        return 0;
    }
  
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
idx3_eq (idx3 const n1, idx3 const n2)
{
  #if HAVE_SSE2
  __m128i const *p1 = (__m128i const *) n1, *p2 = (__m128i const *) n2;
  __m128i eq = _mm_cmpeq_epi32 (*p1, *p2);
  int *peq = (int *) &eq;
  return peq[0] && peq[1] && peq[2];

  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (n1[i] != n2[i])
        return 0;
    }
  
  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_ge (vec3 const v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) v1, *p2 =  (__m128 const *) v2;
  __m128 ge = _mm_cmpge_ps (*p1, *p2);
  float *pge = (float *) &ge;
  return pge[0] && pge[1] && pge[2];
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (v1[i] < v2[i])
        return 0;
    }

  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_gt (vec3 const v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) v1, *p2 =  (__m128 const *) v2;
  __m128 gt = _mm_cmpgt_ps (*p1, *p2);
  float *pgt = (float *) &gt;
  return pgt[0] && pgt[1] && pgt[2];

  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (v1[i] <= v2[i])
        return 0;
    }

  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_le (vec3 const v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) v1, *p2 =  (__m128 const *) v2;
  __m128 le = _mm_cmple_ps (*p1, *p2);
  float *ple = (float *) &le;
  return ple[0] && ple[1] && ple[2];
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (v1[i] > v2[i])
        return 0;
    }

  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_lt (vec3 const v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) v1, *p2 =  (__m128 const *) v2;
  __m128 lt = _mm_cmplt_ps (*p1, *p2);
  float *plt = (float *) &lt;
  return plt[0] && plt[1] && plt[2];

  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (v1[i] >= v2[i])
        return 0;
    }

  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
idx3_lt (idx3 const n1, idx3 const n2)
{
  #if HAVE_SSE2
  __m128i const *p1 = (__m128i const *) n1, *p2 = (__m128i const *) n2;
  __m128i lt = _mm_cmplt_epi32 (*p1, *p2);
  int *plt = (int *) &lt;
  return plt[0] && plt[1] && plt[2];
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (n1[i] >= n2[i])
        return 0;
    }
  
  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
idx3_le (idx3 const n1, idx3 const n2)
{
  #if HAVE_SSE2
  __m128i const *p1 = (__m128i const *) n1, *p2 = (__m128i const *) n2;
  __m128i gt = _mm_cmpgt_epi32 (*p1, *p2); 
  int *pgt = (int *) &gt;
  return !(pgt[0] || pgt[1] || pgt[2]);

  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if (n1[i] > n2[i])
        return 0;
    }
  
  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*/

int
vec3_between_ints (vec3 const v, idx3 const lb, idx3 const ub)
{
  #if HAVE_SSE2
  __m128 const *p1 = (__m128 const *) v;
  __m128i const *p2 = (__m128i const *) lb, *p3 = (__m128i const *) ub;
  __m128 lbf = _mm_cvtepi32_ps (*p2);
  __m128 gtlb = _mm_cmpgt_ps (*p1, lbf);
  float *pgtlb = (float *) &gtlb;
  if (pgtlb[0] && pgtlb[1] && pgtlb[2])
    return 1;
  
  __m128 ubf = _mm_cvtepi32_ps (*p3);
  __m128 ltub = _mm_cmplt_ps (*p1, ubf);
  float *pltub = (float *) &ltub;
  if (pltub[0] && pltub[1] && pltub[2])
    return 1;

  #else
  int i;
  for (i = 0; i < 3; i++)
    {
      if ( (v[i] <= lb[i]) || (v[i] >= ub[i]) )
        return 0;
    }
  
  #endif
  return 1;
}

/*-------------------------------------------------------------------------------------------------*
 * Scalar Operations
 *-------------------------------------------------------------------------------------------------*/

void
vec3_set_all (vec3 v, float x)
{
  #if HAVE_SSE
  __m128 *p = (__m128 *) v;
  *p = _mm_set1_ps (x);

  #else
  int i;
  for (i = 0; i < 3; i++)
    v[i] = x;
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
idx3_set_all (idx3 n, int a)
{
  #if HAVE_SSE
  __m128i *p = (__m128i *) n;
  *p = _mm_set1_epi32 (a);

  #else
  int i;
  for (i = 0; i < 3; i++)
    n[i] = a;
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_scale (vec3 v, float a)
{
  #if HAVE_SSE
  __m128 *p = (__m128 *) v;
  __m128 ma = _mm_set1_ps (a);
  *p = _mm_mul_ps (*p, ma);
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    v[i] *= a;
    
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
idx3_scale (idx3 n, int a)
{
  #if HAVE_SSE2 && HAVE_SSE41
  __m128i *p = (__m128i *) n;
  __m128i ma = _mm_set1_epi32 (a);
  *p = _mm_mul_epi32 (*p, ma);
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    n[i] *= a;
    
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*
 * Vector Operations
 *-------------------------------------------------------------------------------------------------*/

void
vec3_copy (vec3 dest, vec3 const src)
{
  #if HAVE_SSE
  __m128 *p1 = (__m128 *) dest;
  __m128 const *p2 = (__m128 const *) src;
  *p1 = *p2;

  #else
  int i;
  for (i = 0; i < 3; i++)
    dest[i] = src[i];
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
idx3_copy (idx3 dest, idx3 const src)
{
  #if HAVE_SSE2
  __m128i *p1 = (__m128i *) dest;
  __m128i const *p2 = (__m128i const *) src;
  *p1 = *p2;

  #else
  int i;
  for (i = 0; i < 3; i++)
    dest[i] = src[i];

  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_axpby (float a, vec3 x, float b, vec3 const y)
{
  if ((b == 0.0) || (y == NULL))
  {
    #if HAVE_SSE
    __m128 ma = _mm_set1_ps (a);
    __m128 *p1 = (__m128 *) x;
    *p1 = _mm_mul_ps (*p1, ma);
    
    #else
    int i;
    for (i = 0; i < 3; i++)
      x[i] = a * x[i];
    #endif
  }
  else
    {
      #if HAVE_SSE
      __m128 ma = _mm_set1_ps (a), mb = _mm_set1_ps (b);
      __m128 *p1 = (__m128 *) x;
      __m128 const *p2 = (__m128 const *) y;
      *p1 = _mm_mul_ps (*p1, ma);
      __m128 by = _mm_mul_ps (*p2, mb);
      *p1 = _mm_add_ps (*p1, by);
    
      #else
      int i;
      for (i = 0; i < 3; i++)
        x[i] = a * x[i] + b * y[i];
      
      #endif
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_mul (vec3 v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 *p1 = (__m128 *) v1;
  __m128 const *p2 = (__m128 const *) v2;
  *p1 = _mm_mul_ps (*p1, *p2);
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    v1[i] *= v2[i];
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_mul_int (vec3 v, idx3 const n)
{
  #if HAVE_SSE && HAVE_SSE2
  __m128 *p1 = (__m128 *) v;
  __m128i const *p2 = (__m128i const *) n;
  __m128 nf = _mm_cvtepi32_ps (*p2);
  *p1 = _mm_mul_ps (*p1, nf);
  
  #else
  int i;
  for (i = 0; i < 3; i++)
    v[i] *= n[i];
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_div (vec3 v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 *p1 = (__m128 *) v1;
  __m128 const *p2 = (__m128 const *) v2;
  __m128 v2rcp = _mm_rcp_ps (*p2);
  *p1 = _mm_mul_ps (*p1, v2rcp);

  #else
  int i;
  for (i = 0; i < 3; i++)
    v1[i] /= v2[i];
  
  #endif
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
vec3_div_int (vec3 v, idx3 const n)
{
  #if HAVE_SSE && HAVE_SSE2
  __m128 *p1 = (__m128 *) v;
  __m128i const *p2 = (__m128i const *) n;
  __m128 v2rcp = _mm_cvtepi32_ps (*p2);
  v2rcp = _mm_rcp_ps (v2rcp);
  *p1 = _mm_mul_ps (*p1, v2rcp);

  #else
  int i;
  for (i = 0; i < 3; i++)
    v[i] /= n[i];
  
  #endif
  return;
}


/*-------------------------------------------------------------------------------------------------*
 * Vector contraction
 *-------------------------------------------------------------------------------------------------*/


size_t
idx3_product (idx3 const n)
{
  size_t prod = n[0];
  prod *= n[1] * n[2];
  
  return prod;
}


/*-------------------------------------------------------------------------------------------------*/

float 
vec3_dot (vec3 const v1, vec3 const v2)
{
  #if HAVE_SSE
  __m128 const *p1 = (__m128 const *) v1, *p2 = (__m128 const *) v2;
  __m128 prod = _mm_mul_ps (*p1, *p2);
  float *pprod = (float *) &prod;
  return pprod[0] + pprod[1] + pprod[2];

  #else
  int i;
  float dot = 0.0;
  for (i = 0; i < 3; i++)
    dot += v1[i] * v2[i];
  return dot;
  #endif
}

/*-------------------------------------------------------------------------------------------------*/

float
vec3_product (vec3 const v)
{
  return v[0] * v[1] * v[2];
}

/*-------------------------------------------------------------------------------------------------*/
