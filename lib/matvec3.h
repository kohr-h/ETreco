/*
 * matvec3.h -- manipulation of 3d vectors and matrices
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

#ifndef __MATVEC3_H__
#define __MATVEC3_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Aligned vector data types. Since SSE uses 128 bit registers, a dummy array element is added. */
typedef float __attribute__ ((__aligned__ (16)))  vec3[4];
typedef int   __attribute__ ((__aligned__ (16)))  idx3[4];


/*-------------------------------------------------------------------------------------------------*/
// Pretty printing
/*-------------------------------------------------------------------------------------------------*/

void
vec3_print (vec3 const v);

void
idx3_print (idx3 const n);

/*-------------------------------------------------------------------------------------------------*/
// Vector properties
/*-------------------------------------------------------------------------------------------------*/

int
idx3_ispos (idx3 const n);

int
vec3_about_eq (vec3 const v1, vec3 const v2, float rel_tol);

int
idx3_eq (idx3 const n1, idx3 const n2);

int
vec3_ge (vec3 const v1, vec3 const v2);

int
vec3_gt (vec3 const v1, vec3 const v2);

int
vec3_le (vec3 const v1, vec3 const v2);

int
vec3_lt (vec3 const v1, vec3 const v2);

int
idx3_lt (idx3 const n1, idx3 const n2);

int
idx3_le (idx3 const n1, idx3 const n2);

inline int
idx3_inside_range (idx3 const n, idx3 const shp);

inline int
vec3_between (vec3 const v, vec3 const lb, vec3 const ub);

inline int
vec3_between_ints (vec3 const v, idx3 const lb, idx3 const ub);

/*-------------------------------------------------------------------------------------------------*
 * Scalar Operations
 *-------------------------------------------------------------------------------------------------*/

void
vec3_set_all (vec3 v, float x);

void
idx3_set_all (idx3 n, int a);

void
vec3_scale (vec3 v, float a);

void
idx3_scale (idx3 n, int a);

/*-------------------------------------------------------------------------------------------------*
 * Vector Operations
 *-------------------------------------------------------------------------------------------------*/

void
vec3_copy (vec3 dest, vec3 const src);

void
idx3_copy (idx3 dest, idx3 const src);

inline void
vec3_axpby (float a, vec3 x, float b, vec3 const y);

inline void
vec3_mul (vec3 v1, vec3 const v2);

inline void
vec3_mul_int (vec3 v, idx3 const n);

inline void
vec3_div (vec3 v1, vec3 const v2);

inline void
vec3_div_int (vec3 v, idx3 const n);

/*-------------------------------------------------------------------------------------------------*
 * Vector contraction
 *-------------------------------------------------------------------------------------------------*/

inline size_t
idx3_flat (idx3 const n, idx3 const shp);

inline size_t
idx3_product (idx3 const n);

inline float 
vec3_dot (vec3 const v1, vec3 const v2);

inline float
vec3_product (vec3 const v);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __MATVEC3_H__ */
