/*
 * vfunc.h -- abstract vector functions
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

#ifndef __VFUNC_H__
#define __VFUNC_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-------------------------------------------------------------------------------------------------*/

/* Vector function struct; to be called with VFUNC_EVAL (always use pointers) */
typedef struct
{
  void (*f) (float const *x, float *pval, void const *params);
  void const *params;

} vfunc;

#define VFUNC_EVAL(pvf, pval, x) ((pvf)->f((x), (pval), (pvf)->params))


/* Aligned float[2] data type for complex data */
typedef float __attribute__ ((__aligned__ (16)))  cplx[2];


/*-------------------------------------------------------------------------------------------------*/

#endif /* __VFUNC_H__ */
