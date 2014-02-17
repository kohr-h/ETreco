/*
 * gfunc3.h -- 3-dimensional grid functions - private header
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

#ifndef __GFUNC3_PRIVATE_H__
#define __GFUNC3_PRIVATE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex.h>

#include "vec3.h"

#include "vfunc.h"


/* Initialize the grid of GF with the provided parameters X0, CS and SHP as well as the function 
 * type. Do not allocate the data array.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_init_gridonly (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type);


/* Approximate the value of GF at PT by linear interpolation on the grid. This is an optimized 2D 
 * version of the generic function.  This version is for REAL functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_GFDIM
 */
float
gfunc3_interp_linear_2d_r (gfunc3 const *gf, vec3 const pt);


/* Approximate the value of GF at PT by linear interpolation on the grid. This is an optimized 2D 
 * version of the generic function.  This version is for COMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_GFDIM
 */
float
gfunc3_interp_linear_2d_c (gfunc3 const *gf, vec3 const pt);



#endif /* __GFUNC_H__ */
