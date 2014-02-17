/*
 * operators.h -- operations on grid functions
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

#ifndef __OPERATORS_H__
#define __OPERATORS_H__


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex.h>

#include "gfunc3.h"

typedef enum { TRAPEZOIDAL } integration_rule; /* TODO: consider other rules */

/*-------------------------------------------------------------------------------------------------
 * X-ray transform operators
 *-------------------------------------------------------------------------------------------------*/

void
xray_projection (gfunc3 const *volume, vec3 const angles_deg, gfunc3 *proj_img);

void
xray_backprojection (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume, float weight);

void
xray_backprojection_single_axis (gfunc3 const *proj_img, float const theta_deg, int axis,
                                 float tilt_axis_par_shift_px, gfunc3 *volume, float weight);


/*-------------------------------------------------------------------------------------------------
 * 2D Image manipulation
 *-------------------------------------------------------------------------------------------------*/

void
image_rotation (gfunc3 *proj_img, float const psi_deg);

void
histogram_normalization (gfunc3 *proj_img, idx3 bg_ix0, idx3 const bg_shp);

/*-------------------------------------------------------------------------------------------------
 * Miscallaneous
 *-------------------------------------------------------------------------------------------------*/

float
lp_integral (gfunc3 const *gf, integration_rule rule);

#endif // __OPERATORS_H__
