/*
 * operators.h -- operations on grid functions
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

#ifndef __OPERATORS_H__
#define __OPERATORS_H__


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gfunc3.h"

#include "et_params.h"

typedef enum { PROJ_ASSUMPTION, BORN_APPROX } scattering_model;
typedef enum { TRAPEZOIDAL } integration_rule; /* TODO: consider other rules */

/*-------------------------------------------------------------------------------------------------
 * ET Forward operator 
 *-------------------------------------------------------------------------------------------------*/

/* Just for testing */
float *
perp_plane_freqs (gfunc3 const *ft_proj_img_grid, vec3 const normal_angles_deg);

float *
ewald_sphere_freqs (gfunc3 const *ft_proj_img_grid, vec3 const normal_angles_deg, float wave_number);


void
et_scattering_projection (gfunc3 const *scatterer_re, gfunc3 const *scatterer_im, 
                          vec3 const angles_deg, RecParams const *rec_p, 
                          gfunc3 *proj_img_re, gfunc3 *proj_img_im, scattering_model sct_model);


/*-------------------------------------------------------------------------------------------------
 * Back-projection operators
 *-------------------------------------------------------------------------------------------------*/

void
xray_backprojection (gfunc3 const *proj_img, vec3 const angles_deg, gfunc3 *volume);

void
xray_backprojection_single_axis (gfunc3 const *proj_img, float const theta_deg, 
                                 RecParams const *rec_p, gfunc3 *volume);


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
