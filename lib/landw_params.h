/*
 * landw_params.h -- functions to handle Landweber input parameters
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

#ifndef __LANDW_PARAMS_H__
#define __LANDW_PARAMS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec3.h"

#include "gfunc3.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

typedef struct
{
  /* Volume grid parameters */
  idx3 vol_shape;
  vec3 vol_csize;
  vec3 vol_shift_px;
  
  /* Geometry parameters */
  int tilt_axis;
  float tilt_axis_rotation;
  float tilt_axis_par_shift_px;

  /* Detector parameters */
  vec3 detector_px_size;

} LandwParams;

/*-------------------------------------------------------------------------------------------------*/

LandwParams *
new_LandwParams (void);

void
LandwParams_free (LandwParams **pparams);

/*-------------------------------------------------------------------------------------------------*/

void
LandwParams_assign_from_file (LandwParams *params, const char *fname_params);

void
LandwParams_print (const LandwParams *params);

void
LandwParams_apply_to_volume (LandwParams const *params, gfunc3 *vol);

void
LandwParams_apply_to_proj_stack (LandwParams const *params, gfunc3 *proj_stack);

/*-------------------------------------------------------------------------------------------------*/

#endif // __LANDW_PARAMS_H__
