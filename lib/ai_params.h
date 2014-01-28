/*
 * ai_params.h -- functions to handle AI specific input parameters
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

#ifndef __AI_PARAMS_H__
#define __AI_PARAMS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec3.h"

#include "gfunc3.h"

#include "et_params.h"

#include "ai_opts.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

#define NUM_CTF_LOBES   8

typedef float (*mollifier_ft_function) (float const *, float const *);

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

  /* Regularization related parameters */
  vec3 gamma;
  float ctf_trunc;
  float xover_pts[NUM_CTF_LOBES][2];
  float xover_cspline_coeff[NUM_CTF_LOBES][4];

  /* The (assumed real-valued) FT of a mollifier */
  mollifier_ft_function moll_ft;

} AiParams;

/*-------------------------------------------------------------------------------------------------*/

AiParams *
new_AiParams (void);

void
AiParams_free (AiParams **pparams);

/*-------------------------------------------------------------------------------------------------*/

void
AiParams_assign_from_AiOpts (AiParams *params, const AiOpts *opts);

void
AiParams_assign_ctftrunc_from_EtParams (AiParams *ai_params, EtParams const *et_params);

void
AiParams_print (const AiParams *params);

void
AiParams_apply_to_volume (AiParams const *params, gfunc3 *vol);

void
AiParams_apply_to_proj_image (AiParams const *params, gfunc3 *proj_img);

/*-------------------------------------------------------------------------------------------------*/

#endif // __AI_PARAMS_H__
