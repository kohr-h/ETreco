/*
 * params.h -- functions to handle input parameter structures
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

#ifndef __PARAMS_H__
#define __PARAMS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "matvec3.h"

#include "ai_options.h"

// TODO: write descriptions

#define NUM_CTF_LOBES   8

/*-------------------------------------------------------------------------------------------------*/

typedef float (*mollifier_ft_function) (float const *, float const *);

typedef struct
{
  /* Reconstrtuction grid parameters */
  idx3 vol_shape;
  vec3 vol_csize;
  
  /* Necessary CTF parameters */
  float acc_voltage;
  float energy_spread;
  float magnification;
  float cs;
  float cc;
  float aperture;
  float focal_length;
  float cond_ap_angle;
  float defocus_nominal;
  vec3 detector_px_size;
  float mtf_a;
  float mtf_b;
  float mtf_c;
  float mtf_alpha;
  float mtf_beta;
  int mtf_p;
  int mtf_q;
  float acr;
  float tilt_axis;

  /* Derived CTF parameters */
  float wave_number;
  float cc1;
  float aper_cutoff;
  float xover_pts[NUM_CTF_LOBES][2];
  float xover_cspline_coeff[NUM_CTF_LOBES][4];
  
  /* Regularization parameters */
  vec3 gamma;
  float ctf_trunc;

  /* The (assumed real-valued) FT of a mollifier */
  mollifier_ft_function moll_ft;

} RecParams;

/*-------------------------------------------------------------------------------------------------*/

RecParams *
new_RecParams (void);

void
RecParams_free (RecParams **prec_p);

/*-------------------------------------------------------------------------------------------------*/

void
RecParams_assign_from_OptionData (RecParams *rec_p, const OptionData *od);

void
RecParams_print (const RecParams *rec_p);

/*-------------------------------------------------------------------------------------------------*/

#endif // __PARAMS_H__
