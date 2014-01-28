/*
 * et_params.h -- functions to handle ET specific parameters
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

#ifndef __ET_PARAMS_H__
#define __ET_PARAMS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gfunc3.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

typedef struct
{
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
  float mtf_a;
  float mtf_b;
  float mtf_c;
  float mtf_alpha;
  float mtf_beta;
  int mtf_p;
  int mtf_q;
  float acr;

  /* Derived CTF parameters */
  float wave_number;
  float cc1;
  float aper_cutoff;
  
} EtParams;

/*-------------------------------------------------------------------------------------------------*/

EtParams *
new_EtParams (void);

void
EtParams_free (EtParams **pparams);

/*-------------------------------------------------------------------------------------------------*/

void
EtParams_assign_from_file (EtParams *params, const char *fname_et_params);

void
EtParams_print (const EtParams *params);

/*-------------------------------------------------------------------------------------------------*/

#endif // __ET_PARAMS_H__
