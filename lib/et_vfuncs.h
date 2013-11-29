/*
 * et_vfuncs.h -- vector functions related to ET
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

#ifndef __ET_VFUNCS_H__
#define __ET_VFUNCS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec3.h"

#include "vfunc.h"

#include "et_params.h"

/* TODO: write descriptions, esp. if function is real or complex */

/*-------------------------------------------------------------------------------------------------*/

float
ft_delta (float const *xi, vec3 const gamma);

float
ft_gaussian (float const *xi, vec3 const gamma);

/*-------------------------------------------------------------------------------------------------*/

// without preceding 1/2pi
/* may throw EXC_NULL */
void
vfunc_init_ctf (vfunc *vf, RecParams const *ctf_p);

/* may throw EXC_NULL */
void
vfunc_init_detector_mtf (vfunc *vf, RecParams const *rec_p);

/* may throw EXC_NULL */
void
vfunc_init_detector_recip_mtf (vfunc *vf, RecParams const *rec_p);

/*-------------------------------------------------------------------------------------------------*/


// without preceding 1/2pi
/* may throw EXC_NULL */
void
vfunc_init_ft_rk_single_axis (vfunc *vf, RecParams const *rec_p);

/*-------------------------------------------------------------------------------------------------*/

/* may throw EXC_NULL */
void
vfunc_init_ft_charfun_ball3 (vfunc *vf, float const *pradius);

/*-------------------------------------------------------------------------------------------------*/

/* may throw EXC_NULL */
void
vfunc_init_ft_charfun_cyl3 (vfunc *vf, float const *length_radius);

/*-------------------------------------------------------------------------------------------------*/

/* may throw EXC_NULL */
void
vfunc_init_ft_laplacian (vfunc *vf);

/*-------------------------------------------------------------------------------------------------*/

/* may throw EXC_NULL */
void
vfunc_init_ft_lambda (vfunc *vf, float const *pa);

/*-------------------------------------------------------------------------------------------------*/

/* may throw EXC_NULL */
void
vfunc_init_ft_lambda_v (vfunc *vf, float const *pa);

/*-------------------------------------------------------------------------------------------------*/


#endif  /* __ET_VFUNCS_H__ */
