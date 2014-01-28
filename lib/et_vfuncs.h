/*
 * et_vfuncs.h -- vector functions related to ET
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

void
vfunc_init_ctf (vfunc *vf, EtParams const *params);


void
vfunc_init_ctf_acr (vfunc *vf, EtParams const *params);


void
vfunc_init_detector_mtf (vfunc *vf, EtParams const *params);


void
vfunc_init_detector_recip_mtf (vfunc *vf, EtParams const *params);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __ET_VFUNCS_H__ */
