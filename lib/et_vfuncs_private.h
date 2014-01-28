/*
 * et_vfuncs_private.h -- additional private vector functions related to ET
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

#ifndef __ET_VFUNCS_PRIVATE_H__
#define __ET_VFUNCS_PRIVATE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "et_params.h"

/*-------------------------------------------------------------------------------------------------*/

float
ctf_acr_unscaled_radial (const float t, const EtParams *params);

float
ctf_acr_scaling_function (float t, EtParams const *params);

/*-------------------------------------------------------------------------------------------------*/


#endif  /* __ET_VFUNCS_PRIVATE_H__ */
