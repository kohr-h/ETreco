/*
 * operators.h -- operations on grid functions - private header
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

#ifndef __OPERATORS_PRIVATE_H__
#define __OPERATORS_PRIVATE_H__


#include "vec3.h"

#include "gfunc3.h"


float *
perp_plane_freqs (gfunc3 const *ft_proj_img_grid, vec3 const normal_angles_deg);


#endif // __OPERATORS_PRIVATE_H__
