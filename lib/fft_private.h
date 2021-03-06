/*
 * fft_private.h -- fast Fourier transform related routines - private functions
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

#ifndef __FFT_PRIVATE_H__
#define __FFT_PRIVATE_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gfunc3.h"


void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf);

void
gfunc3_grid_bwd_reciprocal (gfunc3 *gf);



#endif /* __FFT_PRIVATE_H__ */
