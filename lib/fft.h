/*
 * fft.h -- fast Fourier transform related routines
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

#ifndef __FFT_H__
#define __FFT_H__

// TODO: write descriptions

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fftw3.h>

#include "gfunc3.h"
#include "vfunc.h"

/* TODO: update descriptions */

// Replace the grid of GF by its reciprocal (half-complex shape)
/* may throw EXC_NULL */
void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf);

// Replace the reciprocal grid of GF (half-complex shape) by the original one 
/* may throw EXC_NULL */
void
gfunc3_hc_grid_bwd_reciprocal (gfunc3 *gf_hc);

// The forward Fourier transform of GF, stored as half-complex array
/* may throw EXC_NULL */
/* may throw EXC_GFINIT */
/* may throw EXC_GFTYPE */
void
fft_forward (gfunc3 *gf);

// The backward Fourier transform of the half-complex GF
/* may throw EXC_NULL */
/* may throw EXC_GFINIT */
/* may throw EXC_GFTYPE */
void
fft_backward (gfunc3 *gf_hc);

#endif // __FFT_H__
