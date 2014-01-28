/*
 * fft.h -- fast Fourier transform related routines
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

#ifndef __FFT_H__
#define __FFT_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex.h>
#include <fftw3.h>

#include "gfunc3.h"
#include "vfunc.h"


/* Only for testing purposes */
void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf);

void
gfunc3_grid_bwd_reciprocal (gfunc3 *gf);


/* The forward Fourier transform of GF, stored as half-complex array.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 * - Rethrows
 */
void
fft_forward (gfunc3 *gf);


/* The backward Fourier transform of the half-complex GF.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 * - Rethrows
 */
void
fft_backward (gfunc3 *gf);


/* The 3D forward Fourier transform of a complex grid function GF evaluated at NFREQS arbitrary 
 * spatial frequencies FREQS using the NFFT library. The result of the transform is stored in the 
 * array FT, which is assumed to be allocated with a size of at least NFREQS.
 * CAUTION: before calling this function for the first time, the x and z axes of GF must be 
 * swapped due to NFFT internals!
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 * - EXC_GFDIM
 * - EXC_BADARG
 * - Rethrows
 * 
 */
void
nfft_transform (gfunc3 const *gf, float const *freqs, size_t nfreqs, float complex **pftvals);

#endif /* __FFT_H__ */
