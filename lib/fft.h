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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fftw3.h>

#include "gfunc3.h"
#include "vfunc.h"


/* Standard Fourier Transform functions */


/* Replace the grid of GF by its reciprocal (half-complex shape).
 * 
 * Thrown exceptions:
 * - EXC_NULL
 */
void
gfunc3_grid_fwd_reciprocal (gfunc3 *gf);


/* Replace the reciprocal grid of GF (half-complex shape) by the original one.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 */
void
gfunc3_hc_grid_bwd_reciprocal (gfunc3 *gf_hc);


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
fft_backward (gfunc3 *gf_hc);


/* The 3D forward Fourier transform of a complex function defined by F_RE (real part) and F_IM 
 * (imaginary part) for NFREQS arbitrary spatial frequencies FREQS using the NFFT library. Result
 * of the transform is stored in the arrays FT_RE (real part) and FT_IM (imaginary part), both of 
 * which are assumed to be allocated and of size NFREQS.
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
nfft3_transform (gfunc3 const *f_re, gfunc3 const *f_im, float const *freqs, size_t nfreqs, 
                 float *ft_re, float *ft_im);

#endif /* __FFT_H__ */
