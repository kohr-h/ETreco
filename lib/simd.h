/*
 * simd.h -- conditional inclusion of SIMD headers
 * 
 * Copyright 2014 Holger Kohr <kohr@num.uni-sb.de>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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

#ifndef __SIMD_H__
#define __SIMD_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_SSE
#include <xmmintrin.h>
#endif

#if HAVE_SSE2
#include <emmintrin.h>
#endif

#if HAVE_SSE3
#include <pmmintrin.h>
#endif

#if HAVE_SSSE3
#include <pmmintrin.h>
#endif

#if HAVE_SSE41
#include <smmintrin.h>
#endif


#endif  /* __SIMD_H__ */
