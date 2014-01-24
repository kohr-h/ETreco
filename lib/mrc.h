/*
 * mrc.h -- I/O support for MRC files
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gfunc3.h"

typedef enum { VOLUME, STACK } mrc_ftype;

// TODO: descriptions

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init_mrc (gfunc3 *gf, char const *mrc_fname, FILE **pfp_in, int *pn_img, mrc_ftype type);

void
gfunc3_read_from_stack (gfunc3 *gf, FILE *fp, int stackpos);

void
gfunc3_to_mrc (gfunc3 const *gf, char const *mrc_fname);

void
temp_mrc_out (gfunc3 const *gf, char const *mrc_fbasename, int count);

/*-------------------------------------------------------------------------------------------------*/
