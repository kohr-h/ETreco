/*
 * tiltangles.h -- functions to handle tiltangles.txt files 
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

#ifndef __TILTANGLES_H__
#define __TILTANGLES_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>


typedef struct
{
  int ntilts;
  int nangles;

  float **angles_deg;

} tiltangles;


/*-------------------------------------------------------------------------------------------------*
 * Memory management
 *-------------------------------------------------------------------------------------------------*/


// Allocate memory for tiltangles struct
tiltangles *
new_tiltangles (void);

// Allocate memory for tiltangles members
void
tiltangles_init (tiltangles *ta, int ntilts, int nangles);

void
tiltangles_free (tiltangles **pta);

/*-------------------------------------------------------------------------------------------------*
 * Member access
 *-------------------------------------------------------------------------------------------------*/

void
tiltangles_get_angles (tiltangles *ta, float *angles, int index);

/*-------------------------------------------------------------------------------------------------*
 * I/O
 *-------------------------------------------------------------------------------------------------*/

void
tiltangles_assign_from_file (tiltangles *ta, char const *ta_fname);

void
tiltangles_to_file (tiltangles const *ta, char const *ta_fname);

/*-------------------------------------------------------------------------------------------------*
 * Other stuff
 *-------------------------------------------------------------------------------------------------*/

void
compute_rotated_basis (vec3 const angles_deg, vec3 om_x, vec3 om_y, vec3 om_z);

/*-------------------------------------------------------------------------------------------------*/

#endif  // __TILTANGLES_H__
