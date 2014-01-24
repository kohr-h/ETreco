/*
 * et_operators.h -- operators specific for ET
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

typedef enum { PROJ_ASSUMPTION, BORN_APPROX } scattering_model;

/*-------------------------------------------------------------------------------------------------
 * ET Forward operator 
 *-------------------------------------------------------------------------------------------------*/

void
et_scattering_projection (gfunc3 const *scatterer, vec3 const angles_deg, RecParams const *rec_p, 
                          gfunc3 *proj_img, scattering_model sct_model);


/*-------------------------------------------------------------------------------------------------*/
