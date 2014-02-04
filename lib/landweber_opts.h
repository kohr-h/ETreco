/*
 * landweber_opts.h -- dispatch options for fwd_op* programs via getopt
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

#ifndef __LANDWEBER_OPTS_H__
#define __LANDWEBER_OPTS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec3.h"
#include "misc.h"

#include "et_operators.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

extern int verbosity_level;
extern int invert_contrast_flag;
extern int normalize_flag;


typedef struct {
  
  /* File names */
  char *fname_in;
  char *fname_out;
  char *fname_params;
  char *fname_tiltangles;
  
  /* Model parameter */
  scattering_model model;

  /* Iteration parameters */
  float relax_param;
  int max_iter;

  /* Normalization parameters */
  idx3 bg_patch_ix0;
  idx3 bg_patch_shape;
  
} LandwOpts;

/*-------------------------------------------------------------------------------------------------*/

/* Create a new LandwOpts structure and return a pointer to it.
 * 
 * Thrown exceptions:  
 * - Rethrows
 */
LandwOpts *
new_LandwOpts (void);

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_print (LandwOpts *opts);


void
LandwOpts_free (LandwOpts **popts);


void
LandwOpts_assign_from_args (LandwOpts *opts, int argc, char **argv);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __LANDWEBER_OPTS_H__ */
