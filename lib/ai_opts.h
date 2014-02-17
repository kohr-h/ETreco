/*
 * opts.h -- dispatch options for ai_* programs via getopt
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

#ifndef __AI_OPTS_H__
#define __AI_OPTS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vec3.h"
#include "misc.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

extern int verbosity_level;
extern int use_gamma_flag;
extern int truncate_ctf_flag;
extern int normalize_flag;
extern int invert_contrast_flag;
extern int autocenter_vol_flag;
extern int use_lambda_flag;
extern int fft_padding;


/* START and END are dummy names to mark the beginning and end of the enumerated types */
typedef enum { MO_START, DELTA, GAUSSIAN, MO_END } mollifier_type; 

typedef struct {
  
  /* File names */
  char *fname_in;
  char *fname_out;
  char *fname_params;
  char *fname_tiltangles;
  
  /* Regularization parameters */
  float gamma;
  float ctf_trunc;
  mollifier_type moll_type;
  
  /* Normalization parameters */
  idx3 bg_patch_ix0;
  idx3 bg_patch_shape;
  
  /* Data subset parameters */
  int num_images;
  int start_index;
  
  /* Feature reco parameters */
  float neg_lapl_param;
  float lambda_pow;
  
} AiOpts;

/*-------------------------------------------------------------------------------------------------*/

/* Create a new AiOpts structure and return a pointer to it.
 * 
 * Thrown exceptions:  
 * - Rethrows
 */
AiOpts *
new_AiOpts (void);

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_print (AiOpts *opts);


void
AiOpts_free (AiOpts **popts);


void
AiOpts_assign_from_args (AiOpts *opts, int argc, char **argv);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __AI_OPTS_H__ */
