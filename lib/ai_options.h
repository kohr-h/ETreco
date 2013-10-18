/*
 * ai_options.h -- dispatch options for ai_* programs via getopt
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

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "matvec3.h"
#include "misc.h"

// TODO: write descriptions

/*-------------------------------------------------------------------------------------------------*/

extern int verbosity_level;
extern int use_ctf_flag;
extern int use_mtf_flag;
extern int use_gamma_flag;
extern int truncate_ctf_flag;
extern int normalize_flag;
extern int invert_contrast_flag;
extern int autocenter_vol_flag;
extern int use_lambda_flag;
extern int fft_padding;


/* START and END are dummy names to mark the beginning and end of the enumerated types */
typedef enum { M_START, DELTA, GAUSSIAN, M_END } mollifier_type; 
typedef enum { T_START, SINGLE_AXIS, DOUBLE_AXIS, CONICAL, T_END } tiltscheme;

typedef struct {
  
  /* File names */
  char *fname_in;
  char *fname_out;
  char *fname_reco_params;
  char *fname_tiltangles;
  
  /* Tilt scheme */
  tiltscheme tilting_scheme;
  
  /* Additional files for double axis */
  char *fname_in_axis2;
  char *fname_reco_params_axis2;
  char *fname_tiltangles_axis2;
  
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
  
} OptionData;

/*-------------------------------------------------------------------------------------------------*/

/* Create a new OptionData structure and return a pointer to it.
 * 
 * Thrown exceptions:  
 * - Rethrows
 */
OptionData *
new_OptionData (void);

/*-------------------------------------------------------------------------------------------------*/

void 
print_version_etc (char const *progname);

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_print (OptionData *od);


void
OptionData_free (OptionData **pod);


void
OptionData_assign_from_args (OptionData *od, int argc, char **argv);

/*-------------------------------------------------------------------------------------------------*/

#endif  /* __OPTIONS_H__ */
