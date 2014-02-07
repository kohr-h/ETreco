/*
 * landw_params.c -- functions to handle Landweber input parametes
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "CException.h"

#include "dictionary.h"
#include "iniparser.h"

#include "vec3.h"
#include "misc.h"

#include "gfunc3.h"

#include "landw_params.h"

#define ONE_MICROMETER  1E3   // [nm]

/*-------------------------------------------------------------------------------------------------*/

LandwParams *
new_LandwParams (void)
{
  CEXCEPTION_T e = EXC_NONE;
  LandwParams *params = NULL;
  
  Try { params = (LandwParams *) ali16_malloc (sizeof (LandwParams)); }  CATCH_RETURN (e, NULL);
    
  params->tilt_axis               = 0;
  params->tilt_axis_rotation      = 0.0;
  params->tilt_axis_par_shift_px  = 0.0;

  return params;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwParams_free (LandwParams **pparams)
{
  if (pparams == NULL)
    return;
  
  free (*pparams);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwParams_assign_from_file (LandwParams *params, const char *fname_params)
{
  int itmp;
  float ta_sx, ta_sy;
  double dtmp;
  dictionary *dict;

  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (fname_params);

  if ((dict = iniparser_load (fname_params)) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read parameters from %s.", fname_params);
      return;
    }


  /* VOLUME */

  if ((itmp = iniparser_getint (dict, "volume:nx", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nx' not found in %s.", fname_params);
      return;
    }

  params->vol_shape[0] = itmp;
  
  if ((itmp = iniparser_getint (dict, "volume:ny", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'ny' not found in %s.", fname_params);
      return;
    }

  params->vol_shape[1] = itmp;

  if ((itmp = iniparser_getint (dict, "volume:nz", -1)) == -1)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'nz' not found in %s.", fname_params);
      return;
    }

  params->vol_shape[2] = itmp;

  if ((dtmp = iniparser_getdouble (dict, "volume:voxel_size", FLT_MAX)) == FLT_MAX)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'voxel_size' not found in %s.", fname_params);
      return;
    }

  vec3_set_all (params->vol_csize, (float) dtmp);


  /* These settings override MRC header */

  params->vol_shift_px[0] = (float) iniparser_getdouble (dict, "volume:shift_x", FLT_MAX);
  params->vol_shift_px[1] = (float) iniparser_getdouble (dict, "volume:shift_y", FLT_MAX);
  params->vol_shift_px[2] = (float) iniparser_getdouble (dict, "volume:shift_z", FLT_MAX);

  
  /* GEOMETRY */
  
  /* Single axis: parallel tilt axis shift */
  dtmp = iniparser_getdouble (dict, "geometry:tilt_axis", 0.0);
  if (fabsf (dtmp) < 45.0)  /* Use "x" backprojection variant */
    {
      params->tilt_axis = 0;
      params->tilt_axis_rotation = (float) dtmp;
    }
  else
    {
      params->tilt_axis = 1;
      /* What's missing to +- 90 degrees */
      params->tilt_axis_rotation = (dtmp > 0) ? (float)(dtmp - 90.0) : (float)(-dtmp + 90.0);
    }

  ta_sx = (float) iniparser_getdouble (dict, "geometry:axis_shift_x", 0.0);
  ta_sy = (float) iniparser_getdouble (dict, "geometry:axis_shift_y", 0.0);

  params->tilt_axis_par_shift_px = ta_sx * sinf (params->tilt_axis_rotation * ONE_DEGREE) 
    + ta_sy * cosf (params->tilt_axis_rotation * ONE_DEGREE);


  
  /* DETECTOR */

  /* Overrides MRC header */
  dtmp = iniparser_getdouble (dict, "detector:pixel_size", 0.0);
  if (dtmp != 0.0)
    {
      params->detector_px_size[0] = (float) dtmp * ONE_MICROMETER;
      params->detector_px_size[1] = (float) dtmp * ONE_MICROMETER;
      params->detector_px_size[2] = 1.0;
    }

  iniparser_freedict (dict);

  return;
}

/*-------------------------------------------------------------------------------------------------*/


void
LandwParams_print (LandwParams const *params)
{
  int i;
  
  CAPTURE_NULL_VOID (params);

  printf ("\n");
  puts ("Landweber parameters:");
  puts ("=====================\n");

  puts ("Geometry:");
  puts ("---------\n");
  
  printf ("volume voxel size    : ");
  if (params->vol_csize[0] == FLT_MAX)
    printf ("(from data)\n");
  else
    printf ("% 9.2f [nm]\n", params->vol_csize[0]);
    
  printf ("volume shift         : (");
  for (i = 0; i < 3; i++)
    {
      if (params->vol_shift_px[i] == FLT_MAX)  printf ("(from data)");
      else  printf ("%7.2f", params->vol_shift_px[i]);
      if (i != 2)  printf (", ");
    }
  printf (") [pixels]\n\n");

  printf ("detector pixel size  : % 9.2f [nm]\n", params->detector_px_size[0]);
  printf ("tilt axis rotation   : % 9.2f [degrees]\n", params->tilt_axis_rotation);
  printf ("\n\n");

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwParams_apply_to_volume (LandwParams const *params, gfunc3 *vol)
{
  int i;
  
  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (vol);
  
  for (i = 0; i < 3; i++)
    {
      if (params->vol_shift_px[i] != FLT_MAX)
        vol->x0[i] = params->vol_shift_px[i] * vol->csize[i];
    }
  
  gfunc3_compute_xmin_xmax (vol);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwParams_apply_to_proj_stack (LandwParams const *params, gfunc3 *proj_stack)
{
  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (proj_stack);
  
  if (params->detector_px_size[0] != 0.0)  /* Detector pixel size is set in config file */
    gfunc3_set_csize (proj_stack, params->detector_px_size);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
