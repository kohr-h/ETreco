/*
 * et_params.c -- functions to handle ET specific parameters
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

#include "misc.h"

#include "gfunc3.h"

#include "et_params.h"

// Physical constants and conversion factors
#define EL_REST_ENERGY  510998.928  // [eV]
#define HC              1239.84193  // [eV*nm]

#define ONE_MILLIMETER  1E6   // [nm]
#define ONE_MICROMETER  1E3   // [nm]
#define ONE_KILOVOLT    1E3   // [V]
#define ONE_MILLIRADIAN 1E-3  // [1]

// Amplitude contrast ratio, fixed value seems to be suitable for organic specimens
#define ACR             0.2


/*-------------------------------------------------------------------------------------------------*/

int use_ctf_flag = 0;
int use_mtf_flag = 0;

/*-------------------------------------------------------------------------------------------------*/

EtParams *
new_EtParams (void)
{
  CEXCEPTION_T e = EXC_NONE;
  EtParams *params = NULL;
  
  Try { params = (EtParams *) ali16_malloc (sizeof (EtParams)); }  CATCH_RETURN (e, NULL);
    
  params->acc_voltage             = 0.0;
  params->energy_spread           = 0.0;
  params->magnification           = 0.0;
  params->cs                      = 0.0;
  params->cc                      = 0.0;
  params->aperture                = 0.0;
  params->focal_length            = 0.0;
  params->cond_ap_angle           = 0.0;
  params->defocus_nominal         = 0.0;
  params->mtf_a                   = 0.0;
  params->mtf_b                   = 0.0;
  params->mtf_c                   = 0.0;
  params->mtf_alpha               = 0.0;
  params->mtf_beta                = 0.0;
  params->mtf_p                   = 0;
  params->mtf_q                   = 0;
  params->acr                     = 0.0;
  params->wave_number             = 0.0;
  params->cc1                     = 0.0;
  params->aper_cutoff             = 0.0;

  return params;
}

/*-------------------------------------------------------------------------------------------------*/

void
EtParams_free (EtParams **pparams)
{
  if (pparams == NULL)
    return;
  
  free (*pparams);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
EtParams_assign_from_file (EtParams *params, const char *fname_params)
{
  double dtmp;
  dictionary *dict;

  CAPTURE_NULL_VOID (params);
  CAPTURE_NULL_VOID (fname_params);

  if ((dict = iniparser_load (fname_params)) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read parameters from %s.", fname_params);
      return;
    }


  /* GEOMETRY */

  if ((dtmp = iniparser_getdouble (dict, "optics:magnification", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'magnification' not found in %s.", fname_params);
      return;
    }

  params->magnification = (float) dtmp;
  

  /* CTF */

  /* Ignore CTF if acc_voltage is zero or not present */
  if ((dtmp = iniparser_getdouble (dict, "electronbeam:acc_voltage", 0.0)) == 0.0)
    {
      use_ctf_flag = 0;
      iniparser_freedict (dict);
      return;
    }

  use_ctf_flag = 1;
  params->acc_voltage = (float) dtmp * ONE_KILOVOLT;

  if ((dtmp = iniparser_getdouble (dict, "electronbeam:energy_spread", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'energy_spread' not found in %s.", 
        fname_params);
      return;
    }

  params->energy_spread = (float) dtmp;

  if ((dtmp = iniparser_getdouble (dict, "optics:cs", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cs' not found in %s.", fname_params);
      return;
    }

  params->cs = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:cc", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cc' not found in %s.", fname_params);
      return;
    }

  params->cc = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:aperture", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'aperture' not found in %s.", fname_params);
      return;
    }

  params->aperture = (float) dtmp * ONE_MICROMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:focal_length", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'focal_length' not found in %s.", fname_params);
      return;
    }

  params->focal_length = (float) dtmp * ONE_MILLIMETER;

  if ((dtmp = iniparser_getdouble (dict, "optics:cond_ap_angle", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'cond_ap_angle' not found in %s.", fname_params);
      return;
    }

  params->cond_ap_angle = (float) dtmp * ONE_MILLIRADIAN;

  if ((dtmp = iniparser_getdouble (dict, "optics:defocus_nominal", -1.0)) == -1.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Key 'defocus_nominal' not found in %s.", fname_params);
      return;
    }

  params->defocus_nominal = (float) dtmp * ONE_MICROMETER;

  /* If not found in config file, set to predefined ACR */
  params->acr = (float) iniparser_getdouble (dict, "volume:famp", ACR);

  /* The MTF parameters all have default values */

  use_mtf_flag = 1;
  
  params->mtf_a     = iniparser_getdouble (dict, "detector:mtf_a", 0.0);
  params->mtf_b     = iniparser_getdouble (dict, "detector:mtf_b", 0.0);
  params->mtf_c     = iniparser_getdouble (dict, "detector:mtf_c", 1.0);
  params->mtf_alpha = iniparser_getdouble (dict, "detector:mtf_alpha", 0.0);
  params->mtf_beta  = iniparser_getdouble (dict, "detector:mtf_beta", 0.0);
  params->mtf_p     = iniparser_getint (dict, "detector:mtf_p", 1);  
  params->mtf_q     = iniparser_getint (dict, "detector:mtf_q", 1);  

  /* If a and b are zero, the MTF collapses to a constant. This means 'no MTF'. */
  if ((params->mtf_a == 0.0) && (params->mtf_b == 0.0))
    use_mtf_flag = 0;
   

  iniparser_freedict (dict);


  /* Derived parameters */

  /* Compute relativistic wave number (unit: [1/nm]) */

  /* momentum * c [eV] */
  dtmp = sqrtf (params->acc_voltage * params->acc_voltage + 2 * params->acc_voltage * EL_REST_ENERGY); 
  params->wave_number = (float) 2 * M_PI * dtmp / HC;


  /* Compute constant derived from cc (unit: [nm]) */
  dtmp = 1.0 / (2 * EL_REST_ENERGY); // some factor
  params->cc1 = (float) (1 + 2 * dtmp * params->acc_voltage) 
                / (params->acc_voltage * (1 + dtmp * params->acc_voltage)) * params->cc;

  /* Compute cutoff radius due to aperture */
  params->aper_cutoff = (params->wave_number * params->aperture) / params->focal_length;

  return;
}

/*-------------------------------------------------------------------------------------------------*/


void
EtParams_print (EtParams const *params)
{
  CAPTURE_NULL_VOID (params);

  printf ("\n");
  puts ("ET parameters:");
  puts ("==============\n");

  if (use_ctf_flag)
    {
      puts ("CTF:");
      puts ("----\n");
      
      printf ("acc_voltage     : % 7.2e [V]\n", params->acc_voltage);
      printf ("energy_spread   : % 9.2f [eV]\n", params->energy_spread);
      printf ("cs              : % 7.2e [nm]\n", params->cs);
      printf ("cc              : % 7.2e [nm]\n", params->cc);
      printf ("aperture        : % 9.2f [nm]\n", params->aperture);
      printf ("focal_length    : % 7.2e [nm]\n", params->focal_length);
      printf ("cond_ap_angle   : % 9.5f [rad]\n", params->cond_ap_angle);
      printf ("defocus_nominal : % 9.2f [nm]\n", params->defocus_nominal);
      printf ("mtf_a           : % 9.2f\n", params->mtf_a);
      printf ("mtf_b           : % 9.2f\n", params->mtf_b);
      printf ("mtf_c           : % 9.2f\n", params->mtf_c);
      printf ("mtf_alpha       : % 9.2f\n", params->mtf_alpha);
      printf ("mtf_beta        : % 9.2f\n", params->mtf_beta);
      printf ("mtf_p           : % 9d\n", params->mtf_p);
      printf ("mtf_q           : % 9d\n", params->mtf_q);
      printf ("\n");
      printf ("acr             : % 9.2f\n", params->acr);
      printf ("wave_number     : % 9.2f [1/nm]\n", params->wave_number);
      printf ("cc1             : % 9.2f [nm]\n", params->cc1);
      printf ("aperture cutoff : % 9.2f [1/nm]\n", params->aper_cutoff);
      printf ("\n");
    }
    
  puts ("Geometry:");
  puts ("---------\n");

  printf ("magnification   : % 9.2f\n", params->magnification);
  puts ("\n");
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
