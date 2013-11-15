/*
 * tiltangles.c -- functions to handle tiltangles.txt files 
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CException.h"

#include "misc.h"

#include "tiltangles.h"

/*-------------------------------------------------------------------------------------------------*
 * Memory management
 *-------------------------------------------------------------------------------------------------*/

tiltangles *
new_tiltangles (void)
{
  CEXCEPTION_T _e = EXC_NONE;
  tiltangles *ta = NULL;
  
  Try { ta = (tiltangles *) ali16_malloc (sizeof (tiltangles));  }  CATCH_RETURN (_e, NULL);
    
  ta->ntilts      = 0;
  ta->angles_deg  = NULL;
  
  ta->angle_means = {0.0, 0.0, 0.0};
  ta->angle_vars  = {0.0, 0.0, 0.0};
  ta->nvarying    = 0;
  
  ta->tiltscheme  = GENERIC;
  return ta;
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_init (tiltangles *ta, int ntilts)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i;

  CAPTURE_NULL_VOID (ta);

  if (ntilts <= 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "ntilts must be positive.");
      return;
    }

  ta->ntilts = ntilts;

  Try { ta->angles_deg = (vec3 *) ali16_malloc (ntilts * sizeof (vec3)); }  CATCH_RETURN_VOID (_e);

  for (i = 0; i < ntilts; i++)
    ta->angles_deg[i] = {0.0, 0.0, 0.0};
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_free (tiltangles **pta)
{
  int i;

  if (pta == NULL)
    return;

  if (*pta == NULL)
    return;
  
  if ((*pta)->angles_deg != NULL)
    free ((*pta)->angles_deg);
  
  free (*pta);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*
 * I/O
 *-------------------------------------------------------------------------------------------------*/

int
skip_comments (FILE *fp)
{
  char c;

  while (1)
    {
      if ((c = fgetc (fp)) != '#')
        {
          fseek (fp, -1, SEEK_CUR);
          break;
        }
      else
        {
          /* Skip comment line */
          for (c = fgetc (fp); c != '\n'; c = fgetc (fp))
            if (feof (fp))
              return 1;
        }
    }

  return feof (fp);
}

/*-------------------------------------------------------------------------------------------------*/

void tiltangles_determine_scheme (tiltangles *ta)
{
  
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_assign_from_file (tiltangles *ta, char const *ta_fname)
{
  CEXCEPTION_T _e = EXC_NONE;
  int i, j, nangles, ntilts;
  FILE *fp;

  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (ta_fname);

  if ((fp = fopen (ta_fname, "r")) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read-only open file %s", ta_fname);
      return;
    }

  if (skip_comments (fp) != 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "File %s contains no data", ta_fname);
      fclose (fp);
      return;
    }

  if (fscanf (fp, " %d %d ", &ntilts, &nangles) != 2)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Missing info about number of tilts and angles in %s", 
        ta_fname);
      fclose (fp);
      return;
    }

  if (skip_comments (fp) != 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "File %s contains no data", ta_fname);
      fclose (fp);
      return;
    }

  Try { tiltangles_init (ta, ntilts); }  
  Catch (_e) { EXC_RETHROW_REPRINT (_e);  fclose (fp);  return; }

  for (i = 0; i < ntilts; i++)
    {
      if (fscanf (fp, " %f %f %f \n", &ta->angles_deg[i][0], &ta->angles_deg[i][1], 
        &ta->angles_deg[i][2]) == EOF)
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "End of file %s reached before all data was read!", 
            ta_fname);
          fclose (fp);
          return;
        }
    }

  /* Mean values */
  for (i = 0; i < ntilts; i++)
    {
      for (j = 0; j < 3; j++)
        ta->angle_means[j] += ta->angles_deg[i][j] / ntilts;
    }
    
  /* Mean squares */
  for (i = 0; i < ntilts; i++)
    {
      for (j = 0; j < 3; j++)
        ta->angle_vars[j] += ta->angles_deg[i][j] * ta->angles_deg[i][j] / ntilts;
    }

  /* var = N/(N-1) * (mean squares - squared mean) */
  for (j = 0; j < 3; j++)
    {
      ta->angle_vars[j] = (ta->angle_vars[j] - ta->angle_means[j] * ta->angle_means[j]) 
        * ntilts / (ntilts - 1);
      if (ta->angle_vars[j] > 1E-2)
        ta->nvarying++;
    }
  
  PRINT_VERBOSE ("Mean angles: (%f, %f, %f)\n", ta->angle_means[0], ta->angle_means[1], 
    ta->angle_means[2]);
  PRINT_VERBOSE ("Angle variances: %d\n", ta->angle_vars[0], ta->angle_vars[1], 
    ta->angle_vars[2]);

  fclose (fp);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_to_file (tiltangles const *ta, char const *ta_fname)
{
  int i;
  FILE *fp;

  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (ta_fname);

  if ((fp = fopen (ta_fname, "w")) == NULL)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write-only open %s.", ta_fname);
      return;
    }

  if (fputs ("# Generated by ETreco (tiltangles module)\n", fp) != 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
      fclose (fp);
      return;
    }
    
  if (fprintf (fp, "%d 3\n", ta->ntilts) == 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
      fclose (fp);
      return;
    }

  if (fputs ("#          psi         theta          phi\n", fp) != 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
      fclose (fp);
      return;
    }

  for (i = 0; i < ta->ntilts; i++)
    {
      if (fprintf(fp, "       % .4f      % .4f      % .4f\n", 
            ta->angles_deg[i][0], ta->angles_deg[i][1], ta->angles_deg[i][2])
            == 0)
        {
          EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
          fclose (fp);
          return;
        }
    }

  fclose (fp);
  return;
}

/*-------------------------------------------------------------------------------------------------*
 * Member access
 *-------------------------------------------------------------------------------------------------*/

void
tiltangles_get_angles (tiltangles *ta, vec3 angles, int index)
{
  int i;
  
  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (angles);
    
  if ((index < 0) || (index >= ta->ntilts))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "index must be between 0 and ntilts-1 (%d)",
        ta->ntilts - 1);
      return;
    }

  vec3_copy (angles, ta->angles_deg[index]);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/
