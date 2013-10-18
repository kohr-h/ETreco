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
  CEXCEPTION_T e = EXC_NONE;
  tiltangles *ta;
  
  Try
  {
    ta = (tiltangles *) ali16_malloc (sizeof (tiltangles));
    
    ta->ntilts     = 0;
    ta->nangles    = 0;
    ta->angles_deg = NULL;
    
    return ta;
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return NULL;
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_init (tiltangles *ta, int ntilts, int nangles)
{
  CEXCEPTION_T e = EXC_NONE;
  int i, j;

  CAPTURE_NULL_VOID (ta);

  Try
  {
    if (ntilts <= 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "ntilts must be positive.");

    if (nangles <= 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "nangles must be positive.");

    if (nangles > 3)
      {
        EXC_THROW_CUSTOMIZED_PRINT (EXC_ID1, "Warning: Only the first 3 angles are evaluated!");
        ta->nangles = 3;
      }
    else
      ta->nangles = nangles;

    ta->ntilts = ntilts;

    ta->angles_deg = (float **) ali16_malloc (ntilts * sizeof (float *));

    for (i = 0; i < ntilts; i++)
      {
        ta->angles_deg[i] =  ali16_malloc (nangles * sizeof (float));
        for (j = 0; j < nangles; j++)
          ta->angles_deg[i][j] = 0.0;
      }
  }
  Catch (e)
  {
    /* EXC_ID1 passes with a warning only */
    if (e != EXC_ID1)
      EXC_RETHROW_REPRINT (e);
  }
  
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
    {
      for (i = 0; i < (*pta)->ntilts; i++)
        free ((*pta)->angles_deg[i]);

      free ((*pta)->angles_deg);
    }
  
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

int
scan_line (FILE *fp, int nangles, float *angles)
{
  int i;
  char fmt[13];

  strcpy (fmt, " %f");

  for (i = 1; i < nangles; i++)
    strcat (fmt, " %f");

  strcat (fmt, " \n");

  return fscanf (fp, fmt, &angles[0], &angles[1], &angles[2]);
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_assign_from_file (tiltangles *ta, char const *ta_fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int i, nangles, ntilts;
  FILE *fp;

  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (ta_fname);

  Try
  {
    if ((fp = fopen (ta_fname, "r")) == NULL)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read-only open file %s", ta_fname);

    if (skip_comments (fp) != 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "File %s contains no data", ta_fname);

    if (fscanf (fp, " %d %d ", &ntilts, &nangles) != 2)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Missing info about number of tilts and angles in %s", 
        ta_fname);

    if (skip_comments (fp) != 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "File %s contains no data", ta_fname);

    tiltangles_init (ta, ntilts, nangles);

    for (i = 0; i < ntilts; i++)
      {
        if (scan_line (fp, nangles, ta->angles_deg[i]) == EOF)
          EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "End of file %s reached before all data was read!", 
            ta_fname);
      }

    fclose (fp);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
tiltangles_to_file (tiltangles const *ta, char const *ta_fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int i;
  FILE *fp;

  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (ta_fname);

  Try
  {
    if ((fp = fopen (ta_fname, "w")) == NULL)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write-only open %s.", ta_fname);

    if (fputs ("# Generated by ETreco (tiltangles module)\n", fp) != 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
      
    if (fprintf (fp, "%d %d\n", ta->ntilts, ta->nangles) == 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);

    if (fputs ("#          psi         theta          phi\n", fp) != 0)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);

    for (i = 0; i < ta->ntilts; i++)
      {
        if (fprintf(fp, "       % .4f      % .4f      % .4f\n", 
              ta->angles_deg[i][0], ta->angles_deg[i][1], ta->angles_deg[i][2])
            == 0)
          EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write to %s.", ta_fname);
      }

    fclose (fp);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*
 * Member access
 *-------------------------------------------------------------------------------------------------*/

void
tiltangles_get_angles (tiltangles *ta, float *angles, int index)
{
  int i;
  
  CAPTURE_NULL_VOID (ta);
  CAPTURE_NULL_VOID (angles);
    
  if ((index < 0) || (index >= ta->ntilts))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "index must be between 0 and ntilts-1 (%d)",
      ta->ntilts - 1);

  for (i = 0; i < ta->nangles; i++)
    angles[i] = ta->angles_deg[index][i];
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/
