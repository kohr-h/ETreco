/*
 * misc.c -- miscellaneous helper functions
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
#include <errno.h>
#include <error.h>
#include <sys/stat.h>

#include "CException.h"
#include "misc.h"
#include "gfunc3.h"
#include "mrc.h"

#define TEMPDIR_STR  "temp/"
#define MAX_COUNT_DIGITS  3

/*-------------------------------------------------------------------------------------------------*/

void *
ali16_malloc (size_t nbytes)
{
  void *alloc_arr;
  
  if (posix_memalign ((void **) &alloc_arr, 16, nbytes))
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_NOMEMORY, "Unable to allocate %lu bytes of memory.", nbytes);
      return NULL;
    }
    
  return alloc_arr;
}

/*-------------------------------------------------------------------------------------------------*/

void
temp_mrc_out (gfunc3 const *gf, char const *mrc_fbasename, int count)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (mrc_fbasename);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  CEXCEPTION_T e = EXC_NONE;
  size_t flen = strlen (mrc_fbasename), templen = strlen (TEMPDIR_STR), extlen = strlen (".mrc");
  char *mrc_fname, *p;


  Try
  {
    /* Try to create the tempdir. If it exists, continue. */
    if (mkdir (TEMPDIR_STR, 0764) && (errno != EEXIST))
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to create directory '%s'.", TEMPDIR_STR);

    mrc_fname = (char *) ali16_malloc (templen + flen + MAX_COUNT_DIGITS + extlen + 1);

    p = mrc_fname;
    sprintf (p, TEMPDIR_STR);
    p += templen;
    sprintf (p, mrc_fbasename);
    p += flen;
    sprintf (p, "%03d", count);
    p += MAX_COUNT_DIGITS;
    sprintf (p, ".mrc");

    printf ("Writing %s\n", mrc_fname);
    
    gfunc3_to_mrc (gf, mrc_fname);
    
    free (mrc_fname);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

#if !GNULIB_DIRNAME
char const *
base_name (char const *filename)
{
  char const *s;
  
  if ((s = strrchr (filename, '/')) == NULL)
    s = filename;
  else
    s++;
    
  return s;
}
#endif

/*-------------------------------------------------------------------------------------------------*/
