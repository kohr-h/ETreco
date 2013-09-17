/*
 * mrc.c -- I/O support for MRC files
 * 
 * Copyright 2013 Holger Kohr <kohr@num.uni-sb.de>
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stddef.h>
#include <byteswap.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>

#include "CException.h"
#include "gfunc3.h"
#include "matvec3.h"
#include "misc.h"
#include "mrc.h"


#define MRC_HEADER_BYTES  1024

static int bswap_flag = 0;

/*-------------------------------------------------------------------------------------------------*/

typedef union {
  int32_t _int;
  float _float;
} d32;

/*-------------------------------------------------------------------------------------------------*/

void
read_uchar (unsigned char *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (pval, sizeof (unsigned char), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_int16 (int16_t *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (pval, sizeof (int16_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    *pval = bswap_16 (*pval);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_int16_arr (int16_t *vals, size_t nmemb, FILE *fp, size_t pos)
{
  size_t i;
  
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (vals, nmemb * sizeof (int16_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    {
      for (i = 0; i < nmemb; i++)
        vals[i] = bswap_16 (vals[i]);
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_int32 (int32_t *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (pval, sizeof (int32_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    *pval = bswap_32 (*pval);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_int32_arr (int32_t *vals, size_t nmemb, FILE *fp, size_t pos)
{
  size_t i;
  
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (vals, nmemb * sizeof (int32_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    {
      for (i = 0; i < nmemb; i++)
        vals[i] = bswap_32 (vals[i]);
    }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_float (float *pval, FILE *fp, size_t pos)
{
  d32 swapval;
  
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (&swapval._float, sizeof (float), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    swapval._int = bswap_32 (swapval._int);
  
  *pval = swapval._float;
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
read_float_arr (float *vals, size_t nmemb, FILE *fp, size_t pos)
{
  size_t i;
  d32 swapval;
  
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fread (vals, nmemb * sizeof (float), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  if (bswap_flag)
    {
      for (i = 0; i < nmemb; i++)
        {
          swapval._float = vals[i];
          swapval._int = bswap_32 (swapval._int);
          vals[i] = swapval._float;
        }
    }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_int16 (const int16_t *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (pval, sizeof (int16_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_int32 (const int32_t *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (pval, sizeof (int32_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_int32_arr (const int32_t *vals, size_t nmemb, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (vals, nmemb * sizeof (int32_t), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_float (float const *pval, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (pval, sizeof (float), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_float_arr (float const *vals, size_t nmemb, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (vals, nmemb * sizeof (float), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
      
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_char_arr (char const *str, size_t nmemb, FILE *fp, size_t pos)
{
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  if (fwrite (str, nmemb * sizeof (char), 1, fp) == 0)
    EXC_THROW_PRINT (EXC_IO);
      
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
write_zero_bytes (size_t num, FILE *fp, size_t pos)
{
  const char zero = '\0';
  size_t i;
  
  if (fseek (fp, pos, SEEK_SET) != 0) 
    EXC_THROW_PRINT (EXC_IO);
  
  for (i = 0; i < num; i++)
  {
    if (fwrite (&zero, 1, 1, fp) == 0)
      EXC_THROW_PRINT (EXC_IO);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_init_mrc (gfunc3 *gf, char const *mrc_fname, FILE **pfp_in, int *pn_img, mrc_ftype type)
{
  CEXCEPTION_T e = EXC_NONE;
  
  /* Parameters defined by the MRC specification */
  unsigned char endian;
  int32_t mode, next;
  float amin, amax, amean;

  /* Auxiliary variables */
  float dmin, dmax, a, b;
  size_t i, dtype_size, total_bytes, filesize;

  /* Variables for temporary storage */
  int16_t *i16_arr;
  int32_t ibuf_arr[3];
  
  FILE *fp;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (mrc_fname);

  /* Try opening the file */
  fp = fopen (mrc_fname, "r");
  if (fp == NULL)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read-only open %s.", mrc_fname);

  /* Test if file contains at least enough Bytes for a standard MRC header */
  if (fseek (fp, MRC_HEADER_BYTES, SEEK_SET) != 0) 
    EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read standard %d Byte MRC header from %s", 
      MRC_HEADER_BYTES, mrc_fname);

  /* Test endianness and set bswap_flag accordingly. 
   * FIXME: condition for big endianness is probably wrong. Test this with true data!
   */
  
  read_uchar (&endian, fp, 212);
  if (endian == 68)       /* Little endian */
    bswap_flag = 0;
  else if (endian == 65)  /* Big endian (???) */
    bswap_flag = 1;
  else
    fprintf (stderr, "Warning: Unable to determine endianness. Assuming little-endian.\n");

  /* Read (nx,ny,nz), mode and next (length of extended header).
   * Try if file has exactly 1024 + next + nx*ny*nz * sizeof(<data type>) Bytes, where 
   * <data type> is int16 or float16
   */
   
  read_int32_arr (ibuf_arr, 3, fp,  0);    // Bytes    0 -- 12: shape
  read_int32     (&mode,       fp, 12);    // Bytes   12 -- 16: mode
  read_int32     (&next,       fp, 92);    // Bytes   92 -- 96: next

  switch (mode)
    {
      case 1: dtype_size = 2; break;  /* 16 bit integers (short) */
      case 2: dtype_size = 4; break;  /* 32 bit floating point numbers (float) */
      case 4: dtype_size = 8; break;  /* 2 * 32 bit floats (complex) */
      default: 
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Data mode %d found in %s not supported.", mode, mrc_fname);
      return;
    }

  total_bytes = MRC_HEADER_BYTES + next + idx3_product (ibuf_arr) * dtype_size;

  if (fseek (fp, total_bytes, SEEK_SET) != 0) 
    EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to read data from %s", mrc_fname);
  
  fseek (fp, 0, SEEK_END);
  filesize = ftell (fp);  
  
  if (total_bytes != filesize)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Size of file %s (%lu Byte) does not match header + "
      "data size according to MRC header (%lu Byte).", mrc_fname, filesize, total_bytes);


  /* Now the real work can start. */
  
  if (gf->is_initialized)
    free (gf->fvals);

  Try
  {
    read_int32_arr (gf->shape, 3, fp,  0);  // Bytes    0 -- 12: shape
  
    if (type == STACK)
      {
        if (pn_img != NULL)
          *pn_img = gf->shape[2];
          
        gf->shape[2] = 1;
      }
    
    read_float_arr (gf->csize, 3, fp, 40);  // Bytes   40 -- 52: total grid size
    vec3_div_int (gf->csize, gf->shape);
    
    read_float     (&amin,        fp, 76);  // Bytes   76 -- 80: amin
    read_float     (&amax,        fp, 80);  // Bytes   80 -- 84: amax
    read_float     (&amean,       fp, 84);  // Bytes   84 -- 88: amean
    

    /* TODO: find out how origin is really handled (if meaningful at all) */
    /* In the meantime, we take zero 
     */
    read_float_arr (gf->x0, 3,    fp, 196); // Bytes 196 -- 208: image origin (pixels)
    vec3_mul (gf->x0, gf->csize);
    vec3_set_all (gf->x0, 0.0);
  
    gf->ntotal = idx3_product (gf->shape);
  
    /* Initialize values from data section in the MRC file */
    switch (mode)
      {
        case 1: /* Read 16-bit integers and copy to float */
        i16_arr = (int16_t *) ali16_malloc (gf->ntotal * sizeof (int16_t));
        read_int16_arr (i16_arr, gf->ntotal, fp, MRC_HEADER_BYTES + next);
        
        gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));
        for (i = 0; i < gf->ntotal; i++)
          gf->fvals[i] = i16_arr[i];
  
        free (i16_arr);
        gf->is_initialized = 1;
        gf->is_halfcomplex = 0;
  
        break;
  
        
        case 2: /* Read float values directly */
        gf->fvals = (float *) ali16_malloc (gf->ntotal * sizeof (float));
        read_float_arr (gf->fvals, gf->ntotal, fp, MRC_HEADER_BYTES + next);
        
        gf->is_initialized = 1;
        gf->is_halfcomplex = 0;
        
        break;
  
  
        case 4: /* Read complex float values directly */
        gf->fvals = (float *) ali16_malloc (2 * gf->ntotal * sizeof (float));
        read_float_arr (gf->fvals, 2 * gf->ntotal, fp, MRC_HEADER_BYTES + next);
        
        gf->is_initialized = 1;
        gf->is_halfcomplex = 1;
        
        break;
      
        
        default: /* Never reached since the same check already happened, but anyway.. */
        EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Data mode %d found in %s not supported.", 
          mode, mrc_fname);
      }
    
    if (pfp_in != NULL)
      *pfp_in = fp;
    else
      fclose (fp);
  
  
    gfunc3_compute_xmin_xmax (gf);
  
    if (mode != 4)
      {
        dmin = gfunc3_min (gf);
        dmax = gfunc3_max (gf);
        PRINT_VERBOSE ("Data minimum: %f\n", dmin);
        PRINT_VERBOSE ("Data maximum: %f\n", dmax);
        
        if ((dmin == FLT_MAX) || (dmax == -FLT_MAX))
          fputs ("Warning: Data may be corrupted!\n", stderr);
  
        a = ((dmax - dmin) > EPS_DENOM) ? (amax - amin) / (dmax - dmin) : 1.0;
        b = amax - a * dmax;
  
        gfunc3_scale (gf, a);
        gfunc3_add_constant (gf, b);
      }
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
      
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_read_from_stack (gfunc3 *gf, FILE *fp, int stackpos)
{
  CEXCEPTION_T e = EXC_NONE;
  int16_t *i16_arr;
  int32_t nz, mode, next;
  size_t i, fpos = MRC_HEADER_BYTES;
  float amin, amax, amean, dmin, dmax, a, b;
  
  CAPTURE_NULL (gf);
  CAPTURE_NULL (fp);
  GFUNC_CHECK_INIT_STATUS (gf);
    
  if (stackpos < 0)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "stackpos must be nonnegative.");

  if (!GFUNC_IS_2D(gf))
    EXC_THROW_CUSTOMIZED_PRINT (EXC_GFDIM, "Only stacks of 2d functions supported.");
  
  Try
  {
    /* Re-read shape and origin and re-compute min and max */
    read_int32_arr (gf->shape, 3, fp,  0);  // Bytes    0 -- 12: shape
    gf->shape[2] = 1;
    gf->ntotal = idx3_product (gf->shape);
    // read_float_arr (gf->x0, 3,    fp, 196); // Bytes 196 -- 208: image origin (pixels)
    // vec3_mul (gf->x0, gf->csize);
    gfunc3_compute_xmin_xmax (gf);
    
    read_int32 (&nz  , fp,  8);
    read_int32 (&mode, fp, 12);
    read_int32 (&next, fp, 92);

    if (stackpos >= nz)
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "stackpos (%d) must be smaller than number of images"
        " in stack (%d).", stackpos, nz);

    read_float (&amin,  fp, 76);
    read_float (&amax,  fp, 80);
    read_float (&amean, fp, 84);
    
    switch (mode)
      {
        case 1:  /* Read 16-bit integers and copy to float */
        fpos += next + stackpos * gf->ntotal * sizeof (int16_t);
        i16_arr = (int16_t *) ali16_malloc (gf->ntotal * sizeof (int16_t));
        read_int16_arr (i16_arr, gf->ntotal, fp, fpos);
        
        for (i = 0; i < gf->ntotal; i++)
          gf->fvals[i] = i16_arr[i];

        free (i16_arr);

        break;
          

        case 2:  /* Read float values directly */
        fpos += next + stackpos * gf->ntotal * sizeof (float);
        read_float_arr (gf->fvals, gf->ntotal, fp, fpos);
        
        break;


        case 4:  /* Read complex float values directly */
        fpos += next + stackpos * 2 * gf->ntotal * sizeof (float);
        read_float_arr (gf->fvals, 2 * gf->ntotal, fp, fpos);
        
        break;


        default:
        EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Data mode %d not supported.", mode);

      }

    if (mode != 4)
      {
        dmin = gfunc3_min (gf);
        dmax = gfunc3_max (gf);
        PRINT_VERBOSE ("Data minimum: %f\n", dmin);
        PRINT_VERBOSE ("Data maximum: %f\n", dmax);
        
        if ((dmin == FLT_MAX) || (dmax == -FLT_MAX))
          fputs ("Warning: Data may be corrupted!", stderr);

        a = ((dmax - dmin) > EPS_DENOM) ? (amax - amin) / (dmax - dmin) : 1.0;
        b = amax - a * dmax;

        gfunc3_scale (gf, a);
        gfunc3_add_constant (gf, b);
      }
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_to_mrc (gfunc3 const *gf, char const *mrc_fname)
{
  CEXCEPTION_T e = EXC_NONE;
  const int16_t nreal = 32, endian_stp = 0x4144;
  const int32_t next = 0, nlabl = 1;
  int32_t mode, nx_ny_nz[3];
  float amin = 0.0, amax = 0.0, amean = 0.0;
  char labl[81];
  
  size_t ntotal_flt;
  float fbuf[3];
  FILE *fp;

  CAPTURE_NULL (gf);
  CAPTURE_NULL (mrc_fname);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  if ((fp = fopen (mrc_fname, "w")) == NULL)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to write-only open %s.", mrc_fname);
  
  if (gf->is_halfcomplex)
    {
      mode = 4;
      ntotal_flt = 2 * gf->ntotal;
    }
  else
    {
      mode = 2;
      ntotal_flt = gf->ntotal;
    }

  Try
  {
    /* Write the MRC header */
    
    /*     0 -- 12: nx, ny, nz */
    /* Copy shape from int to int_32_t first */
    nx_ny_nz[0] = gf->shape[0];  nx_ny_nz[1] = gf->shape[1];  nx_ny_nz[2] = gf->shape[2];
    write_int32_arr (nx_ny_nz, 3, fp,  0);

    /*    12 -- 16: mode */
    write_int32     (&mode,   fp, 12);
    
    /*    16 -- 28: (zeros) */
    write_zero_bytes (12, fp, 16);

    /*    28 -- 40: mx, my, mz (set to nx, ny, nz) */
    write_int32_arr (nx_ny_nz, 3, fp,  28);
    
    /*    40 -- 52: total cell size */
    vec3_copy (fbuf, gf->csize);
    vec3_mul_int (fbuf, gf->shape);
    write_float_arr (fbuf, 3,     fp,  40);
    
    /*    52 -- 76: (zeros) */
    write_zero_bytes (24, fp, 52);
    
    /*    76 -- 80: minimum */
    if (mode != 4)
      amin = gfunc3_min (gf);
    write_float     (&amin,       fp,  76);
    
    /*    80 -- 84: maximum */
    if (mode != 4)
      amax = gfunc3_max (gf);
    write_float     (&amax,       fp,  80);
    
    /*    84 -- 88: average */
    if (mode != 4)
      amean = gfunc3_mean (gf);
    write_float     (&amean,      fp,  84);

    /*    88 -- 92: (zeros) */
    write_zero_bytes (4, fp, 88);

    /*    92 -- 96: next */
    write_int32     (&next,       fp,  92);

    /*   96 -- 130: (zeros) */
    write_zero_bytes (34, fp, 96);

    /*  130 -- 132: nreal */
    write_int16     (&nreal,      fp, 130);
    
    /*  132 -- 196: (zeros) */
    write_zero_bytes (64, fp, 132);

    /*  196 -- 208: origin (in grid units) */
    vec3_copy (fbuf, gf->x0);
    vec3_div  (fbuf, gf->csize);
    write_float_arr (fbuf, 3,     fp, 196);
    
    /*  208 -- 212: (zeros) */
    write_zero_bytes (4, fp, 208);
    
    /*  212 -- 216: endianness stamp */
    write_int16     (&endian_stp, fp, 212);
    
    /*  216 -- 220: (zeros) */
    write_zero_bytes (4, fp, 216);

    /*  220 -- 224: nlabl */
    write_int32     (&nlabl,      fp, 220);
    
    /* 224 -- 1024: 10 char[80] labels (unused ones filled with zeros) */
    snprintf (labl, 80, "Created by %s on %s at %s", PACKAGE_STRING, __DATE__, __TIME__);
    write_char_arr  (labl, 80,    fp, 224);
    write_zero_bytes (720, fp, 304);
    

    /* Write the data */
    
    write_float_arr (gf->fvals, ntotal_flt, fp, MRC_HEADER_BYTES);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  fclose (fp);

  return;
}

/*-------------------------------------------------------------------------------------------------*/

#define XSTR(s) STR(s)
#define STR(s) #s
#define CONCAT(s, t, u) s t u
#define XCONCAT(s, t, u) CONCAT (s, t, u)

#define TEMPDIR_STR  "temp/"
#define COUNT_DIGITS  3
#define DIGIT_FMT XCONCAT ("%0", XSTR(COUNT_DIGITS), "d")


void
temp_mrc_out (gfunc3 const *gf, char const *mrc_fbasename, int count)
{
  CAPTURE_NULL (gf);
  CAPTURE_NULL (mrc_fbasename);
  GFUNC_CHECK_INIT_STATUS (gf);
  
  CEXCEPTION_T e = EXC_NONE;
  size_t flen = strlen (mrc_fbasename), templen = strlen (TEMPDIR_STR), extlen = strlen (".mrc");
  char *mrc_fname;


  Try
  {
    /* Try to create the tempdir. If it exists, continue. */
    if (mkdir (TEMPDIR_STR, 0764) && (errno != EEXIST))
      EXC_THROW_CUSTOMIZED_PRINT (EXC_IO, "Unable to create directory '%s'.", TEMPDIR_STR);

    mrc_fname = (char *) ali16_malloc (templen + flen + COUNT_DIGITS + extlen + 1);

    strncpy (mrc_fname, TEMPDIR_STR, templen);
    strncat (mrc_fname, mrc_fbasename, flen);
    
    if (count != 0)
      snprintf (mrc_fname + templen + flen, COUNT_DIGITS + 1, DIGIT_FMT, count);
      
    strncat (mrc_fname, ".mrc", extlen);

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
