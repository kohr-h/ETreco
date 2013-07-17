/*
 * options.c -- dispatch options via getopt
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
#include <getopt.h>
#include <string.h>
#include <ctype.h>

#if GNULIB_DIRNAME
#include "dirname.h"
#endif

#include "CException.h"
#include "misc.h"
#include "options.h"
#include "matvec3.h"

#define REC_STR   "_rec"  /* append this before the file extension for the reco */

#define BG_PATCH_SIZE  50 /* Default size of patch to compute background stats */
#define FFT_PADDING    64 /* Default padding of functions before FFT */

/*-------------------------------------------------------------------------------------------------*/

#define SHORT_OPTS "o:g:c:f:p:n:m:s:t:NILhvq"
#define LONG_OPTS \
         /* These options set a flag. */                    \
         {"verbose",          no_argument, &verbosity_level, VERB_LEVEL_VERBOSE},  \
         {"quiet",            no_argument, &verbosity_level, VERB_LEVEL_QUIET}, \
         {"normalize",        no_argument, &normalize_flag, 1}, \
         {"invert-contrast",  no_argument, &invert_contrast_flag, 1}, \
         /* These options don't set a flag. We distinguish them by their indices. */\
         {"help",             no_argument, 0, 'h'}, \
         {"version",          no_argument, 0, 'V'}, \
         {"output-file",      required_argument, 0, 'o'}, \
         {"gamma",            required_argument, 0, 'g'}, \
         {"ctf-cutoff",       required_argument, 0, 'c'}, \
         {"reco-params-file", required_argument, 0, 'p'}, \
         {"num-images",       required_argument, 0, 'n'}, \
         {"start-index",      required_argument, 0, 's'}, \
         {"tiltangles-file",  required_argument, 0, 't'}, \
         {"mollifier",        required_argument, 0, 'm'}, \
         {"fft-padding",      required_argument, 0, 'F'}, \
         {"background",       required_argument, 0, 'P'}, \
         {"lambda-pow",       required_argument, 0, 'L'}, \
         {0, 0, 0, 0}


/*-------------------------------------------------------------------------------------------------*/

int verbosity_level      = VERB_LEVEL_NORMAL;
int use_gamma_flag       = 0;
int truncate_ctf_flag    = 0;
int normalize_flag       = 0;
int invert_contrast_flag = 0;
int use_lambda_flag      = 0;
int fft_padding          = 0;

char const *mollifiers[] = {"delta", "gaussian"};

/*-------------------------------------------------------------------------------------------------*/

OptionData *
new_OptionData (void)
{
  CEXCEPTION_T e = EXC_NONE;
  OptionData *od;
  
  Try
  {
    od = (OptionData *) ali16_malloc (sizeof (OptionData));
    
    od->fname_in          = NULL;
    od->fname_out         = NULL;
    od->fname_reco_params = NULL;
    od->fname_tiltangles  = NULL;
    od->gamma             = 0.0;
    od->ctf_trunc         = 0.0;
    od->moll_type         = DELTA;
    od->num_images        = 0;
    od->start_index       = 0;
    od->lambda_pow        = 0.0;
    od->bg_patch_ix0[2]   = -1;
    od->bg_patch_shape[0] = BG_PATCH_SIZE;
    od->bg_patch_shape[1] = BG_PATCH_SIZE;
    od->bg_patch_shape[2] = 1;
    
    return od;
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return NULL;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_free (OptionData **pod)
{
  if (pod == NULL)
    return;
  
  if ((*pod) == NULL)
    return;
    
  free ((*pod)->fname_in);
  free ((*pod)->fname_out);
  free ((*pod)->fname_reco_params);
  free ((*pod)->fname_tiltangles);

  free (*pod);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_print (OptionData *od)
{
  CAPTURE_NULL (od);
  
  /* TODO: make dependent on verbosity */
  printf ("\n\n");
  puts ("Options from command line:");
  puts ("==========================\n");

  printf ("Input file       : %s\n", od->fname_in);
  printf ("Output file      : %s\n", od->fname_out);
  printf ("Reco params file : %s\n", od->fname_reco_params);
  printf ("Tiltangles file  : %s\n", od->fname_tiltangles);

  printf ("gamma            : ");
  if (use_gamma_flag)
    printf ("%e\n", od->gamma);
  else
    printf ("(guessed from data)\n");

  printf ("1/CTF cutoff     : ");
  if (truncate_ctf_flag)
    printf ("%e\n", od->ctf_trunc);
  else
    printf ("(ignored)\n");
    
  printf ("Mollifier        : %s\n", mollifiers[od->moll_type]); 

  printf ("Normalization    : " );
  if (normalize_flag)
    printf ("yes\n");
  else
    printf ("no\n");
  
  printf ("Background patch : ");
  if (!normalize_flag)
    printf ("(unused, no normalization)\n");
  // else if (od->bg_patch_ix0[2] == -1)
    // printf ("size %dx%d (position guessed)\n", od->bg_patch_shape[0], od->bg_patch_shape[1]);
  else
    printf ("size %dx%d at (%d,%d)\n", od->bg_patch_shape[0], od->bg_patch_shape[1], 
    od->bg_patch_ix0[0], od->bg_patch_ix0[1]);

  printf ("Number of images : ");
  if (od->num_images == -1)
    printf ("(determined from data)\n");
  else
    printf ("%d\n", od->num_images);

  printf ("Start index      : %d\n", od->start_index);

  printf ("Lambda           : ");
  if (use_lambda_flag)
    printf ("Lambda^a with a=%f\n", od->lambda_pow);
  else  
    printf ("(not used)\n");
  
  printf ("FFT zero-padding : %d\n", fft_padding);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
print_help (char const *progname)
{
  mollifier_type iter;
  
  printf ("Usage: %s -t <tiltangles-file> [options] tiltseries_file\n", progname);
  puts ("Options:");
  
  puts ("  -t file, --tiltangles-file=file");
  puts ("                 read tilt angles from this text file. This option is required.");
  puts ("  -p file, --reco-params-file=file");
  puts ("                 read CTF and reconstruction parameters from this ini-style file.");
  puts ("                 This option is required.");
  puts ("  -o file, --output-file=file");
  puts ("                 write reconstruction to the specified file; if not given,");
  puts ("                 output file is determined from tiltseries_file by appending");
  puts ("                 `_rec' before its extension.");
  puts ("  -g value, --gamma=value");
  puts ("                 set value for regularization parameter gamma; magnitude");
  puts ("                 should be in the order of detector pixel size [microns].");
  puts ("  -c value, --ctf-cutoff=value");
  puts ("                 set value for reciprocal CTF cutoff; in intervals where ");
  puts ("                 1/CTF exceeds this value, it is replaced by a diff'able");
  puts ("                 transition spline; must be > 1.0 and should not exceed ~5.0.");
  puts ("  -m name, --mollifier=name");
  puts ("                 use the specified mollifier instead of delta; currently ");
  printf ("                 supported are:");
  for (iter = DELTA; iter < LAST; iter++)
    printf ("  %s", mollifiers[iter]);
  printf ("\n");
  puts ("  -N, --normalize");
  puts ("                 normalize the projection images based on their histograms.");
  puts ("  --background=index_x,index_y,size_x,size_y");
  puts ("                 compute background statistics using a patch of size_x x size_y");
  puts ("                 pixels with lower-left indices index_x, index_y.");
  printf ("                 Default is 0,0,%d,%d in the case that this option is omitted.\n", 
    BG_PATCH_SIZE, BG_PATCH_SIZE);
  puts ("  -I, --invert-contrast");
  puts ("                 invert the contrast during normalization; use this option if");
  puts ("                 objects are darker than the background.");
  puts ("  -L value, --lambda-pow=value");
  puts ("                 feature reconstruction: compute Lambda^a(f) instead of f,");
  puts ("                 where Lambda = sqrt(-Laplacian) and a = value.");
  puts ("  -n N, --num-images=N");
  puts ("                 override image number determined from input file by N; useful");
  puts ("                 for recostruction with only a subset of the data.");
  puts ("  -s N, --start-index=N");
  puts ("                 start with N'th image instead of 0; useful for");
  puts ("                 recostruction with only a subset of the data.");
  puts ("  --fft-padding[=N]");
  puts ("                 continue grid functions by N zero pixels in each direction");
  puts ("                 to compute Fourier transforms ('zero-padding').");
  printf ("                 Default: N=%d\n", FFT_PADDING);
  puts ("  -v, --verbose");
  puts ("                 display more information during execution");
  puts ("  -q, --quiet");
  puts ("                 display less information during execution");
  puts ("  -h, --help");
  puts ("                 display this help and exit");
  puts ("      --version");
  puts ("                 output version information and exit");
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

enum { COPYRIGHT_YEAR = 2013 };

void 
print_version_etc (char const *progname)
{
  printf ("%s (%s) %s\n", progname, PACKAGE_NAME, PACKAGE_VERSION);
  printf ("Copyright (C) %d %s\n", COPYRIGHT_YEAR, AUTHORS);
  
  puts ("\
\n\
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.");
  
  printf ("\nReport bugs to: %s\n", PACKAGE_BUGREPORT);
  #ifdef PACKAGE_URL
  printf ("%s home page: <%s>\n", PACKAGE_NAME, PACKAGE_URL);
  #endif

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
fname_rsplit_at_dot (char const *fname, char **pbase_str, char **pext_str)
{
  CEXCEPTION_T e = EXC_NONE;
  int fname_len, base_len, ext_len;
  char const *pext = NULL;
  
  fname_len = strlen (fname);

  /* pext points to the first dot from the right or one char beyond the string */
  if ((pext = strrchr (fname, '.')) == NULL)
    pext = &fname[fname_len];
    
  base_len = pext - fname;
  ext_len  = fname_len - base_len;
  
  Try
  {
    *pbase_str = (char *) ali16_malloc (base_len + 1);
    strncpy (*pbase_str, fname, base_len);  (*pbase_str)[base_len] = '\0';
    
    *pext_str = (char *) ali16_malloc (ext_len + 1);
    strncpy (*pext_str, pext, ext_len);  (*pext_str)[ext_len] = '\0';
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }

  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_fname_in (OptionData *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try
  {
    od->fname_in = (char *) ali16_malloc (len + 1);
    strncpy (od->fname_in, fname, len);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }  
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_assemble_fname_out (OptionData *od,  char const *base, char const *ext)
{
  CEXCEPTION_T e = EXC_NONE;
  int base_len, rec_len, ext_len;
  char *p_tmp;
  
  base_len = strlen (base);
  rec_len  = strlen (REC_STR);
  ext_len  = strlen (ext);
  
  Try
  {
    od->fname_out = (char *) ali16_malloc (base_len + rec_len + ext_len + 1);
    
    p_tmp = od->fname_out;
    strncpy (p_tmp, base, base_len);
    
    p_tmp += base_len;
    strncpy (p_tmp, REC_STR, rec_len);
    
    p_tmp += rec_len;
    strncpy (p_tmp, ext, ext_len);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_determine_fname_out (OptionData *od, char *fname_in)
{
  CEXCEPTION_T e = EXC_NONE;
  char *bs, *ex;
  
  Try
  {
    fname_rsplit_at_dot (fname_in, &bs, &ex);
    OptionData_assemble_fname_out (od, bs, ex);
    
    free (bs);
    free (ex);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }  
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_fname_out (OptionData *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
 
  Try
  {
    od->fname_out = (char *) ali16_malloc (len + 1);
    strncpy (od->fname_out, fname, len);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_fname_reco_params (OptionData *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);

  Try
  {
    od->fname_reco_params = (char *) ali16_malloc (len + 1);
    strncpy (od->fname_reco_params, fname, len);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_fname_tiltangles (OptionData *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try
  {
    od->fname_tiltangles = (char *) ali16_malloc (len + 1);
    strncpy (od->fname_tiltangles, fname, len);
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

/* TODO: change fputs..exit to exceptions */
void
OptionData_set_gamma (OptionData *od, float gamma)
{
  if (gamma <= 0.0)
    {
      fputs ("Parameter of `-g' (`--gamma') must be positive.", stderr);
      exit (EXIT_FAILURE);
    }
  
  od->gamma = gamma;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_ctf_trunc (OptionData *od, float ctf_cut)
{
  if (ctf_cut <= 0.0)
    {
      fputs ("Parameter of `-c' (`--ctf-cutoff') must be positive.", stderr);
      exit (EXIT_FAILURE);
    }
  
  od->ctf_trunc = ctf_cut;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_num_images (OptionData *od, int n_images)
{
  if (n_images <= 0)
    {
      fputs ("Parameter of `-n' (`--num-images') must be positive.", stderr);
      exit (EXIT_FAILURE);
    }
  
  od->num_images = n_images;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_start_index (OptionData *od, int nstart)
{
  if (nstart < 0)
    {
      fputs ("Parameter of `-s' (`--start-index') must be nonnegative.", stderr);
      exit (EXIT_FAILURE);
    }
  
  od->start_index = nstart;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_moll_type (OptionData *od, char *moll_str)
{
  char *p;
  mollifier_type iter;
  
  for (p = moll_str; *p; p++) 
    *p = tolower (*p);
  
  for (iter = DELTA; iter < LAST; iter++)
    {
      if (strcmp (moll_str, mollifiers[iter]) == 0)
        {
          od->moll_type = iter;
          return;
        }
    }
    
  fprintf (stderr, "Unknown mollifier `%s'\n", moll_str);
  exit (EXIT_FAILURE);
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_bg_params (OptionData *od, char const *params_str)
{
  int ix, iy, iz = 0, sx, sy;
  int n = sscanf (params_str, "%d,%d,%d,%d", &ix, &iy, &sx, &sy);
  
  if (n != 4)  /* Nothing provided or wrong format -> defaults */
    {
      ix = 0;
      iy = 0;
      iz = -1;
      sx = BG_PATCH_SIZE;
      sy = BG_PATCH_SIZE;
    }
  
  od->bg_patch_ix0[0]   = ix;
  od->bg_patch_ix0[1]   = iy;
  od->bg_patch_ix0[2]   = iz;
  od->bg_patch_shape[0] = sx;
  od->bg_patch_shape[1] = sy;
  od->bg_patch_shape[2] = 1;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_set_lambda_pow (OptionData *od, float l_a)
{
  od->lambda_pow = l_a;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
OptionData_assign_from_args (OptionData *od, int argc, char **argv)
{
  CAPTURE_NULL (od);

  /* Aux variables */
  CEXCEPTION_T e = EXC_NONE;
  int c;
  char const *progname = base_name(argv[0]);
  float gamma;
  float ctf_cut;
  int n_images;
  int nstart;
  float l_a;
  
  /* Internal flags for the short options */
  static int o_flag = 0;
  static int n_flag = 0;
  static int s_flag = 0;
  static int t_flag = 0;
  static int m_flag = 0;
  static int P_flag = 0;
  static int p_flag = 0;
  static int F_flag = 0;
 
  if (argc == 1)
    {
      print_version_etc (progname);
      puts ("\n");
      print_help (progname);
      exit (EXIT_SUCCESS);
    }
 
  Try
  {
    while (1)
      {
        static struct option long_options[] = { LONG_OPTS };

        /* getopt_long stores the option index here. */
        int option_index = 0;
   
        c = getopt_long (argc, argv, SHORT_OPTS, long_options, &option_index);
   
        /* Detect the end of the options. */
        if (c == -1)
          break;
   
        switch (c)
          {
          case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0) 
              break;
            // This case should never occur!
            printf ("This should never occur! Option %s", long_options[option_index].name);
            if (optarg)
              printf (" with arg %s", optarg);
            printf ("\n");
            break;

          case 'c':
            if (truncate_ctf_flag != 0)
              {
                fputs ("Invalid multiple use of `-c' (`--ctf-cutoff') option.", stderr);
                exit (EXIT_FAILURE);
              }

            ctf_cut = atof (optarg);
            OptionData_set_ctf_trunc (od, ctf_cut);
            truncate_ctf_flag = 1;
            break;

          case 'F':
            if (F_flag != 0)
              {
                fputs ("Invalid multiple use of `--fft-padding' option.", stderr);
                exit (EXIT_FAILURE);
              }
            
            if (optarg)
              {
                fft_padding = atoi (optarg);
                if (fft_padding < 0)
                  {
                    fputs ("Parameter of `--fft-padding' must be nonnegative.", stderr);
                    exit (EXIT_FAILURE);
                  }
              }
            else
              fft_padding = FFT_PADDING;
              
            F_flag = 1;
            break;
   
          case 'g':
            if (use_gamma_flag != 0)
              {
                fputs ("Invalid multiple use of `-g' (`--gamma') option.", stderr);
                exit (EXIT_FAILURE);
              }

            gamma = atof (optarg);
            OptionData_set_gamma (od, gamma);
            use_gamma_flag = 1;
            break;
   
          case 'I':
            /* invert_contrast_flag = 1; */
            break;

          case 'L':
            if (use_lambda_flag != 0)
              {
                fputs ("Invalid multiple use of `-L' (`--lambda-pow') option.", stderr);
                exit (EXIT_FAILURE);
              }

            l_a = atof (optarg);
            OptionData_set_lambda_pow (od, l_a);
            use_lambda_flag = 1;
            break;
   
          case 'm':
            if (m_flag != 0)
              {
                fputs ("Invalid multiple use of `-m' (`--mollifier') option.", stderr);
                exit (EXIT_FAILURE);
              }

            OptionData_set_moll_type (od, optarg);
            m_flag = 1;
            break;
   
          case 'n':
            if (n_flag != 0)
              {
                fputs ("Invalid multiple use of `-n' (`--num-images') option.", stderr);
                exit (EXIT_FAILURE);
              }

            n_images = atoi (optarg);
            OptionData_set_num_images (od, n_images);
            n_flag = 1;
            break;

          case 'N':
            /* normalize_flag = 1; */
            break;
   
          case 'o':
            if (o_flag != 0)
              {
                fputs ("Invalid multiple use of `-o' (`--output-file') option.", stderr);
                exit (EXIT_FAILURE);
              }

            OptionData_set_fname_out (od, optarg);
            o_flag = 1;
            break;
   
          case 'p':
            if (p_flag != 0)
              {
                fputs ("Invalid multiple use of `-p' (`--reco-params-file') option.", stderr);
                exit (EXIT_FAILURE);
              }

            OptionData_set_fname_reco_params (od, optarg);
            p_flag = 1;
            break;
   
          case 'P':
            if (P_flag != 0)
              {
                fputs ("Invalid multiple use of `--background' option.", stderr);
                exit (EXIT_FAILURE);
              }

            OptionData_set_bg_params (od, optarg);
            P_flag = 1;
            break;
   
          case 'q':
            if (verbosity_level != VERB_LEVEL_NORMAL)
              {
                fputs ("Invalid multiple use of `-q' (`--quiet') or `-v' (`--verbose') options.", 
                  stderr);
                exit (EXIT_FAILURE);
              }
              
            verbosity_level = VERB_LEVEL_QUIET;
            break;
   
          case 's':
            if (s_flag != 0)
              {
                fputs ("Invalid multiple use of `-s' (`--start-index') option.", stderr);
                exit (EXIT_FAILURE);
              }

            nstart = atoi (optarg);
            OptionData_set_start_index (od, nstart);
            s_flag = 1;
            break;
   
          case 't':
            if (t_flag != 0)
              {
                fputs ("Invalid multiple use of `-t' (`--tiltangles-file') option.", stderr);
                exit (EXIT_FAILURE);
              }

            OptionData_set_fname_tiltangles (od, optarg);
            t_flag = 1;
            break;
   
          case 'v':
            if (verbosity_level != VERB_LEVEL_NORMAL)
              {
                fputs ("Invalid multiple use of `-q' (`--quiet') or `-v' (`--verbose') options.", 
                  stderr);
                exit (EXIT_FAILURE);
              }
              
            verbosity_level = VERB_LEVEL_VERBOSE;
            break;
   
          case 'V':
            print_version_etc (progname);
            exit (EXIT_SUCCESS);
            break;
   
          case '?':
            /* getopt_long printed an error message. */
            fprintf(stderr, "`%s --help' provides further information.\n", progname);
            exit (EXIT_FAILURE);

          case 'h':

          default:
            print_help (progname);
            exit (EXIT_FAILURE);
          }
      }

    /* Handle various option conflicts */
    
    if (!t_flag)
      {
        fprintf (stderr, "No tiltangles file given.\n\n" 
          "`%s --help' provides further information.", progname);
        exit (EXIT_FAILURE);
      }

    if (!p_flag)
      {
        fprintf (stderr, "No reco parameters file given.\n\n" 
          "`%s --help' provides further information.", progname);
        exit (EXIT_FAILURE);
      }

    if ((!normalize_flag) && P_flag)
      {
        fprintf (stderr, "The `--backgound' option can only be used if the "
          "`-N' (`--normalize') option is enabled.\n\n" 
          "`%s --help' provides further information.", progname);
        exit (EXIT_FAILURE);
      }

      
    /* Process any remaining command line arguments (not options). */
    if ((argc - optind) != 1 )
      {
        print_help (progname);
        exit (EXIT_FAILURE);
      }

    OptionData_set_fname_in (od, argv[optind]);

    
    if (o_flag == 0)
      OptionData_determine_fname_out (od, argv[optind]);
    
  }
  Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
  }
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/
