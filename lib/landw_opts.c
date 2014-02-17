/*
 * landw_opts.c -- dispatch options for landweber_* programs via getopt
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
#include <getopt.h>
#include <string.h>
#include <ctype.h>

#include "dirname.h"

#include "CException.h"

#include "vec3.h"
#include "misc.h"

#include "et_operators.h"

#include "landw_opts.h"


/*-------------------------------------------------------------------------------------------------*/

#define RELAX_PARAM 1.0F
#define MAX_ITER 25
#define BG_PATCH_SIZE  50 /* Default size of patch to compute background stats */
#define REC_STR   "lw_rec_"  /* prepend this to the output tiltseries file name */

/*-------------------------------------------------------------------------------------------------*/

#define SHORT_OPTS "o:f:p:n:m:s:t:r:M:INhvq"
#define LONG_OPTS \
         /* These options set a flag. */                    \
         {"normalize",        no_argument, &normalize_flag, TRUE}, \
         {"verbose",          no_argument, &verbosity_level, VERB_LEVEL_VERBOSE},  \
         {"quiet",            no_argument, &verbosity_level, VERB_LEVEL_QUIET}, \
         {"invert-contrast",  no_argument, &invert_contrast_flag, 1}, \
         /* These options don't set a flag. We distinguish them by their indices. */\
         {"help",             no_argument, 0, 'h'}, \
         {"version",          no_argument, 0, 'V'}, \
         {"output-file",      required_argument, 0, 'o'}, \
         {"params-file",      required_argument, 0, 'p'}, \
         {"model",            required_argument, 0, 'm'}, \
         {"max-iter",         required_argument, 0, 'M'}, \
         {"background",       required_argument, 0, 'P'}, \
         {"relax-param",      required_argument, 0, 'r'}, \
         {"tiltangles-file",  required_argument, 0, 't'}, \
         {0, 0, 0, 0}


/*-------------------------------------------------------------------------------------------------*/

int verbosity_level      = VERB_LEVEL_NORMAL;
int invert_contrast_flag = TRUE;
int fft_padding          = 0;
int normalize_flag       = TRUE;

char const *models[]      = {"", "proj-assumption", "born-approx", ""};

/*-------------------------------------------------------------------------------------------------*/

LandwOpts *
new_LandwOpts (void)
{
  CEXCEPTION_T e = EXC_NONE;
  LandwOpts *opts = NULL;
  
  Try { opts = (LandwOpts *) ali16_malloc (sizeof (LandwOpts)); } CATCH_RETURN (e, NULL);
  
  opts->fname_in          = NULL;
  opts->fname_out         = NULL;
  opts->fname_params      = NULL;
  opts->fname_tiltangles  = NULL;
  opts->model             = PROJ_ASSUMPTION;
  opts->relax_param       = RELAX_PARAM;
  opts->max_iter          = MAX_ITER;
  opts->bg_patch_ix0[2]   = -1;
  opts->bg_patch_shape[0] = BG_PATCH_SIZE;
  opts->bg_patch_shape[1] = BG_PATCH_SIZE;
  opts->bg_patch_shape[2] = 1;
  
  return opts;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_free (LandwOpts **popts)
{
  if (popts == NULL)
    return;
  
  if ((*popts) == NULL)
    return;
    
  free ((*popts)->fname_in);
  free ((*popts)->fname_out);
  free ((*popts)->fname_params);
  free ((*popts)->fname_tiltangles);

  free (*popts);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_print (LandwOpts const *opts)
{
  CAPTURE_NULL_VOID (opts);
  
  if (verbosity_level == VERB_LEVEL_QUIET)
    return;
  
  printf ("\n\n");
  puts ("Options from command line:");
  puts ("==========================\n");

  printf ("Input file          : %s\n", opts->fname_in);
  printf ("Output file         : %s\n", opts->fname_out);
  printf ("ET parameter file   : %s\n", opts->fname_params);
  printf ("Tiltangles file     : %s\n", opts->fname_tiltangles);
  printf ("Model               : %s\n", models[opts->model]); 

  printf ("Contrast inversion  : " );
  if (invert_contrast_flag)
    printf ("yes\n");
  else
    printf ("no\n");

  printf ("Relaxation parameter: %f\n", opts->relax_param);
  printf ("Maximum iterations  : %d\n", opts->max_iter);

  printf ("Normalization         : " );
  if (normalize_flag)
    printf ("yes\n");
  else
    printf ("no\n");
  
  printf ("Background patch      : ");
  if (!normalize_flag)
    printf ("(unused, no normalization)\n");

  else
    printf ("size %dx%d at (%d,%d)\n", opts->bg_patch_shape[0], opts->bg_patch_shape[1], 
    opts->bg_patch_ix0[0], opts->bg_patch_ix0[1]);

  printf ("\n\n");
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
print_help (char const *progname)
{
  scattering_model iter_m;
  
  printf ("Usage: %s [options] TILTSERIES_FILE\n", progname);
  puts ("");
  puts ("Required parameters:");
  puts ("");
  puts ("  -t FILE, --tiltangles-file=FILE");
  puts ("                 read tilt angles from FILE.");
  puts ("  -m NAME, --model=NAME");
  puts ("                 set forward model to NAME. Possible values are:");
  printf ("               ");
  for (iter_m = MD_START + 1; iter_m < MD_END; iter_m++)
    printf ("  %s", models[iter_m]);
  printf ("\n");
  puts ("                 (Default: proj-assumption)");
  puts ("  -p FILE, --params-file=FILE");
  puts ("                 read CTF parameters from FILE (INI-style).");
  puts ("");
  puts ("Options:");
  puts ("");
  puts ("  -o FILE, --output-file=FILE");
  puts ("                 write reconstruction to FILE; if no parameter is given, the");
  puts ("                 output file is determined from TILTSERIES_FILE by prepending");
  printf ("                 `%s'.\n", REC_STR);
  puts ("  -N, --normalize");
  puts ("                 normalize the projection images based on their histograms (enabled");
  puts ("                 by default).");
  puts ("  --background=index_x,index_y,size_x,size_y");
  puts ("                 compute background statistics using a patch of SIZE_X x SIZE_Y");
  puts ("                 pixels with lower-left indices INDEX_X, INDEX_Y.");
  printf ("                 Default: 0,0,%d,%d\n", 
    BG_PATCH_SIZE, BG_PATCH_SIZE);
  puts ("  -I, --invert-contrast");
  puts ("                 invert the contrast of the images; use this option if projections of");
  puts ("                 dense regions are brighter than the background (enabled by default).");
  puts ("  -r MU, --relax-param=MU");
  printf ("                 use MU as relaxation parameter in the iteration (default: %.2f)\n", 
    RELAX_PARAM);
  puts ("  -M NUM, --max-iter=NUM");
  printf ("                 do at most NUM iteration steps (default: %d)\n", MAX_ITER);
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

void
LandwOpts_set_fname_in (LandwOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_in = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_in, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_determine_fname_out (LandwOpts *opts, char const *fname_in)
{
  CEXCEPTION_T e = EXC_NONE;
  int dir_len, base_len, rec_len;
  char *dirname, *basename, *p_tmp;
  
  dirname  = dir_name (fname_in);
  dir_len  = strlen (dirname);
  basename = base_name (fname_in);
  base_len = strlen (basename);
  rec_len   = strlen (REC_STR);
  
  Try { 
    opts->fname_out = (char *) ali16_malloc (dir_len + base_len + rec_len + 2); 
  } CATCH_RETURN_VOID (e);
    
  p_tmp = opts->fname_out;
  strncpy (p_tmp, dirname, dir_len);
  
  p_tmp += dir_len;
  *(p_tmp++) = '/';
  strncpy (p_tmp, REC_STR, rec_len);
  
  p_tmp += rec_len;
  strncpy (p_tmp, basename, base_len);
  
  p_tmp += base_len;
  *p_tmp = '\0';
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_fname_out (LandwOpts *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
 
  Try { od->fname_out = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (od->fname_out, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_fname_params (LandwOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);

  Try { opts->fname_params = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_params, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_fname_tiltangles (LandwOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_tiltangles = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_tiltangles, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_model (LandwOpts *opts, char *model_str)
{
  char *p;
  scattering_model iter;
  
  for (p = model_str; *p; p++) 
    *p = tolower (*p);
  
  for (iter = MD_START + 1; iter < MD_END; iter++)
    {
      if (strcmp (model_str, models[iter]) == 0)
        {
          opts->model = iter;
          return;
        }
    }
    
  EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Unknown model `%s'\n", model_str);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_relax_param (LandwOpts *opts, float tau)
{
  if (tau < 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Relaxation parameter must be nonnegative.");
      return;
    }
    
  opts->relax_param = tau;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_max_iter (LandwOpts *opts, int iter)
{
  if (iter <= 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Number of iterations must be positive.");
      return;
    }
    
  opts->max_iter = iter;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_set_bg_params (LandwOpts *opts, char const *params_str)
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
  
  opts->bg_patch_ix0[0]   = ix;
  opts->bg_patch_ix0[1]   = iy;
  opts->bg_patch_ix0[2]   = iz;
  opts->bg_patch_shape[0] = sx;
  opts->bg_patch_shape[1] = sy;
  opts->bg_patch_shape[2] = 1;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
LandwOpts_assign_from_args (LandwOpts *opts, int argc, char **argv)
{
  CAPTURE_NULL_VOID (opts);

  /* Aux variables */
  int c;
  char const *progname = base_name(argv[0]);
  
  /* Internal flags for the short options */
  int m_flag = 0;
  int M_flag = 0;
  int o_flag = 0;
  int p_flag = 0;
  int P_flag = 0;
  int r_flag = 0;
  int t_flag = 0;
 
  if (argc == 1)
    {
      print_version_etc (progname);
      puts ("\n");
      print_help (progname);
      exit (EXIT_SUCCESS);
    }
 
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
          /* This case should never occur! */
          printf ("This should never occur! Option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'I':
          invert_contrast_flag = 1;  /* Long option already sets the flag, but short one doesn't */
          break;

        case 'm':
          if (m_flag != 0)
            {
              fputs ("Invalid multiple use of `-m' (`--model') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_model (opts, optarg);
          m_flag = 1;
          break;
 
        case 'M':
          if (M_flag != 0)
            {
              fputs ("Invalid multiple use of `-M' (`--max-iter') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_max_iter (opts, atoi(optarg));
          M_flag = 1;
          break;
 
        case 'N':
          normalize_flag = TRUE;  /* Long option already sets the flag, but short one doesn't */
          break;
 
        case 'o':
          if (o_flag != 0)
            {
              fputs ("Invalid multiple use of `-o' (`--output-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_fname_out (opts, optarg);
          o_flag = 1;
          break;
 
        case 'p':
          if (p_flag != 0)
            {
              fputs ("Invalid multiple use of `-p' (`--et-params-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_fname_params (opts, optarg);
          p_flag = 1;
          break;
 
        case 'P':
          if (P_flag != 0)
            {
              fputs ("Invalid multiple use of `--background' option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_bg_params (opts, optarg);
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
 
        case 'r':
          if (r_flag != 0)
            {
              fputs ("Invalid multiple use of `-r' (`--relax-param') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_relax_param (opts, atof (optarg));
          r_flag = 1;
          break;
 
        case 't':
          if (t_flag != 0)
            {
              fputs ("Invalid multiple use of `-t' (`--tiltangles-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          LandwOpts_set_fname_tiltangles (opts, optarg);
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
          fprintf(stderr, "`%s --help' provides further information.\n\n", progname);
          exit (EXIT_FAILURE);

        case 'h':

        default:
          print_help (progname);
          exit (EXIT_FAILURE);
        }
    }

  /* Handle various option conflicts */
  
  if (!m_flag)
    {
      fprintf (stderr, "Error: `--model` option missing.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }

  if (!t_flag)
    {
      fprintf (stderr, "Error: `--tiltangles-file` option missing.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }

  if (!p_flag)
    {
      fprintf (stderr, "Error: `--et-params-file` option missing.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }

  if ((!normalize_flag) && P_flag)
    {
      fprintf (stderr, "The `--backgound' option can only be used if the "
        "`-N' (`--normalize') option is enabled.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }


  /* Process any remaining command line arguments (not options). */
  if ((argc - optind) != 1 )
    {
      print_help (progname);
      exit (EXIT_FAILURE);
    }

  LandwOpts_set_fname_in (opts, argv[optind]);

  
  if (!o_flag)
    LandwOpts_determine_fname_out (opts, argv[optind]);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
