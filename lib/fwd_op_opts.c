/*
 * fwd_opts.c -- dispatch options for forward_op_* programs via getopt
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

#include "fwd_op_opts.h"


#define TSERIES_STR   "tiltseries_"  /* prepend this to the output tiltseries file name */

/*-------------------------------------------------------------------------------------------------*/

#define SHORT_OPTS "o:f:p:n:m:s:t:Ihvq"
#define LONG_OPTS \
         /* These options set a flag. */                    \
         {"verbose",          no_argument, &verbosity_level, VERB_LEVEL_VERBOSE},  \
         {"quiet",            no_argument, &verbosity_level, VERB_LEVEL_QUIET}, \
         {"invert-contrast",  no_argument, &invert_contrast_flag, 1}, \
         /* These options don't set a flag. We distinguish them by their indices. */\
         {"help",             no_argument, 0, 'h'}, \
         {"version",          no_argument, 0, 'V'}, \
         {"output-file",      required_argument, 0, 'o'}, \
         {"params-file",      required_argument, 0, 'p'}, \
         {"model",            required_argument, 0, 'm'}, \
         {"tiltangles-file",  required_argument, 0, 't'}, \
         {0, 0, 0, 0}


/*-------------------------------------------------------------------------------------------------*/

int verbosity_level      = VERB_LEVEL_NORMAL;
int invert_contrast_flag = 1;
int fft_padding          = 0;

char const *models[]      = {"", "proj-assumption", "born-approx", ""};

/*-------------------------------------------------------------------------------------------------*/

FwdOpts *
new_FwdOpts (void)
{
  CEXCEPTION_T e = EXC_NONE;
  FwdOpts *opts = NULL;
  
  Try { opts = (FwdOpts *) ali16_malloc (sizeof (FwdOpts)); } CATCH_RETURN (e, NULL);
  
  opts->fname_in                = NULL;
  opts->fname_out               = NULL;
  opts->fname_params            = NULL;
  opts->fname_tiltangles        = NULL;
  opts->model                   = PROJ_ASSUMPTION;
  
  return opts;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_free (FwdOpts **popts)
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
FwdOpts_print (FwdOpts *opts)
{
  CAPTURE_NULL_VOID (opts);
  
  if (verbosity_level == VERB_LEVEL_QUIET)
    return;
  
  printf ("\n\n");
  puts ("Options from command line:");
  puts ("==========================\n");

  printf ("Input file            : %s\n", opts->fname_in);
  printf ("Output file           : %s\n", opts->fname_out);
  printf ("ET parameter file     : %s\n", opts->fname_params);
  printf ("Tiltangles file       : %s\n", opts->fname_tiltangles);
  printf ("Model                 : %s\n", models[opts->model]); 

  printf ("Contrast inversion    : " );
  if (invert_contrast_flag)
    printf ("yes\n");
  else
    printf ("no\n");

  printf ("\n\n");
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
print_help (char const *progname)
{
  scattering_model iter_m;
  
  printf ("Usage: %s [options] volume_file\n", progname);
  puts ("");
  puts ("Required parameters:");
  puts ("");
  puts ("  -t file, --tiltangles-file=file");
  puts ("                 read tilt angles from FILE.");
  puts ("  -m name, --model=name");
  puts ("                 set forward model to NAME. Possible values are:");
  printf ("               ");
  for (iter_m = MD_START + 1; iter_m < MD_END; iter_m++)
    printf ("  %s", models[iter_m]);
  printf ("\n");
  puts ("                 (Default: proj-assumption)");
  puts ("  -p file, --params-file=file");
  puts ("                 read CTF parameters from FILE (INI-style).");
  puts ("");
  puts ("Options:");
  puts ("");
  puts ("  -o file, --output-file=file");
  puts ("                 write reconstruction to FILE; if no parameter is given, the");
  puts ("                 output file is determined from VOLUME_FILE by prepending");
  printf ("                 `%s'.\n", TSERIES_STR);
  puts ("  -I, --invert-contrast");
  puts ("                 invert the contrast of the images; use this option if projections of");
  puts ("                 dense regions are brighter than the background (enabled by default).");
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
FwdOpts_set_fname_in (FwdOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_in = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_in, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_determine_fname_out (FwdOpts *opts, char const *fname_in)
{
  CEXCEPTION_T e = EXC_NONE;
  int dir_len, base_len, ts_len;
  char *dirname, *basename, *p_tmp;
  
  dirname  = dir_name (fname_in);
  dir_len  = strlen (dirname);
  basename = base_name (fname_in);
  base_len = strlen (basename);
  ts_len   = strlen (TSERIES_STR);
  
  Try { 
    opts->fname_out = (char *) ali16_malloc (dir_len + base_len + ts_len + 2); 
  } CATCH_RETURN_VOID (e);
    
  p_tmp = opts->fname_out;
  strncpy (p_tmp, dirname, dir_len);
  
  p_tmp += dir_len;
  *(p_tmp++) = '/';
  strncpy (p_tmp, TSERIES_STR, ts_len);
  
  p_tmp += ts_len;
  strncpy (p_tmp, basename, base_len);
  
  p_tmp += base_len;
  *p_tmp = '\0';
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_set_fname_out (FwdOpts *od, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
 
  Try { od->fname_out = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (od->fname_out, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_set_fname_params (FwdOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);

  Try { opts->fname_params = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_params, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_set_fname_tiltangles (FwdOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_tiltangles = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_tiltangles, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
FwdOpts_set_model (FwdOpts *opts, char *model_str)
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
FwdOpts_assign_from_args (FwdOpts *opts, int argc, char **argv)
{
  CAPTURE_NULL_VOID (opts);

  /* Aux variables */
  int c;
  char const *progname = base_name(argv[0]);
  
  /* Internal flags for the short options */
  int m_flag = 0;
  int o_flag = 0;
  int p_flag = 0;
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

          FwdOpts_set_model (opts, optarg);
          m_flag = 1;
          break;
 
        case 'o':
          if (o_flag != 0)
            {
              fputs ("Invalid multiple use of `-o' (`--output-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          FwdOpts_set_fname_out (opts, optarg);
          o_flag = 1;
          break;
 
        case 'p':
          if (p_flag != 0)
            {
              fputs ("Invalid multiple use of `-p' (`--et-params-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          FwdOpts_set_fname_params (opts, optarg);
          p_flag = 1;
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
 
        case 't':
          if (t_flag != 0)
            {
              fputs ("Invalid multiple use of `-t' (`--tiltangles-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          FwdOpts_set_fname_tiltangles (opts, optarg);
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

    
  /* Process any remaining command line arguments (not options). */
  if ((argc - optind) != 1 )
    {
      print_help (progname);
      exit (EXIT_FAILURE);
    }

  FwdOpts_set_fname_in (opts, argv[optind]);

  
  if (!o_flag)
    FwdOpts_determine_fname_out (opts, argv[optind]);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
