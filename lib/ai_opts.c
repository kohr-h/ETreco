/*
 * ai_opts.c -- dispatch options for ai_* programs via getopt
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

#include "CException.h"

#include "vec3.h"
#include "misc.h"

#include "ai_opts.h"


#define BG_PATCH_SIZE  50 /* Default size of patch to compute background stats */
#define FFT_PADDING    64 /* Default padding of functions before FFT */

/*-------------------------------------------------------------------------------------------------*/

#define SHORT_OPTS "o:g:c:f:p:n:m:s:t:ANILhvq"
#define LONG_OPTS \
         /* These options set a flag. */                    \
         {"verbose",          no_argument, &verbosity_level, VERB_LEVEL_VERBOSE},  \
         {"quiet",            no_argument, &verbosity_level, VERB_LEVEL_QUIET}, \
         {"normalize",        no_argument, &normalize_flag, 1}, \
         {"invert-contrast",  no_argument, &invert_contrast_flag, 1}, \
         {"autocenter-volume",no_argument, &autocenter_vol_flag, 1}, \
         /* These options don't set a flag. We distinguish them by their indices. */\
         {"help",             no_argument, 0, 'h'}, \
         {"version",          no_argument, 0, 'V'}, \
         {"output-file",      required_argument, 0, 'o'}, \
         {"gamma",            required_argument, 0, 'g'}, \
         {"ctf-cutoff",       required_argument, 0, 'c'}, \
         {"params-file",      required_argument, 0, 'p'}, \
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
int truncate_ctf_flag    = 0;
int normalize_flag       = 1;
int autocenter_vol_flag  = 1;
int invert_contrast_flag = 1;
int use_lambda_flag      = 0;
int fft_padding          = FFT_PADDING;

char const *mollifiers[] = {"", "delta", "gaussian", ""};

/*-------------------------------------------------------------------------------------------------*/

AiOpts *
new_AiOpts (void)
{
  CEXCEPTION_T e = EXC_NONE;
  AiOpts *opts = NULL;
  
  Try { opts = (AiOpts *) ali16_malloc (sizeof (AiOpts)); } CATCH_RETURN (e, NULL);
  
  opts->fname_in                = NULL;
  opts->fname_out               = NULL;
  opts->fname_params       = NULL;
  opts->fname_tiltangles        = NULL;
  opts->gamma                   = 0.0;
  opts->ctf_trunc               = 0.0;
  opts->moll_type               = GAUSSIAN;
  opts->num_images              = 0;
  opts->start_index             = 0;
  opts->lambda_pow              = 0.0;
  opts->bg_patch_ix0[2]         = -1;
  opts->bg_patch_shape[0]       = BG_PATCH_SIZE;
  opts->bg_patch_shape[1]       = BG_PATCH_SIZE;
  opts->bg_patch_shape[2]       = 1;
  
  return opts;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_free (AiOpts **popts)
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
AiOpts_print (AiOpts *opts)
{
  CAPTURE_NULL_VOID (opts);
  
  if (verbosity_level == VERB_LEVEL_QUIET)
    return;
  
  /* TODO: make dependent on verbosity */
  printf ("\n\n");
  puts ("Options from command line:");
  puts ("==========================\n");

  printf ("Input file            : %s\n", opts->fname_in);
  printf ("Output file           : %s\n", opts->fname_out);
  printf ("Reco parameter file   : %s\n", opts->fname_params);
  printf ("Tiltangles file       : %s\n", opts->fname_tiltangles);
  printf ("gamma                 : %e\n", opts->gamma);

  printf ("1/CTF cutoff          : ");
  if (truncate_ctf_flag)
    printf ("%e\n", opts->ctf_trunc);
  else
    printf ("(ignored)\n");
    
  printf ("Mollifier             : %s\n", mollifiers[opts->moll_type]); 

  printf ("Autocenter volume     : " );
  if (autocenter_vol_flag)
    printf ("yes\n");
  else
    printf ("no\n");

  printf ("Contrast inversion    : " );
  if (invert_contrast_flag)
    printf ("yes\n");
  else
    printf ("no\n");

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

  printf ("Number of images      : ");
  if (opts->num_images == -1)
    printf ("(determined from data)\n");
  else
    printf ("%d\n", opts->num_images);

  printf ("Start index           : %d\n", opts->start_index);

  printf ("Lambda                : ");
  if (use_lambda_flag)
    printf ("Lambda^a with a=%f\n", opts->lambda_pow);
  else  
    printf ("(not used)\n");
  
  printf ("FFT zero-padding      : %d\n", fft_padding);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void 
print_help (char const *progname)
{
  mollifier_type iter_m;
  // TODO: write conical tilt help
  
  printf ("Usage: %s [options] tiltseries_file\n", progname);
  puts ("");
  puts ("Required parameters:");
  puts ("");
  puts ("  -t file, --tiltangles-file=file");
  puts ("                 read tilt angles from FILE.");
  puts ("  -g value, --gamma=value");
  puts ("                 set regularization parameter gamma to VALUE; magnitude");
  puts ("                 should be in the order of detector pixel size [microns].");
  puts ("  -p file, --params-file=file");
  puts ("                 read CTF and reconstruction parameters from FILE (INI-style).");
  puts ("");
  puts ("Options:");
  puts ("");
  puts ("  -o file, --output-file=file");
  puts ("                 write reconstruction to FILE; if no parameter is given, the");
  puts ("                 output file is determined from tiltseries_file by appending");
  puts ("                 `_rec' before its extension.");
  puts ("  -c value, --ctf-cutoff=value");
  puts ("                 set value for reciprocal CTF cutoff; in intervals where ");
  puts ("                 1/CTF exceeds VALUE, it is replaced by a differentiable");
  puts ("                 transition spline; must be > 1.0 and should not exceed ~10.0.");
  puts ("  -m name, --mollifier=name");
  puts ("                 use the specified mollifier; Possible values:");
  for (iter_m = MO_START + 1; iter_m < MO_END; iter_m++)
    printf ("  %s", mollifiers[iter_m]);
  printf ("\n");
  puts ("                 (default: gaussian)");
  puts ("  -A, --autocenter-volume");
  puts ("                 automatically center the volume over the tilt-axis (enabled by");
  puts ("                 default).");
  puts ("  -N, --normalize");
  puts ("                 normalize the projection images based on their histograms (enabled");
  puts ("                 by default).");
  puts ("  --background=index_x,index_y,size_x,size_y");
  puts ("                 compute background statistics using a patch of SIZE_X x SIZE_Y");
  puts ("                 pixels with lower-left indices INDEX_X, INDEX_Y.");
  printf ("                 Default: 0,0,%d,%d\n", 
    BG_PATCH_SIZE, BG_PATCH_SIZE);
  puts ("  -I, --invert-contrast");
  puts ("                 invert the contrast of the images; use this option if dense regions");
  puts ("                 are darker than the background (enabled by default).");
  puts ("  -L value, --lambda-pow=value");
  puts ("                 feature reconstruction: compute Lambda^a(f) instead of f,");
  puts ("                 where Lambda = sqrt(-Laplacian) and a = VALUE.");
  puts ("  -n N, --num-images=N");
  puts ("                 override image number determined from input file by N; useful");
  puts ("                 for reconstruction with only a subset of the data.");
  puts ("  -s N, --start-index=N");
  puts ("                 start with N'th image instead of 0; useful for");
  puts ("                 reconstruction with only a subset of the data.");
  puts ("  --fft-padding[=N]");
  puts ("                 continue grid functions by N zero pixels in each direction prior");
  puts ("                 to computing Fourier transforms ('zero-padding').");
  printf ("                 (Default: N=%d)\n", FFT_PADDING);
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
  
  Try {
    *pbase_str = (char *) ali16_malloc (base_len + 1);
    *pext_str = (char *) ali16_malloc (ext_len + 1);
  } Catch (e)
  {
    EXC_RETHROW_REPRINT (e);
    free (*pbase_str); *pbase_str = NULL;
    free (*pext_str);  *pext_str  = NULL;
    return;
  }

  strncpy (*pbase_str, fname, base_len);  (*pbase_str)[base_len] = '\0';
  strncpy (*pext_str, pext, ext_len);  (*pext_str)[ext_len] = '\0';
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_fname_in (AiOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_in = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_in, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

#define REC_STR   "_rec"  /* append this before the file extension for the reco */

void
AiOpts_assemble_fname_out (AiOpts *opts,  char const *base, char const *ext)
{
  CEXCEPTION_T e = EXC_NONE;
  int base_len, rec_len, ext_len;
  char *p_tmp;
  
  base_len = strlen (base);
  rec_len  = strlen (REC_STR);
  ext_len  = strlen (ext);
  
  Try { opts->fname_out = (char *) ali16_malloc (base_len + rec_len + ext_len + 1); }
  CATCH_RETURN_VOID (e);
    
  p_tmp = opts->fname_out;
  strncpy (p_tmp, base, base_len);
  
  p_tmp += base_len;
  strncpy (p_tmp, REC_STR, rec_len);
  
  p_tmp += rec_len;
  strncpy (p_tmp, ext, ext_len);
  
  p_tmp += ext_len;
  *p_tmp = '\0';
    
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_determine_fname_out (AiOpts *opts, char *fname_in)
{
  CEXCEPTION_T e = EXC_NONE;
  char *bs, *ex;
  
  Try { fname_rsplit_at_dot (fname_in, &bs, &ex); } CATCH_RETURN_VOID (e);
  Try { AiOpts_assemble_fname_out (opts, bs, ex); } CATCH_RETURN_VOID (e);
    
  free (bs);
  free (ex);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_fname_out (AiOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
 
  Try { opts->fname_out = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_out, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_fname_params (AiOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);

  Try { opts->fname_params = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_params, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_fname_tiltangles (AiOpts *opts, char const *fname)
{
  CEXCEPTION_T e = EXC_NONE;
  int len = strlen (fname);
  
  Try { opts->fname_tiltangles = (char *) ali16_malloc (len + 1); } CATCH_RETURN_VOID (e);
  strncpy (opts->fname_tiltangles, fname, len);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_gamma (AiOpts *opts, float gamma)
{
  if (gamma <= 0.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Parameter of `-g' (`--gamma') must be positive.");
      return;
    }
    
  opts->gamma = gamma;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_ctf_trunc (AiOpts *opts, float ctf_cut)
{
  if (ctf_cut <= 0.0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Parameter of `-c' (`--ctf-cutoff') must be positive.");
      return;
    }
    
  opts->ctf_trunc = ctf_cut;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_num_images (AiOpts *opts, int n_images)
{
  if (n_images <= 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Parameter of `-n' (`--num-images') must be positive.");
      return;
    }
  
  opts->num_images = n_images;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_start_index (AiOpts *opts, int nstart)
{
  if (nstart < 0)
    {
      EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, 
        "Parameter of `-s' (`--start-index') must be nonnegative.");
      return;
    }
  
  opts->start_index = nstart;
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_moll_type (AiOpts *opts, char *moll_str)
{
  char *p;
  mollifier_type iter;
  
  for (p = moll_str; *p; p++) 
    *p = tolower (*p);
  
  for (iter = MO_START + 1; iter < MO_END; iter++)
    {
      if (strcmp (moll_str, mollifiers[iter]) == 0)
        {
          opts->moll_type = iter;
          return;
        }
    }
    
  EXC_THROW_CUSTOMIZED_PRINT (EXC_BADARG, "Unknown mollifier `%s'\n", moll_str);
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_set_bg_params (AiOpts *opts, char const *params_str)
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
AiOpts_set_lambda_pow (AiOpts *opts, float l_a)
{
  opts->lambda_pow = l_a;
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/

void
AiOpts_assign_from_args (AiOpts *opts, int argc, char **argv)
{
  CAPTURE_NULL_VOID (opts);

  /* Aux variables */
  int c;
  char const *progname = base_name(argv[0]);
  float gamma;
  float ctf_cut;
  int n_images;
  int nstart;
  float l_a;
  
  /* Internal flags for the short options */
  int F_flag = 0;
  int g_flag = 0;
  int m_flag = 0;
  int n_flag = 0;
  int o_flag = 0;
  int p_flag = 0;
  int P_flag = 0;
  int s_flag = 0;
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
          // This case should never occur!
          printf ("This should never occur! Option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'A':
          if (autocenter_vol_flag != 0)
            {
              fputs ("Invalid multiple use of `-A' (`--autocenter-volume') option.", stderr);
              exit (EXIT_FAILURE);
            }
            
          autocenter_vol_flag = 1;
          break;
          
        case 'c':
          if (truncate_ctf_flag != 0)
            {
              fputs ("Invalid multiple use of `-c' (`--ctf-cutoff') option.", stderr);
              exit (EXIT_FAILURE);
            }

          ctf_cut = atof (optarg);
          AiOpts_set_ctf_trunc (opts, ctf_cut);
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
          if (g_flag != 0)
            {
              fputs ("Invalid multiple use of `-g' (`--gamma') option.", stderr);
              exit (EXIT_FAILURE);
            }

          gamma = atof (optarg);
          AiOpts_set_gamma (opts, gamma);
          g_flag = 1;
          break;
 
        case 'I':
          invert_contrast_flag = 1;  /* Long option already sets the flag, but short one doesn't */
          break;

        case 'L':
          if (use_lambda_flag != 0)
            {
              fputs ("Invalid multiple use of `-L' (`--lambda-pow') option.", stderr);
              exit (EXIT_FAILURE);
            }

          l_a = atof (optarg);
          AiOpts_set_lambda_pow (opts, l_a);
          use_lambda_flag = 1;
          break;
 
        case 'm':
          if (m_flag != 0)
            {
              fputs ("Invalid multiple use of `-m' (`--mollifier') option.", stderr);
              exit (EXIT_FAILURE);
            }

          AiOpts_set_moll_type (opts, optarg);
          m_flag = 1;
          break;
 
        case 'n':
          if (n_flag != 0)
            {
              fputs ("Invalid multiple use of `-n' (`--num-images') option.", stderr);
              exit (EXIT_FAILURE);
            }

          n_images = atoi (optarg);
          AiOpts_set_num_images (opts, n_images);
          n_flag = 1;
          break;

        case 'N':
          normalize_flag = 1;  /* Long option already sets the flag, but short one doesn't */
          break;
 
        case 'o':
          if (o_flag != 0)
            {
              fputs ("Invalid multiple use of `-o' (`--output-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          AiOpts_set_fname_out (opts, optarg);
          o_flag = 1;
          break;
 
        case 'p':
          if (p_flag != 0)
            {
              fputs ("Invalid multiple use of `-p' (`--params-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          AiOpts_set_fname_params (opts, optarg);
          p_flag = 1;
          break;
 
        case 'P':
          if (P_flag != 0)
            {
              fputs ("Invalid multiple use of `--background' option.", stderr);
              exit (EXIT_FAILURE);
            }

          AiOpts_set_bg_params (opts, optarg);
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
          AiOpts_set_start_index (opts, nstart);
          s_flag = 1;
          break;
 
        case 't':
          if (t_flag != 0)
            {
              fputs ("Invalid multiple use of `-t' (`--tiltangles-file') option.", stderr);
              exit (EXIT_FAILURE);
            }

          AiOpts_set_fname_tiltangles (opts, optarg);
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
  
  if (!t_flag)
    {
      fprintf (stderr, "Error: `--tiltangles-file` option missing.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }

  if (!g_flag)
    {
      fprintf (stderr, "Error: `--gamma` option missing.\n\n" 
        "`%s --help' provides further information.\n\n", progname);
      exit (EXIT_FAILURE);
    }

  if (!p_flag)
    {
      fprintf (stderr, "Error: `--params-file` option missing.\n\n" 
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

  AiOpts_set_fname_in (opts, argv[optind]);

  
  if (!o_flag)
    AiOpts_determine_fname_out (opts, argv[optind]);
  
  return;
}

/*-------------------------------------------------------------------------------------------------*/
