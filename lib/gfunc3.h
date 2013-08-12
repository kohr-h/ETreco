/*
 * gfunc3.h -- 3-dimensional grid functions
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

#ifndef __GFUNC_H__
#define __GFUNC_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vfunc.h"
#include "matvec3.h"

// TODO: update descriptions

typedef enum {REAL, HALFCOMPLEX} gfunc_type;

/** Structure for a function defined on a 3-dimensional cartesian grid.
 * 
 * Detailed description follows here.
 */
typedef struct
{
  /* Grid */

  /* essential data */
  vec3 x0;         /**< Origin of the grid */
  vec3 csize;      /**< Sizes of a grid cell */
  idx3 shape;      /**< Per-dimension number of grid cells */

  // derived data
  vec3 xmin;       /**< Minimal coordinates */
  vec3 xmax;       /**< Maximal coordinates */
  size_t ntotal;   /**< Total number of grid cells */

  // internal buffers
  vec3 _fbuf;      /**< Temporary storage for a vector */
  int _ntmp;       /**< Temporary storage for an integer */


  // status flags
  int is_initialized;  /**< 1 if grid is initialized and memory for fvals is allocated, 0 otherwise */
  int is_halfcomplex;  /**< If this flag is 1, the x dimension is halved and the values are 
                            interpreted as complex numbers. 0 means standard real-valued. */

  // function values
  float *fvals;    /**< float[ntotal] array holding the function values */
  
} gfunc3;


/*-------------------------------------------------------------------------------------------------*
 * Checks 
 *-------------------------------------------------------------------------------------------------*/

#define GFUNC_IS_2D(_pgf) ((_pgf)->shape[2] == 1)

#define GFUNC_CHECK_INIT_STATUS(_pgf)  \
do { \
  if ((_pgf)->is_initialized == 0) \
    EXC_THROW_PRINT (EXC_GFINIT); \
} while (0)

/*-------------------------------------------------------------------------------------------------*/
// Allocation
/*-------------------------------------------------------------------------------------------------*/

gfunc3 *
new_gfunc3 (void);

// Free memory of gfunc3 structure
void
gfunc3_free (gfunc3 **pgf);

/*-------------------------------------------------------------------------------------------------*/
// Structure initialization
/*-------------------------------------------------------------------------------------------------*/

// 
void
gfunc3_init (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type);

void
gfunc3_init_from_foreign_grid (gfunc3 *gf, gfunc3 const *gf_template);

void
gfunc3_set_csize (gfunc3 *gf, vec3 const cs);

// Initialize xmin and xmax from x0, csize and shape
void
gfunc3_compute_xmin_xmax (gfunc3 *gf);

/*-------------------------------------------------------------------------------------------------*/
// Function value initialization
/*-------------------------------------------------------------------------------------------------*/

void
gfunc3_assign_fvals_from_vfunc (gfunc3 *gf, const vfunc *vf);

/*-------------------------------------------------------------------------------------------------*/
// Screen output
/*-------------------------------------------------------------------------------------------------*/

// Print a grid summary on stdout, preceded by INTRO_TEXT.
void
gfunc3_print_grid (gfunc3 const *gf, char const *intro_text);

/*-------------------------------------------------------------------------------------------------*/
// Attributes
/*-------------------------------------------------------------------------------------------------*/

float
gfunc3_min (gfunc3 const *gf);

float
gfunc3_max (gfunc3 const *gf);

/* If called with PMEAN = NULL, the mean value is computed */
float
gfunc3_variance (gfunc3 const *gf, float const *pmean);

float
gfunc3_mean (gfunc3 const *gf);

// Check if the grids of GF1 and GF2 are equal up to a relative accuracy of EPS_GRID.  If yes,
// return 1, else 0.
int
gfunc3_grids_are_equal (gfunc3 const *gf1, gfunc3 const *gf2);

// Check if the grid of GF_SUB is a subgrid of GF's grid with relative accuracy EPS_GRID. If yes,
// return 1, else 0.
int
gfunc3_grid_is_subgrid (gfunc3 const *gf, gfunc3 const *gf_sub);

/*-------------------------------------------------------------------------------------------------*/
// Operations
/*-------------------------------------------------------------------------------------------------*/

// Copy SRC to DEST. Memory for DEST->DATA is (re)allocated
void
gfunc3_copy (gfunc3 *dest, gfunc3 const *src);

// Set GF1 <- A * GF1 + GF2. GF2 may be defined on a subgrid of GF1->GRID
void
gfunc3_axpy (float a, gfunc3 *gf1, gfunc3 const *gf2);

// Set GF <- A * GF + VF
void
gfunc3_axpy_vfunc (float a, gfunc3 *gf, const vfunc *vf);

// Multiply GF1 with GF2. GF2 may be defined on a subgrid of GF1->GRID
void
gfunc3_mul (gfunc3 *gf1, gfunc3 const *gf2);

// Multiply GF with VF
void
gfunc3_mul_vfunc (gfunc3 *gf, const vfunc *vf);

// Divide GF1 by GF2. GF2 may be defined on a subgrid of GF1->GRID
void
gfunc3_div (gfunc3 *gf1, gfunc3 const *gf2);

// Divide GF by VF
void
gfunc3_div_vfunc (gfunc3 *gf, const vfunc *vf);

// Multiply GF by the constant A
void
gfunc3_scale (gfunc3 *gf, float a);

// Add the constant C to GF
void
gfunc3_add_constant (gfunc3 *gf, float c);

// Translate (=shift the grid of) GF by the vector S
void
gfunc3_translate (gfunc3 *gf, vec3 const s);

// Only scale the grid of GF by the factor A
void
gfunc3_scale_grid (gfunc3 *gf, float a);

// Dilate GF by the constant A, i.e. scale the grid by A and multiply the
// values by |A|^(-d/2), where d is the actual dimension of GF
void
gfunc3_dilate (gfunc3 *gf, float a);

//  Set the value of GF at the given grid index IDX to the value PVAL points to (real or complex)
void
gfunc3_set (gfunc3 *gf, idx3 const idx, float const *pval);

// Set all values of GF to VAL
void
gfunc3_set_all (gfunc3 *gf, float const *pval);

/* Make GF non-negative by setting negative values to 0.0 */
void
gfunc3_make_nonneg (gfunc3 *gf);

// Set gf to zero at the boundary grid points
// void gfunc3_set_boundary_zero (gfunc3 gf);

/*-------------------------------------------------------------------------------------------------*/
// Evaluation
/*-------------------------------------------------------------------------------------------------*/

// Evaluate GF at the given grid index IDX
float
gfunc3_eval (gfunc3 *gf, idx3 const idx);

// Approximate the value of GF at PT by nearest neighbor interpolation
float
gfunc3_interp_nearest (gfunc3 const *gf, vec3 const pt);

// Approximate the value of GF at PT by nearest neighbor interpolation
float
gfunc3_interp_nearest_2d (gfunc3 const *gf, vec3 const pt);

// Approximate the value of GF at V by linear interpolation
float
gfunc3_interp_linear (gfunc3 const *gf, vec3 const pt);

// Approximate the value of GF at V by linear interpolation
float
gfunc3_interp_linear_2d (gfunc3 const *gf, vec3 const pt);

/*-------------------------------------------------------------------------------------------------*/
// Domain change
/*-------------------------------------------------------------------------------------------------*/

size_t *
gfunc3_subgrid_flatidcs (gfunc3 const *gf, gfunc3 const *gf_sub);

/* Pad gf with zeros; the padding vector gives the number of points that are added on both sides 
 * (per dimension) of the original grid
 */
void
gfunc3_zeropad (gfunc3 *gf, idx3 const padding);

/* Revert the effect of gfunc3_zeropad */
void
gfunc3_unpad (gfunc3 *gf, idx3 const padding);

/*
// Restrict gf to the subgrid sub_g. Superfluent values are discarded
void gfunc3_restrict (gfunc3 * gf, grid3 sub_g);

// Continue gf by zero to the bigger grid sup_g
void gfunc3_continue (gfunc3 * gf, grid3 sup_g);

// Set the values of gf equal to the values of the subfunction gf_sub where
// defined
void gfunc3_equate (gfunc3 gf, gfunc3 gf_sub);

void gfunc3_zeropad (gfunc3 * gf, int const *padding);


void gfunc3_unpad (gfunc3 * gf, int const *padding);
*/

#endif /* __GFUNC_H__ */
