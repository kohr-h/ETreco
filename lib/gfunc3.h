/*
 * gfunc3.h -- 3-dimensional grid functions
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

#ifndef __GFUNC3_H__
#define __GFUNC3_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex.h>

#include "vec3.h"

#include "vfunc.h"


typedef enum {REAL, HALFCOMPLEX, COMPLEX} gfunc_type;
typedef enum {EUCLIDEAN, POLAR} grid_type;

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

  /* derived data */
  vec3 xmin;       /**< Minimal coordinates */
  vec3 xmax;       /**< Maximal coordinates */
  size_t ntotal;   /**< Total number of grid cells */

  /* internal buffers */
  vec3 _fbuf;      /**< Temporary storage for a vector */
  int _ntmp;       /**< Temporary storage for an integer */


  /* status flags */
  int is_initialized;  /**< 1 if grid is initialized and memory for fvals is allocated, 0 otherwise */
  gfunc_type type;     /**< If HALFCOMPLEX, the x dimension is halved and the values are 
                            interpreted as complex numbers. COMPLEX means that FVALS holds 
                            2 * NTOTAL entries. */

  /* function values */
  float *fvals;    /**< Array holding the function values. Size depends on TYPE. */
  
} gfunc3;


/*-------------------------------------------------------------------------------------------------*
 * Checks 
 *-------------------------------------------------------------------------------------------------*/

#define GFUNC_IS_2D(_pgf) ((_pgf)->shape[2] == 1)
#define GFUNC_IS_REAL(_pgf) ((_pgf)->type == REAL)
#define GFUNC_IS_COMPLEX(_pgf) (((_pgf)->type == HALFCOMPLEX) || ((_pgf)->type == COMPLEX))

#define GFUNC_CAPTURE_UNINIT(_pgf, _retval)  \
do { \
  if ((_pgf)->is_initialized == 0) { EXC_THROW_PRINT (EXC_GFINIT); return _retval; }\
} while (0)

#define GFUNC_CAPTURE_UNINIT_VOID(_pgf)  \
do { \
  if ((_pgf)->is_initialized == 0) { EXC_THROW_PRINT (EXC_GFINIT); return; }\
} while (0)



/*-------------------------------------------------------------------------------------------------*
 * Allocation
 *-------------------------------------------------------------------------------------------------*/

/* Create a new gfunc3 structure and return a pointer to it.
 * 
 * Thrown exceptions:  
 * - Rethrows
 */
gfunc3 *
new_gfunc3 (void);


/* Free memory of gfunc3 structure.
 * 
 * Thrown exceptions: 
 */
void
gfunc3_free (gfunc3 **pgf);



/*-------------------------------------------------------------------------------------------------*
 * Structure initialization
 *-------------------------------------------------------------------------------------------------*/

/* Initialize the grid of GF with the provided parameters X0, CS and SHP.  The function type 
 * GF_TYPE (REAL or HALFCOMPLEX) determines the layout of the data array, which is allocated 
 * according to the grid shape SHP.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - Rethrows
 */
void
gfunc3_init (gfunc3 *gf, vec3 const x0, vec3 const cs, idx3 const shp, gfunc_type gf_type);


/* Initialize grid and data array of GF from GF_TEMPLATE.  The data array of GF is free'd if 
 * necessary.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_init_from_foreign_grid (gfunc3 *gf, gfunc3 const *gf_template);


/* Set cell size of GF to CS and recompute the grid.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_set_csize (gfunc3 *gf, vec3 const cs);


/* Set the origin of GF to X0 and recompute the grid.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_set_x0 (gfunc3 *gf, vec3 const x0);


/* Compute minimal and maximal coordinates of GF using origin, shape and cell size.
 * 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_compute_xmin_xmax (gfunc3 *gf);



/*-------------------------------------------------------------------------------------------------*
 * Function value initialization
 *-------------------------------------------------------------------------------------------------*/

/* Assign the values in the data array of GF according to the abstract vector function VF.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 */
void
gfunc3_assign_fvals_from_vfunc (gfunc3 *gf, const vfunc *vf);

/*-------------------------------------------------------------------------------------------------*
 * Screen output
 *-------------------------------------------------------------------------------------------------*/

/* Print a grid summary of GF to standard output, preceded by INTRO_TEXT.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_print_grid (gfunc3 const *gf, char const *intro_text);



/*-------------------------------------------------------------------------------------------------*
 * Attributes
 *-------------------------------------------------------------------------------------------------*/

/* Return the minimum value in the data array of GF.  Works for REAL functions only.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 */ 
float
gfunc3_min (gfunc3 const *gf);


/* Return the maximum value in the data array of GF.  Works for REAL functions only.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 */ 
float
gfunc3_max (gfunc3 const *gf);


/* Return the mean value in the data array of GF.  Currently works for REAL functions only.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */ 
float
gfunc3_mean (gfunc3 const *gf);

/* Return the variance of the data array of GF.  If PMEAN is NULL, the mean value is computed.
 * Currently works for REAL functions only.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */ 
float
gfunc3_variance (gfunc3 const *gf, float const *pmean);



/*-------------------------------------------------------------------------------------------------*
 * Comparison
 *-------------------------------------------------------------------------------------------------*/


/* Return TRUE if the grids of GF1 and GF2 are equal up to a relative accuracy of EPS_GRID, 
 * FALSE otherwise.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
int
gfunc3_grids_are_equal (gfunc3 const *gf1, gfunc3 const *gf2);


/* Return TRUE if the grid of GF_SUB is a subgrid of GF's grid with relative accuracy EPS_GRID, 
 * FALSE otherwise.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
int
gfunc3_grid_is_subgrid (gfunc3 const *gf, gfunc3 const *gf_sub);



/*-------------------------------------------------------------------------------------------------*
 * Operations
 *-------------------------------------------------------------------------------------------------*/


/* Copy grid and data of SRC to DEST.  Memory for DEST->DATA is (re)allocated.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - Rethrows
 */
void
gfunc3_copy (gfunc3 *dest, gfunc3 const *src);


/* Store the real part of the COMPLEX type GF in RE and return RE. If RE is NULL, a new gfunc3 is 
 * initialized with the real part values and returned. 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 * - EXC_SUBGRID
 * - Rethrows
 */
gfunc3 *
gfunc3_realpart (gfunc3 const *gf, gfunc3 *re);


/* Store the imaginary part of the COMPLEX type GF in IM and return IM. If IM is NULL, a new gfunc3 
 * is initialized with the imaginary part values and returned.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 * - EXC_SUBGRID
 * - Rethrows
 */
gfunc3 *
gfunc3_imagpart (gfunc3 const *gf, gfunc3 *im);


/* Transfer the REAL type function GF to a COMPLEX function by allocating sufficient space and 
 * use the old values as the real part of the new values.
 */
void
gfunc3_real2complex (gfunc3 *gf);


/* Transfer the REAL type function GF to a COMPLEX function by allocating sufficient space and 
 * use the old values as the imaginary part of the new values.
 */
void
gfunc3_imag2complex (gfunc3 *gf);


/* Swap x and z axes of GF in-place. Not supporting HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 */ 
void 
gfunc3_swapxz (gfunc3 *gf);


/* Set GF1 <- A * GF1 + GF2.  GF2 may be defined on a subgrid of GF1's grid.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - Rethrows
 */ 
void
gfunc3_axpy (float a, gfunc3 *gf1, gfunc3 const *gf2);


/* Set GF <- A * GF + VF 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 */
void
gfunc3_axpy_vfunc (float a, gfunc3 *gf, const vfunc *vf);


/* Set GF1 <- GF1 * GF2.  GF2 may be defined on a subgrid of GF1's GRID.  Mixed multiplication of REAL
 * and HALFCOMPLEX is currently only possible for GF1 being HALFCOMPLEX.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - Rethrows
 */ 
void
gfunc3_mul (gfunc3 *gf1, gfunc3 const *gf2);


/* Set GF <- GF * VF. 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 */
void
gfunc3_mul_vfunc (gfunc3 *gf, const vfunc *vf);


/* Set GF1 <- GF1 / GF2.  GF2 may be defined on a subgrid of GF1's GRID.  Currently, only REAL 
 * multiplication is implemented.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - Rethrows
 */
void
gfunc3_div (gfunc3 *gf1, gfunc3 const *gf2);


/* Set GF <- GF / VF. 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */
void
gfunc3_div_vfunc (gfunc3 *gf, const vfunc *vf);


/* Set GF <- A * GF. If GF is of REAL type and A has nonzero imaginary part, GF is converted to 
 * a COMPLEX function.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 */
void
gfunc3_scale (gfunc3 *gf, float complex a);


/* Set GF <- GF + C.  Currently implemented only for REAL function.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */
void
gfunc3_add_constant (gfunc3 *gf, float c);


/* Translate (= shift the grid of) GF by the vector S.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 */
void
gfunc3_translate (gfunc3 *gf, vec3 const s);


/* Scale the grid of GF by the factor A.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_BADARG
 */
void
gfunc3_scale_grid (gfunc3 *gf, float a);


/* Dilate GF by the constant A, i.e. scale the grid by A and multiply the values by |A|^(-d/2), 
 * where d=2,3 is the dimension of GF.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_BADARG
 * - Rethrows
 */
void
gfunc3_dilate (gfunc3 *gf, float a);


/*  Set the value of GF at the given grid index IDX to the value PVAL points to (real or complex).
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_INDEX
 */
void
gfunc3_set (gfunc3 *gf, idx3 const idx, float complex val);


/* Set all values of GF to VAL. *
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 */
void
gfunc3_set_all (gfunc3 *gf, float complex val);


/* Make GF non-negative by setting negative values to 0.0f *
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_GFTYPE
 */
void
gfunc3_make_nonneg (gfunc3 *gf);



/*-------------------------------------------------------------------------------------------------*
 * Evaluation
 *-------------------------------------------------------------------------------------------------*/


/* Evaluate GF at the given grid index IDX.  Currently not implemented for HALFCOMPLEX functions.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_INDEX
 */
float
gfunc3_eval (gfunc3 *gf, idx3 const idx);


/* Approximate the value of GF at PT by nearest neighbor interpolation on the grid.  Currently not 
 * implemented for HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */
float
gfunc3_interp_nearest (gfunc3 const *gf, vec3 const pt);


/* Approximate the value of GF at PT by nearest neighbor interpolation on the grid. This is an 
 * optimized 2D version of the generic function.  Currently not implemented for HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_GFDIM
 */
float
gfunc3_interp_nearest_2d (gfunc3 const *gf, vec3 const pt);


/* Approximate the value of GF at PT by linear interpolation on the grid.  Currently not implemented 
 * for HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */
float
gfunc3_interp_linear (gfunc3 const *gf, vec3 const pt);


/* Approximate the value of GF at PT by linear interpolation on the grid. This is an optimized 2D 
 * version of the generic function.  Currently not implemented for HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_GFDIM
 */
float
gfunc3_interp_linear_2d (gfunc3 const *gf, vec3 const pt);


/* Store the coordinates of GF's grid in a float* array. This array is allocated in the function and 
 * returned. Interpret the grid according to GRID_TYPE as the transform of an Euclidean grid.
 * 
 * Thrown exceptions:
 * - EXC_NULL
 */

float *
gfunc3_grid_points (gfunc3 const *gf, grid_type gridtype);



/*-------------------------------------------------------------------------------------------------*
 * Domain change
 *-------------------------------------------------------------------------------------------------*/


/* TODO: make this function private? */
/* Return the flattened indices of the subgrid (GF_SUB's grid) inside the large grid (GF's grid).
 * 
 * Thrown exceptions:
 * - EXC_NULL
 * - EXC_SUBGRID
 * - Rethrows
 */
size_t *
gfunc3_subgrid_flatidcs (gfunc3 const *gf, gfunc3 const *gf_sub);


/* Pad GF with zeros; the PADDING vector gives the number of points that are added on both sides 
 * (per dimension) of the original grid of GF.  Memory for the new array is allocated.  Currently 
 * not implemented for HALFCOMPLEX functions.
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 */
void
gfunc3_zeropad (gfunc3 *gf, idx3 const padding);


/* Revert the effect of gfunc3_zeropad.  Currently not implemented for HALFCOMPLEX functions. 
 * 
 * Thrown exceptions: 
 * - EXC_NULL
 * - EXC_GFINIT
 * - EXC_UNIMPL
 * - EXC_BADARG
 */
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

#endif /* __GFUNC3_H__ */
