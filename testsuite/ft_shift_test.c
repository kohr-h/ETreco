/*
 * ft_shift_test.c - test if shifts are handled correctly by Fourier transforms
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

#include <config.h>

#include "ETreco.h"

int main(int argc, char **argv)
{
  CEXCEPTION_T e;
  
  idx3 shp = {20, 31}, pd = {20, 15};
  vec3 cs = {0.3, 3.2}, sft = {-1.4, 21.5};
  
  float ones[2] = {1.0, 1.0};
  
  Try
  {
    gfunc3 *img = new_gfunc3 ();
    
    gfunc3_init (img, NULL, cs, shp, REAL);
    gfunc3_set_all (img, ones);
    gfunc3_print_grid (img, "initial grid");
    temp_mrc_out (img, "img_orig_", 1);

    gfunc3_zeropad (img, pd);
    temp_mrc_out (img, "img_zp_", 1);

    gfunc3_translate (img, sft);

    fft_forward (img);
    temp_mrc_out (img, "ft_shift", 1);
    
    fft_backward (img);
    temp_mrc_out (img, "img_ift_ft_shift_", 1);
  }
  Catch (e)
  {
    EXC_REPRINT;
  }
  
  return 0;
}

