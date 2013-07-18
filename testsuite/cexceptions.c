/*
 * cexceptions.c
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

#include "ETreco.h"


void foo (void *p)
{
  static int count = 0;
  count++;
  if (p == NULL)
    EXC_THROW_CUSTOMIZED_PRINT (EXC_NULL, "Null Pointer %d", count);
  
  puts ("Everything OK.");
  
  return;
}


void bar (void *p)
{
  CEXCEPTION_T e = EXC_NONE;

  Try {foo (p);}
  Catch (e) {EXC_RETHROW_REPRINT (e);}
  
  puts ("bar: no exception, returning normally");
  
  return;
}

void baz (void *p)
{
  CEXCEPTION_T e = EXC_NONE;

  Try 
  {
    bar (p);
  }
  Catch (e)
  {
    puts ("baz: rethrowing..");
    EXC_THROW_PRINT (e);
  }
  
  return;
}

int main(int argc, char **argv)
{
  CEXCEPTION_T e = EXC_NONE;
  int i, *pi = &i;
  void *p = NULL; 
  gfunc3 *gf;
  
  Try
  {
    gf = new_gfunc3 ();
    // bar (pi);
    baz (p);
  }
  Catch (e)
  {
    EXC_REPRINT;
    return e;
  }
  
  return 0;
}

