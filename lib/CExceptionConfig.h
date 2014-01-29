/*
 * CExceptionConfig.h - config file for CException
 * 
 * Copyright 2014 Holger Kohr <kohr@num.uni-sb.de>
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
 
#ifndef __CEXCEPTIONCONFIG_H__
#define __CEXCEPTIONCONFIG_H__

/*-------------------------------------------------------------------------------------------------
 * Exception codes and descriptions
 *-------------------------------------------------------------------------------------------------*/

typedef enum {EXC_NONE, 
              EXC_NOMEMORY,
              EXC_NULL, 
              EXC_BADARG, 
              EXC_INDEX, 
              EXC_GFINIT, 
              EXC_GFTYPE,
              EXC_GFDIM,
              EXC_SUBGRID, 
              EXC_IO,
              EXC_COMPUTE,
              EXC_UNIMPL, 
              EXC_UNSPEC,
              EXC_ID1,
              EXC_ID2} 
              exception_code;

static char const *_exc_names[] __attribute__ ((unused)) = {
  "None", 
  "OutOfMemory exception",
  "NullPointer exception", 
  "Argument exception",
  "IndexRange exception",
  "FunctionInitStatus exception",
  "FunctionType exception",
  "FunctionDimension exception",
  "Subgrid exception",
  "I/O exception",
  "Computational exception",
  "NotImplemented exception",
  "Unspecified exception"
  "",
  ""
  };

/*-------------------------------------------------------------------------------------------------
 * Exception throwing macros
 *-------------------------------------------------------------------------------------------------*/

#define EXC_THROW_PRINT(_exc_id) \
do { \
  fprintf (stderr, "At %s() [%s:%d]: %s\n", __func__, __FILE__, __LINE__, _exc_names[_exc_id]); \
  Throw (_exc_id); \
} while (0)

#define EXC_THROW_CUSTOMIZED_PRINT(_exc_id, _format, ...) \
do { \
  fprintf (stderr, "At %s() [%s:%d]: %s: ", __func__, __FILE__, __LINE__, _exc_names[_exc_id]); \
  fprintf (stderr, _format, ##__VA_ARGS__); \
  fprintf (stderr, "\n"); \
  Throw (_exc_id); \
} while (0)

#define EXC_RETHROW_REPRINT(_exc_id) \
do { \
  fprintf (stderr, "raised at %s() [%s:%d]:\n", __func__, __FILE__, __LINE__); \
  Throw (_exc_id); \
} while (0)

#define EXC_REPRINT \
do { \
  fprintf (stderr, "raised at %s() [%s:%d].\n", __func__, __FILE__, __LINE__); \
} while (0)


#define CATCH_RETURN_VOID(_exc) Catch (_exc) {EXC_RETHROW_REPRINT (_exc); return;}
#define CATCH_RETURN(_exc, _retval) Catch (_exc) {EXC_RETHROW_REPRINT (_exc); return _retval;}
#define CATCH_EXIT_FAIL(_exc) Catch (_exc) {EXC_REPRINT; exit (EXIT_FAILURE);}


/*-------------------------------------------------------------------------------------------------*/

/* Optionally define the exception type (something like an int which can be directly assigned) */
#define CEXCEPTION_T    exception_code

/* Optionally define the reserved value representing NO EXCEPTION */
#define CEXCEPTION_NONE (EXC_NONE)

/* Optionally define a special handler for unhandled exceptions */
#define CEXCEPTION_NO_CATCH_HANDLER(id) fprintf (stderr, "Uncaught exception %d\n", id)

/*-------------------------------------------------------------------------------------------------*/

#endif // __CEXCEPTIONCONFIG_H__
