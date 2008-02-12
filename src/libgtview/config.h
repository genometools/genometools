/*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef CONFIG_H
#define CONFIG_H

#include "lua.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libgtview/color.h"
#include "libgtext/genome_feature_type.h"

/* Represents domination status of an ordered pair.
   Used when two different types collapse into the same parent
   to determine splitting precedence. */
typedef enum
{
  DOMINATES_FIRST,
  DOMINATES_SECOND,
  DOMINATES_EQUAL,
  DOMINATES_NOT_SPECIFIED,
  DOMINATES_UNKNOWN_TYPE
} DominateStatus;

/* Holds configuration info for the libgtview classes. */
typedef struct Config Config;

/* Create a new Config object with given verbosity. If set, warnings will be
   given. */
Config*        config_new(bool verbose, Error*);
/* Create a Config object wich reuses the given Lua state. */
Config*        config_new_with_state(lua_State*);
/* Load and executes a Lua configuration file with given <filename>.
   This file must contain a global table called 'config'. */
int            config_load_file(Config*, Str *filename, Error*);
/* Reload the Lua configuration file. */
void           config_reload(Config*);
/* Retrieve a color value from the configuration for <key> (i.e., feature). */
Color          config_get_color(const Config*, const char *key);
/* Similar to previous function. Necessary for Ruby bindings, because
   apparently 'dl/import' cannot handle returned structs. */
void           config_get_colorptr(const Config*, Color*, const char *key);
/* Sets a color value in the configuration for <key> (i.e., feature) to a
   certain value. */
void           config_set_color(Config*, const char *key, Color*);
/* Retrieve string value of <key> in <section>.
   If not set, <deflt> is returned. */
const char*    config_get_cstr(const Config*, const char *section,
                               const char *key, const char *deflt);
/* Set string <key> in <section> to <value>. */
void           config_set_cstr(Config*, const char *section, const char *key,
                               const char *value);
/* Retrieve numeric value of <key> in <section>.
   If not set, <deflt> is returned.*/
double         config_get_num(const Config*, const char *section,
                              const char *key, double deflt);
/* Set numeric value of <key> in <section> to <number>. */
void           config_set_num(Config*, const char *section, const char *key,
                              double number);
/* Retrieve a list of string (as StrArray) for <key> in <section>, returns NULL
   on error. */
StrArray*      config_get_cstr_list(const Config*, const char *section,
                                    const char *key);
/* Set <key> in <section> to <list> of strings. */
void           config_set_cstr_list(Config*, const char *section,
                                    const char *key, StrArray *list);
/* Check if <checkstr> appears in list of strings named <key> in <section>. */
bool           config_cstr_in_list(const Config*, const char *section,
                                   const char *key, const char *checkstr);
/* Returns verbosity status. */
bool           config_get_verbose(const Config*);

/* Compares two GenomeFeatureTypes <gft1> and <gft2> w.r.t. their splitting
   precendence as defined in the config object.
   If a type dominates, it will be drawn on top of the other in the image. */
DominateStatus config_dominates(Config*, GenomeFeatureType gft1,
                                GenomeFeatureType gft2);
int            config_unit_test(Error*);
/* Deletes a Config object but leaves the internal Lua state intact. */
void           config_delete_without_state(Config*);
void           config_delete(Config*);

#endif
