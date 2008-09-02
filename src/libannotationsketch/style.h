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

#ifndef STYLE_H
#define STYLE_H

#include "lua.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libannotationsketch/color.h"
#include "libgtext/genome_feature_type.h"

/* Holds configuration information. */
typedef struct Style Style;

/* Creates a new Style object with given verbosity. If set, warnings will be
   given. */
Style*        style_new(bool verbose, Error*);
/* Creates a Style object wich reuses the given Lua state. */
Style*        style_new_with_state(lua_State*);
/* Creates a deep copy of the given Style object. */
Style*        style_clone(const Style*, Error*);
/* Loads and executes Lua style file with given <filename>.
   This file must contain a global table called 'Style'. */
int            style_load_file(Style*, Str *filename, Error*);
/* Loads and executes Lua style code from the given String <instr>.
   This code must contain a global table called 'Style'. */
int            style_load_str(Style*, Str *instr, Error*);
/* Generates Lua code which represents the given Style object and
   writes it into the String object <outstr>.*/
int            style_to_str(const Style*, Str *outstr, Error*);
/* Reloads the Lua style file. */
void           style_reload(Style*);
/* Retrieves a color value from the Style for <key>.
   If not set, false is returned and a default color is written. */
bool           style_get_color(const Style*, const char *section,
                                const char *key, Color*);
/* Sets a color value in the Style for <key> (i.e., feature) to a
   certain value. */
void           style_set_color(Style*, const char *section,
                                const char *key, Color*);
/* Retrieve string value of <key> in <section>.
   If not set, false is returned. */
bool           style_get_str(const Style*, const char *section,
                              const char *key, Str*);
/* Set string <key> in <section> to <value>. */
void           style_set_str(Style*, const char *section, const char *key,
                               Str *value);
/* Retrieve numeric value of <key> in <section>.
   If not set, false is returned.*/
bool           style_get_num(const Style*, const char *section,
                              const char *key, double*);
/* Set numeric value of <key> in <section> to <number>. */
void           style_set_num(Style*, const char *section, const char *key,
                              double number);
/* Retrieve boolean value of <key> in <section>.
   If not set, false is returned.*/
bool           style_get_bool(const Style*, const char *section,
                               const char *key, bool*);
/* Set boolean value of <key> in <section> to <number>. */
void           style_set_bool(Style*, const char *section, const char *key,
                               bool);
/* Unset value of <key> in <section>. */
void           style_unset(Style*, const char *section, const char *key);
/* Returns verbosity status. */
bool           style_get_verbose(const Style*);

int            style_unit_test(Error*);
/* Deletes a Style object but leaves the internal Lua state intact. */
void           style_delete_without_state(Style*);
void           style_delete(Style*);

#endif
