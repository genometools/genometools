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

/* Holds configuration info for the libgtview classes. */
typedef struct Config Config;

/* Creates a new Config object with given verbosity. If set, warnings will be
   given. */
Config*        config_new(bool verbose, Error*);
/* Creates a Config object wich reuses the given Lua state. */
Config*        config_new_with_state(lua_State*);
/* Creates a deep copy of the given Config object. */
Config*        config_clone(Config*, Error*);
/* Loads and executes a Lua configuration file with given <filename>.
   This file must contain a global table called 'config'. */
int            config_load_file(Config*, Str *filename, Error*);
/* Loads and executes a Lua configuration code from the given String <instr>.
   This code must contain a global table called 'config'. */
int            config_load_str(Config*, Str *instr);
/* Generates Lua code which represents the given Config object and
   writes it into the String object <outstr>.*/
int            config_to_str(Config*, Str *outstr);
/* Reloads the Lua configuration file. */
void           config_reload(Config*);
/* Retrieves a color value from the configuration for <key>.
   If not set, false is returned and a default color is written. */
bool           config_get_color(const Config*, const char *section,
                                const char *key, Color*);
/* Sets a color value in the configuration for <key> (i.e., feature) to a
   certain value. */
void           config_set_color(Config*, const char *section,
                                const char *key, Color*);
/* Retrieve string value of <key> in <section>.
   If not set, false is returned. */
bool           config_get_str(const Config*, const char *section,
                              const char *key, Str*);
/* Set string <key> in <section> to <value>. */
void           config_set_str(Config*, const char *section, const char *key,
                               Str *value);
/* Retrieve numeric value of <key> in <section>.
   If not set, false is returned.*/
bool           config_get_num(const Config*, const char *section,
                              const char *key, double*);
/* Set numeric value of <key> in <section> to <number>. */
void           config_set_num(Config*, const char *section, const char *key,
                              double number);
/* Retrieve boolean value of <key> in <section>.
   If not set, false is returned.*/
bool           config_get_bool(const Config*, const char *section,
                               const char *key, bool*);
/* Set boolean value of <key> in <section> to <number>. */
void           config_set_bool(Config*, const char *section, const char *key,
                               bool);
/* Unset value of <key> in <section>. */
void           config_unset(Config*, const char *section, const char *key);
/* Returns verbosity status. */
bool           config_get_verbose(const Config*);

int            config_unit_test(Error*);
/* Deletes a Config object but leaves the internal Lua state intact. */
void           config_delete_without_state(Config*);
void           config_delete(Config*);

#endif
