/*
  Copyright (c) 2007-2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#ifndef STYLE_API_H
#define STYLE_API_H

#include <stdbool.h>
#include "annotationsketch/color_api.h"
#include "core/error_api.h"
#include "core/str_api.h"

/* Holds AnnotationSketch style information. */
typedef struct GtStyle GtStyle;

/* Creates a new GtStyle object with given verbosity. If set, warnings will be
   given. */
GtStyle*      gt_style_new(bool verbose, GtError*);
/* Creates a deep copy of the given GtStyle object. */
GtStyle*      gt_style_clone(const GtStyle*, GtError*);
/* Loads and executes Lua style file with given <filename>.
   This file must contain a global table called 'style'. */
int            gt_style_load_file(GtStyle*, const char *filename, GtError*);
/* Loads and executes Lua style code from the given String <instr>.
   This code must contain a global table called 'style'. */
int            gt_style_load_str(GtStyle*, GtStr *instr, GtError*);
/* Generates Lua code which represents the given GtStyle object and
   writes it into the String object <outstr>.*/
int            gt_style_to_str(const GtStyle*, GtStr *outstr, GtError*);
/* Reloads the Lua style file. */
void           gt_style_reload(GtStyle*);
/* Sets a color value in the GtStyle for <key> (i.e., feature) to a
   certain value. */
void           gt_style_set_color(GtStyle*, const char *section,
                                  const char *key, const GtColor*);
/* Set string <key> in <section> to <value>. */
void           gt_style_set_str(GtStyle*, const char *section, const char *key,
                               GtStr *value);
/* Set numeric value of <key> in <section> to <number>. */
void           gt_style_set_num(GtStyle*, const char *section, const char *key,
                              double number);
/* Set boolean value of <key> in <section>. */
void           gt_style_set_bool(GtStyle*, const char *section,
                                 const char *key, bool);
/* Unset value of <key> in <section>. */
void           gt_style_unset(GtStyle*, const char *section, const char *key);
/* Returns verbosity status. */
bool           gt_style_get_verbose(const GtStyle*);

void           gt_style_delete(GtStyle*);

#endif
