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
#include "core/str.h"

/* Holds AnnotationSketch style information. */
typedef struct GT_Style GT_Style;

/* Creates a new GT_Style object with given verbosity. If set, warnings will be
   given. */
GT_Style*      gt_style_new(bool verbose, GT_Error*);
/* Creates a deep copy of the given GT_Style object. */
GT_Style*      gt_style_clone(const GT_Style*, GT_Error*);
/* Loads and executes Lua style file with given <filename>.
   This file must contain a global table called 'style'. */
int            gt_style_load_file(GT_Style*, const char *filename, GT_Error*);
/* Loads and executes Lua style code from the given String <instr>.
   This code must contain a global table called 'style'. */
int            gt_style_load_str(GT_Style*, GT_Str *instr, GT_Error*);
/* Generates Lua code which represents the given GT_Style object and
   writes it into the String object <outstr>.*/
int            gt_style_to_str(const GT_Style*, GT_Str *outstr, GT_Error*);
/* Reloads the Lua style file. */
void           gt_style_reload(GT_Style*);
/* Sets a color value in the GT_Style for <key> (i.e., feature) to a
   certain value. */
void           gt_style_set_color(GT_Style*, const char *section,
                                  const char *key, const GT_Color*);
/* Set string <key> in <section> to <value>. */
void           gt_style_set_str(GT_Style*, const char *section, const char *key,
                               GT_Str *value);
/* Set numeric value of <key> in <section> to <number>. */
void           gt_style_set_num(GT_Style*, const char *section, const char *key,
                              double number);
/* Set boolean value of <key> in <section>. */
void           gt_style_set_bool(GT_Style*, const char *section,
                                 const char *key, bool);
/* Unset value of <key> in <section>. */
void           gt_style_unset(GT_Style*, const char *section, const char *key);
/* Returns verbosity status. */
bool           gt_style_get_verbose(const GT_Style*);

void           gt_style_delete(GT_Style*);

#endif
