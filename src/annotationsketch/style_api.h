/*
  Copyright (c) 2007-2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

/* Objects of the <GtStyle> class hold __AnnotationSketch__ style information
   like colors, margins, collapsing options, and others. The class provides
   methods to set values of various types. Each value is organised into
   a __section__ and is identified by a __key__. That is, a __section__, __key__
   pair must uniquely identify a value. */
typedef struct GtStyle GtStyle;

/* Creates a new <GtStyle> object. */
GtStyle*      gt_style_new(GtError*);
/* Increments the reference count of the given <GtStyle>. */
GtStyle*      gt_style_ref(GtStyle*);
/* Enables unsafe mode (``io'' and ``os'' libraries loaded). */
void           gt_style_unsafe_mode(GtStyle*);
/* Enables safe mode (``io'' and ``os'' libraries not accessible). */
void           gt_style_safe_mode(GtStyle*);
/* Returns true if <sty> is in unsafe mode. */
bool           gt_style_is_unsafe(GtStyle *sty);
/* Creates a independent (``deep'') copy of the given <GtStyle> object. */
GtStyle*      gt_style_clone(const GtStyle*, GtError*);
/* Loads and executes Lua style file with given <filename>.
   This file must define a global table called __style__. */
int            gt_style_load_file(GtStyle*, const char *filename, GtError*);
/* Loads and executes Lua style code from the given <GtStr> <instr>.
   This code must define a global table called __style__. */
int            gt_style_load_str(GtStyle*, GtStr *instr, GtError*);
/* Generates Lua code which represents the given <GtStyle> object and
   writes it into the <GtStr> object <outstr>.*/
int            gt_style_to_str(const GtStyle*, GtStr *outstr, GtError*);
/* Reloads the Lua style file. */
void           gt_style_reload(GtStyle*);
/* Sets a color value in the <GtStyle> for section <section> and <key> to a
   certain <color>. */
void           gt_style_set_color(GtStyle*, const char *section,
                                  const char *key, const GtColor *color);
/* Set string with key <key> in <section> to <value>. */
void           gt_style_set_str(GtStyle*, const char *section, const char *key,
                               GtStr *value);
/* Set numeric value of key <key> in <section> to <number>. */
void           gt_style_set_num(GtStyle*, const char *section, const char *key,
                              double number);
/* Set boolean value of key <key> in <section> to <val>. */
void           gt_style_set_bool(GtStyle*, const char *section,
                                 const char *key, bool val);
/* Unset value of key <key> in <section>. */
void           gt_style_unset(GtStyle*, const char *section, const char *key);
/* Deletes this <style>. */
void           gt_style_delete(GtStyle *style);

#endif
