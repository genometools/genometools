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

#include "annotationsketch/color.h"
#include "core/str.h"
#include "core/strarray.h"
#include "extended/genome_feature_type.h"
#include "extended/genome_node.h"

/* Holds style information. */
typedef struct GT_Style GT_Style;

/* Creates a new GT_Style object with given verbosity. If set, warnings will be
   given. */
GT_Style*      gt_style_new(bool verbose, Error*);
/* Creates a deep copy of the given GT_Style object. */
GT_Style*      gt_style_clone(const GT_Style*, Error*);
/* Loads and executes Lua style file with given <filename>.
   This file must contain a global table called 'style'. */
int            gt_style_load_file(GT_Style*, const char *filename, Error*);
/* Loads and executes Lua style code from the given String <instr>.
   This code must contain a global table called 'style'. */
int            gt_style_load_str(GT_Style*, Str *instr, Error*);
/* Generates Lua code which represents the given GT_Style object and
   writes it into the String object <outstr>.*/
int            gt_style_to_str(const GT_Style*, Str *outstr, Error*);
/* Reloads the Lua style file. */
void           gt_style_reload(GT_Style*);
/* Retrieves a color value from the GT_Style for <key>.
   If not set, false is returned and a default color is written. */
bool           gt_style_get_color(const GT_Style*, const char *section,
                                const char *key, Color*, GenomeNode*);
/* Sets a color value in the GT_Style for <key> (i.e., feature) to a
   certain value. */
void           gt_style_set_color(GT_Style*, const char *section,
                                const char *key, Color*);
/* Retrieve string value of <key> in <section>.
   If not set, false is returned. */
bool           gt_style_get_str(const GT_Style*, const char *section,
                              const char *key, Str*, GenomeNode*);
/* Set string <key> in <section> to <value>. */
void           gt_style_set_str(GT_Style*, const char *section, const char *key,
                               Str *value);
/* Retrieve numeric value of <key> in <section>.
   If not set, false is returned.*/
bool           gt_style_get_num(const GT_Style*, const char *section,
                              const char *key, double*, GenomeNode*);
/* Set numeric value of <key> in <section> to <number>. */
void           gt_style_set_num(GT_Style*, const char *section, const char *key,
                              double number);
/* Retrieve boolean value of <key> in <section>.
   If not set, false is returned.*/
bool           gt_style_get_bool(const GT_Style*, const char *section,
                               const char *key, bool*, GenomeNode*);
/* Set boolean value of <key> in <section>. */
void           gt_style_set_bool(GT_Style*, const char *section,
                                 const char *key, bool);
/* Unset value of <key> in <section>. */
void           gt_style_unset(GT_Style*, const char *section, const char *key);
/* Returns verbosity status. */
bool           gt_style_get_verbose(const GT_Style*);

int            gt_style_unit_test(Error*);
/* Deletes a GT_Style object but leaves the internal Lua state intact. */
void           gt_style_delete_without_state(GT_Style*);
void           gt_style_delete(GT_Style*);

#endif
