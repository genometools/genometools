/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <stdbool.h>
#include "core/io.h"
#include "core/str.h"

typedef struct GtTokenizer GtTokenizer;

GtTokenizer*  gt_tokenizer_new(GtIO*); /* takes ownership */
/* activates the skipping of comment lines in the tokenizer (lines starting
   with '#') */
void          gt_tokenizer_skip_comment_lines(GtTokenizer*);
/* returns the current token */
GtStr*        gt_tokenizer_get_token(GtTokenizer*);
bool          gt_tokenizer_has_token(GtTokenizer*);
bool          gt_tokenizer_line_start(const GtTokenizer*);
void          gt_tokenizer_next_token(GtTokenizer*); /* go to the next token */
GtUword       gt_tokenizer_get_line_number(const GtTokenizer*);
const char*   gt_tokenizer_get_filename(const GtTokenizer*);
int           gt_tokenizer_unit_test(GtError*);
void          gt_tokenizer_delete(GtTokenizer*);

#endif
