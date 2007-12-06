/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/io.h"
#include "libgtcore/str.h"

typedef struct Tokenizer Tokenizer;

Tokenizer*    tokenizer_new(IO*); /* takes ownership */
/* activates the skipping of comment lines in the tokenizer (lines starting
   with '#') */
void          tokenizer_skip_comment_lines(Tokenizer*);
/* returns the current token */
Str*          tokenizer_get_token(Tokenizer*);
bool          tokenizer_has_token(Tokenizer*);
bool          tokenizer_line_start(const Tokenizer*);
void          tokenizer_next_token(Tokenizer*); /* go to the next token */
unsigned long tokenizer_get_line_number(const Tokenizer*);
const char*   tokenizer_get_filename(const Tokenizer*);
int           tokenizer_unit_test(Error*);
void          tokenizer_delete(Tokenizer*);

#endif
