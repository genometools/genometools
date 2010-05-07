/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef SPLITTER_API_H
#define SPLITTER_API_H

#include "core/error_api.h"

/* The <GtSplitter> class defines objects which can split given strings into
   tokens delimited by a given character, allowing for convenient access to
   each token. */
typedef struct GtSplitter GtSplitter;

/* Create a new <GtSplitter> object. */
GtSplitter*   gt_splitter_new(void);

/* Split <string> of given <length> into tokens delimited by 'delimiter'.
   Note that <string> is modified in the splitting process! */
void          gt_splitter_split(GtSplitter*, char *string, unsigned long length,
                                char delimiter);

/* Get all tokens in an array. */
char**        gt_splitter_get_tokens(GtSplitter*);

/* Get token with number <token_num>. */
char*         gt_splitter_get_token(GtSplitter*, unsigned long token_num);

/* Reset the splitter. */
void          gt_splitter_reset(GtSplitter*);

/* Returns the number of tokens. */
unsigned long gt_splitter_size(GtSplitter*);

/* Delete the <GtSplitter> object. */
void          gt_splitter_delete(GtSplitter*);

#endif
