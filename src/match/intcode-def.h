/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef INTCODE_DEF_H
#define INTCODE_DEF_H

#define GT_PREFIXLENBITS   4
#define GT_CODEBITS        (32 - GT_PREFIXLENBITS)
#define GT_MAXPREFIXLENGTH ((GT_CODEBITS) >> 1)
#define GT_MAXCODEVALUE    ((1U << (GT_CODEBITS)) - 1)

typedef struct
{
  unsigned int maxprefixindex:GT_PREFIXLENBITS;
  unsigned int code:GT_CODEBITS;
  GtUword position; /* get rid of this by using information from encseq */
} Codeatposition;

#endif
