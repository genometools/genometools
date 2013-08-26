/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>

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

#ifndef WORD_H
#define WORD_H

/* Define GtUword as an unsigned integer with the machine word size (4 byte on
   32-bit systems and 8 byte on 64-bit systems). */
#if defined (_LP64) || defined (_WIN64)
typedef unsigned long long GtUword;
#else
typedef unsigned int GtUword;
#endif

/* Define GtWord as a signed integer with the machine word size (4 byte on
   32-bit systems and 8 byte on 64-bit systems). */
#if defined (_LP64) || defined (_WIN64)
typedef long long GtWord;
#else
typedef int GtWord;
#endif

#endif
