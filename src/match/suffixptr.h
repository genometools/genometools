/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#ifndef SUFFIXPTR_H
#define SUFFIXPTR_H

#define SUFFIXPTRNEWVERSION
#ifdef  SUFFIXPTRNEWVERSION
typedef struct
{
  unsigned long value;
} Suffixptr;

#define SUFFIXPTRGET(TAB,IDX)     TAB[IDX].value
#define SUFFIXPTRSET(TAB,IDX,VAL) TAB[IDX].value = VAL
#define SUFFIXPTRDEREF(X)         (X)->value
#define SUFFIXPTRDEREFSET(X,VAL)  (X)->value = VAL

#else
typedef unsigned long Suffixptr;

#define SUFFIXPTRGET(TAB,IDX)     TAB[IDX]
#define SUFFIXPTRSET(TAB,IDX,VAL) TAB[IDX] = VAL
#define SUFFIXPTRDEREF(X)         *(X)
#define SUFFIXPTRDEREFSET(X,VAL)  *(X) = VAL
#endif

#endif
