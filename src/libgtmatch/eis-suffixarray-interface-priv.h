/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_SUFFIXARRAY_INTERFACE_PRIV_H
#define EIS_SUFFIXARRAY_INTERFACE_PRIV_H

#include "eis-suffixarray-interface.h"

/**
 * Used to pass suffixarray to read from and alphabet to encode with
 * to readers.
 */
struct suffixarrayFileInterface
{
  struct SASeqSrc baseClass;
  struct saTaggedXltorStateList xltorStates;
  Suffixarray *sa;                      /**< the suffix array to read from */
  char numBWTFileReaders;
};

static inline SuffixarrayFileInterface *
SASS2SAI(SASeqSrc *baseClass)
{
  return (SuffixarrayFileInterface *)
    ((char *)baseClass - offsetof(SuffixarrayFileInterface, baseClass));
}

static inline const SuffixarrayFileInterface *
constSASS2SAI(const SASeqSrc *baseClass)
{
  return (const SuffixarrayFileInterface *)
    ((const char *)baseClass - offsetof(SuffixarrayFileInterface, baseClass));
}

static inline SASeqSrc *
SAI2SASS(SuffixarrayFileInterface *sai)
{
  return &sai->baseClass;
}

static inline const SASeqSrc *
constSAI2SASS(const SuffixarrayFileInterface *sai)
{
  return &sai->baseClass;
}

#endif
