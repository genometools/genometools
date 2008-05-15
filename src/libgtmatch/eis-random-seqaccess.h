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
#ifndef EIS_RANDOM_SEQACCESS_H
#define EIS_RANDOM_SEQACCESS_H

/**
 * @file eis-construction-interface.h
 * @brief Abstract interface for suffix array sources.
 *
 * The function pointer types defined here are essentially the
 * interface used for index construction: a source for suffix array
 * and BWT sequence entries is characterized by an interface
 * conforming to these functions.
 */
#include <stdlib.h>

#include "libgtcore/error.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-mrangealphabet.h"

/**
 * \brief generic method to access the original encoded sequence
 * @return actual number of symbols acquired
 */
typedef size_t (*accessSeqSubStr)(const void *state, Symbol *dest, Seqpos pos,
                                  size_t len);

struct randomSeqAccessor
{
  accessSeqSubStr accessFunc;
  void *state;
};

typedef struct randomSeqAccessor RandomSeqAccessor;

static inline size_t
accessSequence(RandomSeqAccessor accessor, Symbol *dest, Seqpos pos,
               size_t len)
{
  return accessor.accessFunc(accessor.state, dest, pos, len);
}

#endif
