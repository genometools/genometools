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

#ifndef GREEDYFWDMAT_H
#define GREEDYFWDMAT_H

#include <stdbool.h>
#include "core/alphabet.h"
#include "core/error.h"
#include "core/types_api.h"
#include "core/str_array.h"
#include "core/defined-types.h"

#include "core/encseq.h"

typedef unsigned long (*Greedygmatchforwardfunction) (const void *,
                                                      unsigned long offset,
                                                      unsigned long left,
                                                      unsigned long right,
                                                      unsigned long *,
                                                      const GtUchar *,
                                                      const GtUchar *);

int gt_findsubquerygmatchforward(const GtEncseq *encseq,
                              const void *genericindex,
                              unsigned long totallength,
                              Greedygmatchforwardfunction gmatchforward,
                              const GtAlphabet *alphabet,
                              const GtStrArray *queryfilenames,
                              Definedunsignedlong minlength,
                              Definedunsignedlong maxlength,
                              bool showsequence,
                              bool showquerypos,
                              bool showsubjectpos,
                              GtError *err);

int runsubstringiteration(Greedygmatchforwardfunction gmatchforward,
                          const void *genericindex,
                          unsigned long totalwidth,
                          const unsigned long *leftborder,
                          const unsigned long *countspecialcodes,
                          const GtAlphabet *alphabet,
                          unsigned int prefixlength,
                          const GtStrArray *queryfilenames,
                          GtError *err);
#endif
