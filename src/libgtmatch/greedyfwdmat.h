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
#include "libgtcore/strarray.h"
#include "libgtcore/error.h"
#include "defined-types.h"
#include "seqpos-def.h"
#include "ushort-def.h"
#include "alphadef.h"
#include "encseq-def.h"

typedef unsigned long (*Greedygmatchforwardfunction) (const void *,
                                                      unsigned long offset,
                                                      Seqpos left,
                                                      Seqpos right,
                                                      Seqpos *,
                                                      const Uchar *,
                                                      const Uchar *);

int findsubquerygmatchforward(const Encodedsequence *encseq,
                              const void *genericindex,
                              Seqpos totallength,
                              Greedygmatchforwardfunction gmatchforward,
                              const Alphabet *alphabet,
                              const StrArray *queryfilenames,
                              Definedunsignedlong minlength,
                              Definedunsignedlong maxlength,
                              bool showsequence,
                              bool showquerypos,
                              bool showsubjectpos,
                              Error *err);

int runsubstringiteration(Greedygmatchforwardfunction gmatchforward,
                          const void *genericindex,
                          Seqpos totalwidth,
                          const Seqpos *leftborder,
                          const Seqpos *countspecialcodes,
                          const Alphabet *alphabet,
                          unsigned int prefixlength,
                          const StrArray *queryfilenames,
                          Error *err);
#endif
