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

#ifndef EMIMERGEESA_H
#define EMIMERGEESA_H

#include "seqpos-def.h"
#include "sarr-def.h"
#include "merger-trie.h"

#define SIZEOFMERGERESULTBUFFER BUFSIZ

/* make the type opaque */

typedef struct
{
  unsigned int idx;  /* index of genome in list of all genomes */
  Seqpos startpos;   /* in the range [0..totallength single index] */
} Indexedsuffix;

typedef struct
{
  unsigned int nextaccessidx,  /* in the range [0..SIZEOFMERGERESULTBUFFER] */
               nextstoreidx;   /* in the range [0..SIZEOFMERGERESULTBUFFER] */
  Seqpos lcptabstore[SIZEOFMERGERESULTBUFFER];
  Indexedsuffix suftabstore[SIZEOFMERGERESULTBUFFER];
  bool lastpage;
} Suflcpbuffer;

typedef struct
{
  uint64_t ident;              /* can be arbitrary large */
  unsigned int numofentries,   /* in the range [0..numofindexes-1] */
               numofindexes;   /* number of indexes */
  Seqpos *nextpostable;        /* in the range [0..totallength single index] */
  Suflcpbuffer buf;
  Mergertrierep trierep;
  Suffixarray *suffixarraytable;
  Alphabet *alpha;
} Emissionmergedesa;

#endif
