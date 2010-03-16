/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
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

#ifndef ENCSEQ_SPECIALSRANK_H
#define ENCSEQ_SPECIALSRANK_H

#include "match/encodedsequence.h"
#include "match/seqpos-def.h"

typedef struct specialsRankLookup SpecialsRankLookup;

extern SpecialsRankLookup *
newSpecialsRankLookup(const GtEncodedsequence *encseq, Readmode readmode,
                     unsigned sampleIntervalLog2);

extern void
deleteSpecialsRankLookup(SpecialsRankLookup *table);

static inline Seqpos
specialsRank(const SpecialsRankLookup *rankTable, Seqpos pos);

extern const GtEncodedsequence *
SPRTGetOrigEncseq(const SpecialsRankLookup *rankTable);

#include "match/encseq-specialsrank-priv.h"

#endif
