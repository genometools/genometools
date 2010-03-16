/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
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

#ifndef REPEATS_H
#define REPEATS_H

#include "core/unused_api.h"
#include "match/querymatch.h"
#include "repeattypes.h"
#include "ltrharvest-opt.h"

void showrepeats (RepeatInfo *repeatinfo,
                  unsigned long seedminlength);

int simpleexactselfmatchstore(void *info,
                              const GtEncodedsequence *encseq,
                              Seqpos len,
                              Seqpos pos1,
                              Seqpos pos2,
                              GtError *err);

int subsimpleexactselfmatchstore(void *info,
                                 GT_UNUSED const GtEncodedsequence *encseq,
                                 const Querymatch *querymatch,
                                 GT_UNUSED GtError *err);

#endif
