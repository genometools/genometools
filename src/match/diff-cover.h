/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef DIFF_COVER_H
#define DIFF_COVER_H

#include "seqpos-def.h"
#include "encodedsequence.h"
#include "readmode-def.h"
#include "core/logger.h"

typedef struct Differencecover Differencecover;

void differencecovers_check(const GtEncodedsequence *encseq,Readmode readmode);

Differencecover *differencecover_new(unsigned int vparam,
                                     const GtEncodedsequence *encseq,
                                     Readmode readmode,
                                     GtLogger *logger);

int differencecover_vparamverify(const Differencecover *dcov,GtError *err);

void differencecover_sortsample(Differencecover *dcov,bool cmpcharbychar,
                                bool withcheck);

void differencecover_delete(Differencecover *dcov);

void dc_sortunsortedbucket(void *data,
                           Seqpos *left,
                           Seqpos *right,
                           Seqpos depth);

#endif
