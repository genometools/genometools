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

#ifndef SFX_APFXLEN_H
#define SFX_APFXLEN_H

#include "core/error_api.h"

#define GT_RECOMMENDED_MULTIPLIER_DEFAULT 0.25

unsigned int gt_recommendedprefixlength(unsigned int numofchars,
                                        unsigned long totallength,
                                        double recommendedmultiplier,
                                        bool withspecialsuffixes);

unsigned int gt_whatisthemaximalprefixlength(unsigned int numofchars,
                                          unsigned long totallength,
                                          unsigned int prefixlenbits,
                                          bool withspecialsuffixes);

int gt_checkprefixlength(unsigned int maxprefixlen,
                         unsigned int prefixlength,
                         GtError *err);

void gt_showmaximalprefixlength(GtLogger *logger,
                                unsigned int maxprefixlen,
                                unsigned int recommended);

#endif
