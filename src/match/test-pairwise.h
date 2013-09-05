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

#ifndef TEST_PAIRWISE_H
#define TEST_PAIRWISE_H

#include <stdbool.h>
#include "core/types_api.h"

typedef void (*Checkcmppairfuntype)(bool,
                                    const GtUchar *,GtUword,
                                    const GtUchar *,GtUword);

void gt_runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                               const char *file1,
                               const char *file2);

GtUword gt_runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                                     const char *text);

GtUword gt_applycheckfunctiontotext(const GtUchar *text,
                                       GtUword textlen,
                                       void *info);

GtUword gt_runcheckfunctiononalphalen(Checkcmppairfuntype checkfunction,
                                         const char *charlist,
                                         GtUword len);

GtUword gt_computegreedyunitedist(const GtUchar *useq,
                                        GtUword ulen,
                                        const GtUchar *vseq,
                                        GtUword vlen);

void gt_checkgreedyunitedist(bool forward,
                          const GtUchar *useq,
                          GtUword ulen,
                          const GtUchar *vseq,
                          GtUword vlen);

#endif
