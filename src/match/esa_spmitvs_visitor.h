/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef ESA_SPMITVS_VISITOR_H
#define ESA_SPMITVS_VISITOR_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/logger.h"
#include "match/esa_visitor.h"

typedef struct GtESASpmitvsVisitor GtESASpmitvsVisitor;

const GtESAVisitorClass* gt_esa_spmitvs_visitor_class(void);
GtESAVisitor*            gt_esa_spmitvs_visitor_new(const GtEncseq *encseq,
                                                    GtReadmode readmode,
                                                    unsigned int prefixlength,
                                                    GtError *err);
void                     gt_esa_spmitvs_visitor_print_results(
                                                     GtESASpmitvsVisitor*,
                                                     unsigned long nonspecials);

#endif
