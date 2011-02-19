/*
  Copyright (c) 2008-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#ifndef SA_VISITOR_REP_H
#define SA_VISITOR_REP_H

#include "gth/sa_visitor.h"

/* the ``spliced alignment visitor'' interface */
struct GthSAVisitorClass {
  size_t size;
  void (*free)(GthSAVisitor*);
  void (*preface)(GthSAVisitor*);
  void (*visit_sa)(GthSAVisitor*, GthSA*);
  void (*trailer)(GthSAVisitor*, unsigned long num_of_sas);
};

struct GthSAVisitor {
  const GthSAVisitorClass *c_class;
};

GthSAVisitor* gth_sa_visitor_create(const GthSAVisitorClass*);
void*         gth_sa_visitor_cast(const GthSAVisitorClass*, GthSAVisitor*);

#endif
