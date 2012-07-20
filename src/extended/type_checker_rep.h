/*
  Copyright (c) 2008, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008       Center for Bioinformatics, University of Hamburg

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

#ifndef TYPE_CHECKER_REP_H
#define TYPE_CHECKER_REP_H

#include <stdio.h>
#include "core/cstr_table.h"
#include "extended/feature_type.h"
#include "extended/type_checker.h"

struct GtTypeCheckerClass {
  size_t size;
  const char* (*description)(GtTypeChecker*);
  bool        (*is_valid)(GtTypeChecker*, const char *type);
  bool        (*is_partof)(GtTypeChecker*, const char *parent_type,
                           const char *child_type);
  void        (*free)(GtTypeChecker*);
};

struct GtTypeChecker {
  const GtTypeCheckerClass *c_class;
  unsigned int reference_count;
};

GtTypeChecker* gt_type_checker_create(const GtTypeCheckerClass*);
void*          gt_type_checker_cast(const GtTypeCheckerClass*, GtTypeChecker*);

#endif
