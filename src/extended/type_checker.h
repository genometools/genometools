/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TYPE_CHECKER_H
#define TYPE_CHECKER_H

/* The <GtTypeChecker> interface. */

typedef struct GtTypeCheckerClass GtTypeCheckerClass;
typedef struct GtTypeChecker GtTypeChecker;

#include "core/str_array.h"
#include "extended/feature_type.h"

/* Return a new reference to <type_checker>. */
GtTypeChecker* gt_type_checker_ref(GtTypeChecker *type_checker);
/* Returns <true> if <type> is a valid type, <false> otherwise. */
bool           gt_type_checker_is_valid(GtTypeChecker*, const char *type);
void           gt_type_checker_delete(GtTypeChecker*);

#endif
