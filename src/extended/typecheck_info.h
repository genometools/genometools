/*
  Copyright (c) 2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef TYPECHECK_INFO_H
#define TYPECHECK_INFO_H

#include "core/option_api.h"
#include "extended/type_checker_api.h"

typedef struct GtTypecheckInfo GtTypecheckInfo;

GtTypecheckInfo* gt_typecheck_info_new(void);
void             gt_typecheck_info_delete(GtTypecheckInfo*);
void             gt_typecheck_info_register_options(GtTypecheckInfo*,
                                                    GtOptionParser*);
bool             gt_typecheck_info_option_used(const GtTypecheckInfo*);
GtTypeChecker*   gt_typecheck_info_create_type_checker(const GtTypecheckInfo*,
                                                       GtError *err);

#endif
