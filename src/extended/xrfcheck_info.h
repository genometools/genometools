/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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

#ifndef XRFCHECK_INFO_H
#define XRFCHECK_INFO_H

#include "core/option_api.h"
#include "extended/xrf_checker_api.h"

typedef struct GtXRFCheckInfo GtXRFCheckInfo;

GtXRFCheckInfo* gt_xrfcheck_info_new(void);
void            gt_xrfcheck_info_register_options(GtXRFCheckInfo*,
                                                    GtOptionParser*);
bool            gt_xrfcheck_info_option_used(const GtXRFCheckInfo*);
GtXRFChecker*   gt_xrfcheck_info_create_xrf_checker(const GtXRFCheckInfo*,
                                                    GtError *err);
void            gt_xrfcheck_info_delete(GtXRFCheckInfo*);

#endif
