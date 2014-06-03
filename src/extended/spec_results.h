/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#ifndef SPEC_RESULTS_H
#define SPEC_RESULTS_H

typedef struct GtSpecResults GtSpecResults;

typedef enum {
  GT_SPEC_SUCCESS,
  GT_SPEC_FAILURE,
  GT_SPEC_RUNTIME_ERROR
} GtSpecResultStatus;

#include "core/file_api.h"

GtSpecResults* gt_spec_results_new(void);
void           gt_spec_results_add_result(GtSpecResults *sr,
                                          const char *aspect,
                                          GtGenomeNode *node,
                                          GtSpecResultStatus status,
                                          const char *error_string);
void           gt_spec_results_add_cc(GtSpecResults *sr);
void           gt_spec_results_report(GtSpecResults *sr, GtFile *outfile,
                                      const char *specfile, bool details,
                                      bool colored);
void           gt_spec_results_delete(GtSpecResults *spec_results);

#endif
