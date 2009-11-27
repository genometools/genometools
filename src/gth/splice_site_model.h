/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICE_SITE_MODEL_H
#define SPLICE_SITE_MODEL_H

#include "core/error.h"

/* generic splice site model */
typedef struct GthSpliceSiteModel GthSpliceSiteModel;

GthSpliceSiteModel* gth_splice_site_model_new(void);
void                gth_splice_site_model_delete(GthSpliceSiteModel*);
int                 gth_splice_site_model_load_bssm(GthSpliceSiteModel*,
                                                    const char *bssmfile,
                                                    GtError*);
void                gth_splice_site_model_U12intronmodel_set_usage(
                                                        GthSpliceSiteModel*,
                                                        bool useU12intronmodel);
void                gth_splice_site_model_set_U12typedonorprob(
                                                        GthSpliceSiteModel*,
                                                        GthFlt prob);
void                gth_splice_site_model_set_U12typedonorprob_one_mismatch(
                                                        GthSpliceSiteModel*,
                                                        GthFlt prob);

#endif
