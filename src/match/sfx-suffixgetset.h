/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_SUFFIXGETSET_H
#define SFX_SUFFIXGETSET_H

#include <stdio.h>
#include <inttypes.h>
#include "core/logger_api.h"
#include "core/error_api.h"
#include "core/encseq_api.h"
#include "core/encseq.h"
#include "core/bitbuffer.h"

#define GT_SUFFIXSORTSPACE_EXPORT_SET(SSSP,EXPORTPTR,INDEX,POS)\
        if ((EXPORTPTR)->ulongtabsectionptr != NULL)\
        {\
          (EXPORTPTR)->ulongtabsectionptr[INDEX] = POS;\
        } else\
        {\
          gt_assert((EXPORTPTR)->uinttabsectionptr != NULL &&\
                    (POS) <= (GtUword) UINT_MAX);\
          (EXPORTPTR)->uinttabsectionptr[INDEX] = (uint32_t) POS;\
        }\
        if (POS == 0)\
        {\
          gt_suffixsortspace_updatelongest(SSSP,INDEX);\
        }

typedef struct GtSuffixsortspace GtSuffixsortspace;

typedef struct
{
  uint32_t *uinttabsectionptr;
  GtUword *ulongtabsectionptr;
} GtSuffixsortspace_exportptr;

typedef void (*GtProcessunsortedsuffixrange)(void *,
                                             GtUword,
                                             GtUword,
                                             GtUword);

GtSuffixsortspace *gt_suffixsortspace_new(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint,
                                          GtLogger *logger);

void gt_suffixsortspace_delete(GtSuffixsortspace *suffixsortspace,
                               bool checklongestdefined);

void gt_suffixsortspace_showrange(const GtSuffixsortspace *sssp,
                                  GtUword subbucketleft,
                                  GtUword width);

GtUword gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,
                                     GtUword idx);

void gt_suffixsortspace_nooffsets(const GtSuffixsortspace *sssp);

void gt_suffixsortspace_updatelongest(GtSuffixsortspace *sssp,
                                      GtUword idx);

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                  GtUword idx,
                                  GtUword value);

void gt_suffixsortspace_init_seqstartpos(GtSuffixsortspace *sssp,
                                         const GtEncseq *encseq);

void gt_suffixsortspace_init_identity(GtSuffixsortspace *sssp,
                                      GtUword numofsuffixes);

GtSuffixsortspace_exportptr *gt_suffixsortspace_exportptr(
                                  GtSuffixsortspace *sssp,
                                  GtUword subbucketleft);

void gt_suffixsortspace_export_done(GtSuffixsortspace *sssp);

GtUword gt_suffixsortspace_get(const GtSuffixsortspace *sssp,
                               GtUword subbucketleft,
                               GtUword idx);

const GtUword *gt_suffixsortspace_getptr_ulong(const GtSuffixsortspace *sssp,
                                               GtUword subbucketleft);

const uint32_t *gt_suffixsortspace_getptr_uint32(const GtSuffixsortspace *sssp,
                                                 GtUword subbucketleft);

void gt_suffixsortspace_set(GtSuffixsortspace *sssp,
                            GtUword subbucketleft,
                            GtUword idx,
                            GtUword value);

GtUword gt_suffixsortspace_bucketleftidx_get(const GtSuffixsortspace *sssp);

void gt_suffixsortspace_bucketleftidx_set(GtSuffixsortspace *sssp,
                                          GtUword value);

void gt_suffixsortspace_sortspace_delete(GtSuffixsortspace *sssp);

void gt_suffixsortspace_partoffset_set(GtSuffixsortspace *sssp,
                                       GtUword partoffset);

const GtUword *gt_suffixsortspace_ulong_get(const GtSuffixsortspace *sssp);

GtUword gt_suffixsortspace_longest(const GtSuffixsortspace *sssp);

uint64_t gt_suffixsortspace_requiredspace(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint);

void gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                 const GtSuffixsortspace *sssp,
                                 GtUword numberofsuffixes);

void gt_suffixsortspace_compressed_to_file (const GtSuffixsortspace *sssp,
                                            GtBitbuffer *bb,
                                            GtUword numberofsuffixes);

typedef struct GtSSSPbuf GtSSSPbuf;

GtSSSPbuf *gt_SSSPbuf_new(GtUword size);

GtUword gt_SSSPbuf_filled(const GtSSSPbuf *sssp_buf);

void gt_SSSPbuf_delete(GtSSSPbuf *sssp_buf);

bool gt_SSSPbuf_fillspecialnextpage(GtSuffixsortspace *sssp,
                                    GtReadmode readmode,
                                    GtSpecialrangeiterator *sri,
                                    GtUword totallength,
                                    GtSSSPbuf *sssp_buf);

#endif
