/*
  Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2013 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRDIGEST_FILE_OUT_STREAM_H
#define LTRDIGEST_FILE_OUT_STREAM_H

#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"
#include "ltr/ltrdigest_def.h"

/* implements the ``node_stream'' interface */
typedef struct GtLTRdigestFileOutStream GtLTRdigestFileOutStream;

const GtNodeStreamClass* gt_ltrdigest_file_out_stream_class(void);

GtNodeStream* gt_ltrdigest_file_out_stream_new(GtNodeStream *in_stream,
                                        int tests_to_run,
                                        GtRegionMapping *rmap,
                                        char *file_prefix,
                                        unsigned int seqnamelen,
                                        GtError *err);

int           gt_ltrdigest_file_out_stream_write_metadata(
                                              GtLTRdigestFileOutStream *ls,
                                              int tests_to_run,
                                              const char *trnafilename,
                                              const char *gfffilename,
                                              GtRange ppt_len,
                                              GtRange ubox_len,
                                              unsigned int ppt_radius,
                                              GtRange alilen,
                                              unsigned int max_edist,
                                              GtRange offsetlen,
                                              GtRange trnaoffsetlen,
                                              unsigned int pbs_radius,
                                              GtStrArray *hmm_files,
                                              unsigned int chain_max_gap_length,
                                              double evalue_cutoff,
                                              GtError *err);

void          gt_ltrdigest_file_out_stream_enable_pdom_alignment_output(
                                                                 GtNodeStream*);
void          gt_ltrdigest_file_out_stream_enable_aa_sequence_output(
                                                                 GtNodeStream*);

#endif
