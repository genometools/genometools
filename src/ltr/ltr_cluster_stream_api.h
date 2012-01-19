/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
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

#ifndef LTR_CLUSTER_STREAM_API_H
#define LTR_CLUSTER_STREAM_API_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "extended/node_stream_api.h"

typedef struct GtLTRClusterStream GtLTRClusterStream;

/* Implements the <GtNodeStream> interface. <GtLTRClusterStream> first extracts
   sequences for all features within <GtFeatureNode>. After that
   <GtMatchIteratorBlast> with BLASTN-process is used to match the sequences for
   each feature against themself. These matches are stored in <GtMatchBlast>.
   The match information is used to perform a single linkage clustering based
   on overlapping of the involved sequences. Finally the features are annotated
   with clusterids. */
GtNodeStream* gt_ltr_cluster_stream_new(GtNodeStream *in_stream,
                                        GtEncseq *encseq,
                                        GtStr *file_prefix,
                                        unsigned long plarge,
                                        unsigned long psmall,
                                        double evalue,
                                        bool dust,
                                        int word_size,
                                        int gapopen,
                                        int gapextend,
                                        int penalty,
                                        int reward,
                                        int num_threads,
                                        double xdrop,
                                        double identity,
                                        const char *moreblast,
                                        bool from_file,
                                        char **current_state,
                                        GtError *err);

#endif
