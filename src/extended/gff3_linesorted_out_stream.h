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

#ifndef GFF3_LINESORTED_OUT_STREAM_H
#define GFF3_LINESORTED_OUT_STREAM_H

#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

/* implements the ``genome stream'' interface */
typedef struct GtGFF3LinesortedOutStream GtGFF3LinesortedOutStream;

const GtNodeStreamClass* gt_gff3_linesorted_out_stream_class(void);
GtNodeStream* gt_gff3_linesorted_out_stream_new(GtNodeStream *in_stream,
                                                GtFile *outfp);
/* Set the width with which the FASTA sequences of <GtSequenceNode>s passed
   through <gff3_out_stream> are shown to <fasta_width>.
   Per default, each FASTA entry is shown on a single line. */
void          gt_gff3_linesorted_out_stream_set_fasta_width(
                                     GtGFF3LinesortedOutStream *gff3_out_stream,
                                     GtUword fasta_width);
/* If this method is called upon <gff3_out_stream>, use the original ID
   attributes provided in the input (instead of creating new ones, which
   is the default). Memory consumption for <gff3_out_stream> is raised from O(1)
   to O(<input_size>), because bookkeeping of used IDs becomes necessary to
   avoid ID collisions. */
void          gt_gff3_linesorted_out_stream_retain_id_attributes(
                                    GtGFF3LinesortedOutStream *gff3_out_stream);

#endif
