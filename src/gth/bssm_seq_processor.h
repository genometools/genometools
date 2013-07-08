/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#ifndef BSSM_SEQ_PROCESSOR_H
#define BSSM_SEQ_PROCESSOR_H

#include "core/file.h"
#include "extended/region_mapping_api.h"

typedef struct GthBSSMSeqProcessor GthBSSMSeqProcessor;

GthBSSMSeqProcessor* gth_bssm_seq_processor_new(const char *outdir,
                                                GtFileMode file_mode,
                                                bool force, bool gcdonor,
                                                GtError *err);
void                 gth_bssm_seq_processor_delete(GthBSSMSeqProcessor*);
void                 gth_bssm_seq_processor_proc_exon(GthBSSMSeqProcessor*,
                                                      unsigned int phase,
                                                      GtStr *seqid,
                                                      const GtRange *range,
                                                      bool reverse,
                                                      const GtStr *seq);
void                 gth_bssm_seq_processor_proc_intron(GthBSSMSeqProcessor*,
                                                        unsigned int phase,
                                                        GtStr *seqid,
                                                        const GtRange *range,
                                                        bool reverse,
                                                        const GtStr *seq);
void                 gth_bssm_seq_processor_squash(GthBSSMSeqProcessor*);
int                  gth_bssm_seq_processor_find_true_sites(
                                                           GthBSSMSeqProcessor*,
                                                            GtRegionMapping
                                                            *region_mapping,
                                                            GtError *err);
int                  gth_bssm_seq_processor_find_false_sites(
                                                           GthBSSMSeqProcessor*,
                                                             GtRegionMapping
                                                             *region_mapping,
                                                             GtError *err);
int                  gth_bssm_seq_processor_write_intermediate(
                                                          GthBSSMSeqProcessor*,
                                                               GtError *err);
void                 gth_bssm_seq_processor_sample(GthBSSMSeqProcessor*,
                                                   bool verbose, GtFile *logfp);
void                 gth_bssm_seq_processor_write(GthBSSMSeqProcessor*);

#endif
