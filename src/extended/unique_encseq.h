/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef UNIQUE_ENCSEQ_H
#define UNIQUE_ENCSEQ_H

#include "core/encseq_api.h"
#include "core/hashmap_api.h"
#include "unique_encseq_rep.h"
#include "core/array_api.h"
#include "match/xdrop.h"
#include "match/greedyedist.h"
#include "extended/alignment.h"
#include "core/seq.h"
#include "core/encseq.h"
#include "extended/swalign.h"

typedef struct GtUniqueEncseq GtUniqueEncseq;

GtUniqueEncseq *gt_unique_encseq_new(void);
void gt_unique_encseq_delete(GtUniqueEncseq *unique_encseq);

GtUniqueEncseqDB *gt_unique_encseq_new_db(GtEncseq *encseq);
void gt_unique_encseq_delete_db(GtUniqueEncseqDB *uedb);

void gt_unique_encseq_database_stats_coarse(GtUniqueEncseqDB *uedb, FILE *fp);
void gt_unique_encseq_database_stats_fine(GtUniqueEncseqDB *uedb, FILE *fp);
int  gt_unique_encseq_check_db(GtUniqueEncseqDB *uedb,
     GtLogger *debug_logger,
     GtError *err);
int  gt_unique_encseq_get_sequence_from_range(GtRange *range,
     GtEncseq *unique_encseq, GtUniqueEncseqDB *uedb, FILE *fp, GtError *err);
int  gt_unique_encseq_get_sequence_from_idx(GtUword idx,
     GtEncseq *unique_encseq, GtUniqueEncseqDB *uedb, FILE *fp, GtError *err);

void gt_unique_encseq_print_editscripts(GtUniqueEncseqDB *uedb,
    GtEncseq *unique_encseq, FILE *fp);
void gt_unique_encseq_print_db(GtUniqueEncseqDB *uedb);

GtUniqueEncseqDB *gt_unique_encseq_uedb_read(FILE *fp);
int gt_unique_encseq_encseq2uniqueencseq(GtUniqueEncseqDB *uedb,
    const GtEncseq *encseq, const char *indexname,
    GtError *err);

void getencseqkmers_only_regular_kmers2(GtUniqueEncseqInfo *ueinfo);

#endif
