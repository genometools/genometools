/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef READS_LIBRARIES_TABLE_H
#define READS_LIBRARIES_TABLE_H

#include "core/error_api.h"

typedef struct GtReadsLibrariesTable GtReadsLibrariesTable;

GtReadsLibrariesTable* gt_reads_libraries_table_new(GtUword noflibraries);

void                   gt_reads_libraries_table_delete(
                                            GtReadsLibrariesTable *rlt);

void                   gt_reads_libraries_table_add(GtReadsLibrariesTable *rlt,
                                                    GtUword first_seqnum,
                                                    GtUword insertlength,
                                                    GtUword stdev,
                                                    bool paired);

void                   gt_reads_libraries_table_get_library(
                                            GtReadsLibrariesTable *rlt,
                                            GtUword libnum,
                                            GtUword *first_seqnum,
                                            GtUword *insertlength,
                                            GtUword *stdev);

void                   gt_reads_libraries_table_save(GtReadsLibrariesTable *rlt,
                                                     FILE *rlt_fp);

GtUword          gt_reads_libraries_table_noflibraries(
                                            GtReadsLibrariesTable *rlt);

GtUword          gt_reads_libraries_table_firstunpaired(
                                            GtReadsLibrariesTable *rlt);

GtReadsLibrariesTable* gt_reads_libraries_table_load(FILE *rlt_fp,
                                                     GtError *err);

#endif
