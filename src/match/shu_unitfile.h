/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#ifndef SHU_UNITFILE_H
#define SHU_UNITFILE_H

#include "core/encseq_api.h"
#include "core/str_array_api.h"
#include "core/unused_api.h"

typedef struct {
  unsigned long num_of_genomes,
                num_of_files;
  GtStrArray *genome_names;
  const GtStrArray *file_names;
  /*array holding the mapping of file to genome*/
  unsigned long *map_files;
} GtShuUnitFileInfo;

/*
  reads the Unitfile and collects the names of the Genomes, and mapping of
  file number to genome, file_names and num_of_files has to be set befor
  calling this function. will allocate space for genome_names and map_files
  use gt_unit_info_delete to free unit_info.
  sets err and returns 1 on error
*/
int gt_read_genomediff_unitfile(const GtStr *unitfile,
                                GtShuUnitFileInfo *unit_info,
                                GT_UNUSED GtLogger *logger,
                                GtError *err);

void gt_shu_unit_info_files_as_units(GtShuUnitFileInfo *unit_info);

GtShuUnitFileInfo *gt_shu_unit_info_new(const GtEncseq *encseq);

/*
  frees memory of unit_info struct
*/
void gt_shu_unit_info_delete(GtShuUnitFileInfo *unit_info);

#endif
