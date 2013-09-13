/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "spmsuftab.h"

GtSpmsuftab *gt_spmsuftab_new(GtUword numofentries,
                              GtUword maxvalue,
                              unsigned int bitsforseqnumrelpos,
                              GtLogger *logger)
{
  GtSpmsuftab *spmsuftab;
  unsigned int bitsforpositions;
  GtUword required;

  spmsuftab = gt_malloc(sizeof (*spmsuftab));
  required = (GtUword) gt_spmsuftab_requiredspace(numofentries,
                                                        maxvalue,
                                                        bitsforseqnumrelpos);
  bitsforpositions = gt_determinebitspervalue(maxvalue);
  if (bitsforpositions < bitsforseqnumrelpos)
  {
    gt_logger_log(logger,"use "GT_WU" bitpackarray-entries for all positions "
                         "(%u bits each,"GT_WU" bytes total)",
                         numofentries,bitsforpositions,required);
    spmsuftab->bitpackarray
      = gt_compact_ulong_store_new(numofentries,bitsforpositions);
    spmsuftab->usebitsforpositions = true;
    spmsuftab->maxvalue = maxvalue;
  } else
  {
    gt_logger_log(logger,"use "GT_WU" bitpackarray-entries for all "
                  "seqnum/relpos-pairs (%u bits each,"GT_WU" bytes total)",
                  numofentries,bitsforseqnumrelpos,required);
    spmsuftab->bitpackarray
      = gt_compact_ulong_store_new(numofentries,bitsforseqnumrelpos);
    spmsuftab->usebitsforpositions = false;
    spmsuftab->maxvalue = (1UL << bitsforseqnumrelpos) - 1;
  }
  spmsuftab->partoffset = 0;
  spmsuftab->numofentries = numofentries;
  return spmsuftab;
}

void gt_spmsuftab_delete(GtSpmsuftab *spmsuftab)
{
  if (spmsuftab != NULL)
  {
    gt_compact_ulong_store_delete(spmsuftab->bitpackarray);
    gt_free(spmsuftab);
  }
}

bool gt_spmsuftab_usebitsforpositions(const GtSpmsuftab *spmsuftab)
{
  return spmsuftab->usebitsforpositions;
}

void gt_spmsuftab_partoffset(GtSpmsuftab *spmsuftab,GtUword offset)
{
  spmsuftab->partoffset = offset;
}

size_t gt_spmsuftab_requiredspace(GtUword numofentries,
                                  GtUword maxvalue,
                                  unsigned int bitsforseqnumrelpos)
{
  unsigned int bitsforpositions = gt_determinebitspervalue(maxvalue);

  return sizeof (GtSpmsuftab) +
         gt_compact_ulong_store_size(numofentries,
                                     MIN(bitsforpositions,bitsforseqnumrelpos));
}
