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
#include "spmsuftab.h"

GtSpmsuftab *gt_spmsuftab_new(unsigned long numofentries,
                              unsigned long maxvalue,
                              GtLogger *logger)
{
  GtSpmsuftab *spmsuftab;
  unsigned int bitspervalue;
  unsigned long required = (unsigned long)
                           gt_spmsuftab_requiredspace(numofentries,maxvalue);

  bitspervalue = gt_determinebitspervalue((uint64_t) maxvalue);
  gt_logger_log(logger,"use %lu bitpackarray-entries of %u bits (%lu bytes)",
                numofentries,bitspervalue,required);
  spmsuftab = gt_malloc(sizeof (*spmsuftab));
  spmsuftab->bitpackarray
    = gt_GtCompactulongstore_new(numofentries,bitspervalue);
  spmsuftab->partoffset = 0;
  spmsuftab->numofentries = numofentries;
  spmsuftab->maxvalue = maxvalue;
  return spmsuftab;
}

void gt_spmsuftab_delete(GtSpmsuftab *spmsuftab)
{
  if (spmsuftab != NULL)
  {
    gt_GtCompactulongstore_delete(spmsuftab->bitpackarray);
    gt_free(spmsuftab);
  }
}

void gt_spmsuftab_partoffset(GtSpmsuftab *spmsuftab,unsigned long offset)
{
  spmsuftab->partoffset = offset;
}

size_t gt_spmsuftab_requiredspace(unsigned long numofentries,
                                  unsigned long maxvalue)
{
  unsigned int bitspervalue = gt_determinebitspervalue((uint64_t) maxvalue);

  return sizeof (GtSpmsuftab) +
         gt_GtCompactulongstore_size(numofentries,bitspervalue);
}
