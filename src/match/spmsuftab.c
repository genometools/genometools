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
  GtSpmsuftab *spmsuftab = gt_malloc(sizeof (*spmsuftab));
  unsigned int bitspervalue = gt_determinebitspervalue((uint64_t) maxvalue);
  unsigned long required = (unsigned long)
                           gt_spmsuftab_requiredspace(numofentries,maxvalue);

  if (bitspervalue > 32U)
  {
    gt_logger_log(logger,"use %lu bitpackarray-entries of %u bits (%lu bytes)",
                  numofentries,bitspervalue,required);
    spmsuftab->bitpackarray
#ifdef THOMASBITPACK
      = bitpackarray_new(bitspervalue,(BitOffset) numofentries,true);
    gt_logger_log(logger,"thomas' version");
#else
      = gt_GtCompactulongstore_new(numofentries,bitspervalue);
    gt_logger_log(logger,"sk version");
#endif
    spmsuftab->suftab = NULL;
  } else
  {
    gt_logger_log(logger,"use %lu uint32_t-entries (%lu bytes)",
                         numofentries,required);
    spmsuftab->bitpackarray = NULL;
    spmsuftab->suftab = gt_malloc(sizeof (*spmsuftab->suftab) * numofentries);
  }
  spmsuftab->partoffset = 0;
  spmsuftab->numofentries = numofentries;
  spmsuftab->maxvalue = maxvalue;
  return spmsuftab;
}

void gt_spmsuftab_delete(GtSpmsuftab *spmsuftab)
{
  if (spmsuftab->bitpackarray != NULL)
  {
    gt_free(spmsuftab->suftab);
#ifdef THOMASBITPACK
    bitpackarray_delete(spmsuftab->bitpackarray);
#else
    gt_GtCompactulongstore_delete(spmsuftab->bitpackarray);
#endif
  } else
  {
    gt_free(spmsuftab->suftab);
  }
  gt_free(spmsuftab);
}

void gt_spmsuftab_partoffset(GtSpmsuftab *spmsuftab,unsigned long offset)
{
  spmsuftab->partoffset = offset;
}

size_t gt_spmsuftab_requiredspace(unsigned long numofentries,
                                  GT_UNUSED unsigned long maxvalue)
{
  unsigned int bitspervalue = gt_determinebitspervalue((uint64_t) maxvalue);

  if (bitspervalue > 32U)
  {
    return sizeof (GtSpmsuftab) +
#ifdef THOMASBITPACK
           sizeofbitarray(bitspervalue,(BitOffset) numofentries);
#else
           gt_GtCompactulongstore_size(numofentries,bitspervalue);
#endif
  } else
  {
    return sizeof (GtSpmsuftab) + numofentries * sizeof (uint32_t);
  }
}
