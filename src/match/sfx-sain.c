/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "sfx-sain.h"
#include "core/intbits.h"

#define GT_SSTARLENGTH_MAX 100

struct GtSainlabels
{
  GtBitsequence *starseq;
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                totallength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
};

void gt_sain_labels_show(const GtSainlabels *sainlabels)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",sainlabels->countStype,
                      (double) sainlabels->countStype/sainlabels->totallength);
  printf("Sstar-type: %lu (%.2f)\n",sainlabels->countSstartype,
                    (double) sainlabels->countSstartype/sainlabels->countStype);
  printf("Sstar-type.length: %lu (%.2f)\n",sainlabels->totalSstarlength,
              (double) sainlabels->totalSstarlength/sainlabels->countSstartype);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (sainlabels->lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,sainlabels->lendist[idx],
               (double) sainlabels->lendist[idx]/sainlabels->countSstartype);
    }
  }
  if (sainlabels->longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,sainlabels->longerthanmax,
                 (double) sainlabels->longerthanmax/sainlabels->countSstartype);
  }
}

GtSainlabels *gt_sain_labels_new(const GtEncseq *encseq)
{
  unsigned long position,
                idx,
                nextcc = GT_UNIQUEINT(SEPARATOR),
                nextSstartypepos;
  bool nextisStype = true;
  GtEncseqReader *esr;
  GtSainlabels *sainlabels;

  esr = gt_encseq_create_reader_with_readmode(encseq,GT_READMODE_REVERSE,0);
  sainlabels = gt_malloc(sizeof *sainlabels);
  sainlabels->totalSstarlength = 0;
  sainlabels->countStype = 0;
  sainlabels->countSstartype = 0;
  sainlabels->longerthanmax = 0;
  sainlabels->totallength = gt_encseq_total_length(encseq);
  GT_INITBITTAB(sainlabels->starseq,sainlabels->totallength+1);
  nextSstartypepos = sainlabels->totallength;
  for (idx = 0; idx<=(unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    sainlabels->lendist[idx] = 0;
  }
  for (position = sainlabels->totallength-1; /* Nothing */; position--)
  {
    GtUchar cc = gt_encseq_reader_next_encoded_char(esr);
    bool currentisStype;
    unsigned long currentcc
      = ISSPECIAL(cc) ? GT_UNIQUEINT(cc) : (unsigned long) cc;

    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      sainlabels->countStype++;
      currentisStype = true;
    } else
    {
      currentisStype = false;
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      sainlabels->countSstartype++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position + 1;
      sainlabels->totalSstarlength += currentlen;
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        sainlabels->lendist[currentlen]++;
      } else
      {
        sainlabels->longerthanmax++;
      }
      nextSstartypepos = position;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_encseq_reader_delete(esr);
  return sainlabels;
}

void gt_sain_labels_delete(GtSainlabels *sainlabels)
{
  if (sainlabels != NULL)
  {
    gt_free(sainlabels->starseq);
    gt_free(sainlabels);
  }
}
