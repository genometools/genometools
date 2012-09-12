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

#include "sfx-sais.h"
#include "core/intbits.h"

struct GtSainlabels
{
  GtBitsequence *starseq;
};

#define GT_SSTARLENGTH_MAX 100

static void showSstarlendist(const unsigned long *lendist,
                             unsigned long longerthanmax,
                             unsigned long countSstar)
{
  unsigned long idx;

  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,lendist[idx],
                                (double) lendist[idx]/countSstar);
    }
  }
  if (longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,longerthanmax,
                              (double) longerthanmax/countSstar);
  }
}

GtSainlabels *gt_sain_labels_new(const GtEncseq *encseq)
{
  unsigned long position,
                countS = 0,
                countSstar = 0,
                nextcc = GT_UNIQUEINT(SEPARATOR),
                nextSstartypepos,
                totalSstarlength = 0,
                lendist[GT_SSTARLENGTH_MAX+1] = {0},
                longerthanmax = 0,
                totallength = gt_encseq_total_length(encseq);
  bool nextisStype = true;
  GtEncseqReader *esr;
  GtSainlabels *sainlabels;

  esr = gt_encseq_create_reader_with_readmode(encseq,GT_READMODE_REVERSE,0);
  nextSstartypepos = totallength;
  sainlabels = gt_malloc(sizeof *sainlabels);
  GT_INITBITTAB(sainlabels->starseq,totallength+1);
  for (position = totallength-1; /* Nothing */; position--)
  {
    GtUchar cc = gt_encseq_reader_next_encoded_char(esr);
    bool currentisStype;
    unsigned long currentcc
      = ISSPECIAL(cc) ? GT_UNIQUEINT(cc) : (unsigned long) cc;

    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      countS++;
      currentisStype = true;
    } else
    {
      currentisStype = false;
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      countSstar++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position + 1;
      totalSstarlength += currentlen;
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        lendist[currentlen]++;
      } else
      {
        longerthanmax++;
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
  printf("S-type: %lu (%.2f)\n",countS,(double) countS/totallength);
  printf("Sstar-type: %lu (%.2f)\n",countSstar,(double) countSstar/countS);
  printf("Sstar-type.length: %lu (%.2f)\n",totalSstarlength,
                   (double) totalSstarlength/countSstar);
  showSstarlendist(lendist,longerthanmax,countSstar);
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
