/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "core/error.h"
#include "core/str.h"
#include "match/sarr-def.h"
#include "match/spacedef.h"
#include "match/esa-seqread.h"
#include "match/intcode-def.h"

#include "ltrharvest-opt.h"
#include "repeats.h"
#include "searchforLTRs.h"
#include "duplicates.h"
#include "outputstd.h"
#include "outputfasta.h"
#include "outputgff3.h"

#include "match/esa-maxpairs.pr"

static int runltrharvest(LTRharvestoptions *lo, GtError *err)
{
  Sequentialsuffixarrayreader *ssar; /* suffix array */
  bool had_err = false;
  const Encodedsequence *encseq;

  gt_error_check(err);

  ssar = newSequentialsuffixarrayreaderfromfile(lo->str_indexname,
                                  SARR_LCPTAB | SARR_SUFTAB |
                                  SARR_ESQTAB | SARR_DESTAB |
                                  SARR_SSPTAB,
                                  SEQ_mappedboth,
                                  err);
  if (ssar == NULL)
  {
    return -1;
  }
  encseq = encseqSequentialsuffixarrayreader(ssar);

  /* test if motif is valid and encode motif */
  if (testmotifandencodemotif (&lo->motif, encseq, err) != 0)
  {
    had_err = true;
  }

  /* show defined option and values */
  if (!had_err && lo->verbosemode)
  {
    showuserdefinedoptionsandvalues(lo);
  }

  lo->markpos = getencseqssptab(encseq);

  /* init array for maximal repeats */
  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.ssarptr = ssar;

  /* search for maximal repeats */
  if (!had_err && enumeratemaxpairs(ssar,encseq,
                                    readmodeSequentialsuffixarrayreader(ssar),
                                    (unsigned int)lo->minseedlength,
                                    (void*) simpleexactselfmatchstore,
                                    lo,
                                    NULL,
                                    err) != 0)
  {
    had_err = true;
  }

  /* init array for candidate pairs */
  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if (!had_err && searchforLTRs (ssar, lo, lo->markpos, err) != 0)
  {
    had_err = true;
  }

  /* free array for maximal repeats */
  FREEARRAY(&lo->repeatinfo.repeats, Repeat);

  /* remove exact duplicates */
  if (!had_err)
  {
    removeduplicates(&lo->arrayLTRboundaries);
  }

  /* remove overlapping predictions if desired */
  if (!had_err && (lo->nooverlapallowed || lo->bestofoverlap))
  {
    removeoverlapswithlowersimilarity(&lo->arrayLTRboundaries,
                                      lo->nooverlapallowed);
  }

  /* print multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutput)
  {
    if (showpredictionsmultiplefasta(lo,
                                     false,
                                     60U,
                                     ssar,
                                     true,
                                     err) != 0)
    {
      had_err = true;
    }
  }

  /* print inner region multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutputinnerregion)
  {
    if (showpredictionsmultiplefasta(lo,
                                     true,
                                     60U,
                                     ssar,
                                     true,
                                     err) != 0)
    {
      had_err = true;
    }
  }

  /* print GFF3 format file of predictions */
  if (!had_err && lo->gff3output)
  {
    if (printgff3format(lo, encseq,err) != 0)
    {
      had_err = true;
    }
  }

  /* print predictions to stdout */
  if (!had_err)
  {
    if (showinfoiffoundfullLTRs(lo, ssar) != 0)
    {
      had_err = true;
    }
  }

  FREEARRAY(&lo->arrayLTRboundaries, LTRboundaries);
  freeSequentialsuffixarrayreader(&ssar);

  return had_err ? -1 : 0;
}

int parseargsandcallltrharvest(int argc,const char *argv[],GtError *err)
{
  LTRharvestoptions lo;
  int had_err = 0;

  if (ltrharvestoptions(&lo,argc,argv,err) != 0)
  {
    had_err = -1;
  } else
  {
    printargsline(argv,argc);
    had_err = runltrharvest(&lo,err);
  }
  wrapltrharvestoptions(&lo);
  return had_err;
}
