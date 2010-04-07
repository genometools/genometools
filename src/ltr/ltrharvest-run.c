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
#include "match/esa-seqread.h"
#include "match/esa-maxpairs.h"

#include "ltrharvest-opt.h"
#include "ltrharvest-run.h"
#include "repeats.h"
#include "searchforLTRs.h"
#include "duplicates.h"
#include "outputstd.h"
#include "outputfasta.h"
#include "outputgff3.h"

static int bdptrcompare(const void *a, const void *b)
{
  const LTRboundaries **bda, **bdb;

  bda = (const LTRboundaries **) a;
  bdb = (const LTRboundaries **) b;
  if ((*bda)->contignumber < (*bdb)->contignumber)
  {
    return -1;
  }
  if ((*bda)->contignumber > (*bdb)->contignumber)
  {
    return 1;
  }
  if ((*bda)->leftLTR_5 < (*bdb)->leftLTR_5)
  {
    return -1;
  }
  if ((*bda)->leftLTR_5 > (*bdb)->leftLTR_5)
  {
    return 1;
  }
  return 0;
}

/*
  XXX: better directly sort the ArrayLTRboundaries
*/

static const LTRboundaries **sortedltrboundaries(unsigned long *numofboundaries,
                                                 const GtArrayLTRboundaries
                                                 *ltr)
{
  unsigned long countboundaries = 0, nextfill = 0;
  const LTRboundaries *bd, **bdptrtab;

  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      countboundaries++;
    }
  }
  bdptrtab = gt_malloc(sizeof (LTRboundaries *) * countboundaries);
  nextfill = 0;
  for (bd = ltr->spaceLTRboundaries; bd < ltr->spaceLTRboundaries +
                                          ltr->nextfreeLTRboundaries; bd++)
  {
    if (!bd->skipped)
    {
      bdptrtab[nextfill++] = bd;
    }
  }
  qsort(bdptrtab,(size_t) countboundaries, sizeof (LTRboundaries *),
        bdptrcompare);
  *numofboundaries = countboundaries;
  return bdptrtab;
}

static int runltrharvest(LTRharvestoptions *lo, GtError *err)
{
  Sequentialsuffixarrayreader *ssar;
  bool had_err = false;
  unsigned long numofboundaries;
  const LTRboundaries **bdptrtab = NULL;
  GtArrayLTRboundaries arrayLTRboundaries;  /* stores all predicted */
                                            /*   LTR elements */
  const GtEncodedsequence *encseq;

  gt_error_check(err);

  ssar = gt_newSequentialsuffixarrayreaderfromfile(lo->str_indexname,
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB | SARR_DESTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    return -1;
  }
  encseq = gt_encseqSequentialsuffixarrayreader(ssar);

  /* test if motif is valid and encode motif */
  if (gt_testmotifandencodemotif (&lo->motif, encseq, err) != 0)
  {
    had_err = true;
  }

  /* show defined option and values */
  if (!had_err && lo->verbosemode)
  {
    gt_showuserdefinedoptionsandvalues(lo);
  }

  /* init array for maximal repeats */
  GT_INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.encseq = encseq;

  /* search for maximal repeats */
  if (!had_err && gt_enumeratemaxpairs(ssar,
                                   encseq,
                                   gt_readmodeSequentialsuffixarrayreader(ssar),
                                   (unsigned int) lo->minseedlength,
                                   gt_simpleexactselfmatchstore,
                                   &lo->repeatinfo,
                                   NULL,
                                   err) != 0)
  {
    had_err = true;
  }

  /* init array for candidate pairs */
  GT_INITARRAY(&arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if (!had_err && gt_searchforLTRs (lo, &arrayLTRboundaries, encseq, err) != 0)
  {
    had_err = true;
  }

  /* free array for maximal repeats */
  GT_FREEARRAY(&lo->repeatinfo.repeats, Repeat);

  /* remove exact duplicates */
  if (!had_err)
  {
    gt_removeduplicates(&arrayLTRboundaries);
  }

  /* remove overlapping predictions if desired */
  if (!had_err && (lo->nooverlapallowed || lo->bestofoverlap))
  {
    gt_removeoverlapswithlowersimilarity(&arrayLTRboundaries,
                                         lo->nooverlapallowed);
  }

  if (!had_err)
  {
    bdptrtab = sortedltrboundaries(&numofboundaries,&arrayLTRboundaries);
  }

  /* print multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutput)
  {
    if (gt_showpredictionsmultiplefasta(lo,
                                     bdptrtab,
                                     numofboundaries,
                                     false,
                                     60U,
                                     true,
                                     err) != 0)
    {
      had_err = true;
    }
  }

  /* print inner region multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutputinnerregion)
  {
    if (gt_showpredictionsmultiplefasta(lo,
                                     bdptrtab,
                                     numofboundaries,
                                     true,
                                     60U,
                                     true,
                                     err) != 0)
    {
      had_err = true;
    }
  }

  /* print GFF3 format file of predictions */
  if (!had_err && lo->gff3output && numofboundaries > 0)
  {
    if (gt_printgff3format(lo,bdptrtab,numofboundaries,encseq,err) != 0)
    {
      had_err = true;
    }
  }

  /* print predictions to stdout */
  if (!had_err)
  {
    gt_showinfoiffoundfullLTRs(lo,bdptrtab,numofboundaries,encseq);
  }

  GT_FREEARRAY(&arrayLTRboundaries, LTRboundaries);
  gt_freeSequentialsuffixarrayreader(&ssar);
  gt_free(bdptrtab);

  return had_err ? -1 : 0;
}

int gt_parseargsandcallltrharvest(int argc,const char *argv[],GtError *err)
{
  LTRharvestoptions lo;
  int had_err = 0;

  if (ltrharvestoptions(&lo,argc,argv,err) != 0)
  {
    had_err = -1;
  } else
  {
    gt_printargsline(argv,argc);
    had_err = runltrharvest(&lo,err);
  }
  gt_wrapltrharvestoptions(&lo);
  return had_err;
}
