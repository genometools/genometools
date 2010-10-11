/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <math.h>
#include "core/fa.h"
#include "fmindex.h"

#include "fmi-mapspec.pr"

static int writefmascii(const char *indexname,
                        const Fmindex *fm,
                        const GtSpecialcharinfo *specialcharinfo,
                        bool storeindexpos,
                        GtError *err)
{
  FILE *fmafp;

  gt_error_check(err);
  fmafp = gt_fa_fopen_with_suffix(indexname,FMASCIIFILESUFFIX,"wb",err);
  if (fmafp == NULL)
  {
    return -1;
  }
  fprintf (fmafp, "bwtlength=%lu\n",
           fm->bwtlength);
  fprintf (fmafp, "longest=%lu\n",
                   fm->longestsuffixpos);
  fprintf (fmafp, "storeindexpos=%d\n", storeindexpos ? 1 : 0);
  fprintf (fmafp, "log2blocksize=%u\n", fm->log2bsize);
  fprintf (fmafp, "log2markdist=%u\n", fm->log2markdist);
  fprintf (fmafp, "specialcharacters=%lu\n",specialcharinfo->specialcharacters);
  fprintf (fmafp, "specialranges=%lu\n",specialcharinfo->specialranges);
  fprintf (fmafp, "realspecialranges=%lu\n",specialcharinfo->realspecialranges);
  fprintf (fmafp, "lengthofspecialprefix=%lu\n",
           specialcharinfo->lengthofspecialprefix);
  fprintf (fmafp, "lengthofspecialsuffix=%lu\n",
           specialcharinfo->lengthofspecialsuffix);
  fprintf (fmafp, "wildcards=%lu\n",specialcharinfo->wildcards);
  fprintf (fmafp, "wildcardranges=%lu\n",specialcharinfo->wildcardranges);
  fprintf (fmafp, "realwildcardranges=%lu\n",
                   specialcharinfo->realwildcardranges);
  fprintf (fmafp, "lengthofwildcardprefix=%lu\n",
           specialcharinfo->lengthofwildcardprefix);
  fprintf (fmafp, "lengthofwildcardsuffix=%lu\n",
           specialcharinfo->lengthofwildcardsuffix);
  fprintf (fmafp, "suffixlength=%u\n", fm->suffixlength);
  gt_fa_xfclose(fmafp);
  return 0;
}

static int writefmdata(const char *indexname,
                       Fmindex *fm,
                       bool storeindexpos,
                       GtError *err)
{
  FILE *fp;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,FMDATAFILESUFFIX,"wb",err);
  if (fp == NULL)
  {
    return -1;
  }
  if (gt_flushfmindex2file(fp,fm,storeindexpos,err) != 0)
  {
    return -2;
  }
  gt_fa_xfclose(fp);
  return 0;
}

int gt_saveFmindex(const char *indexname,Fmindex *fm,
                   const GtSpecialcharinfo *specialcharinfo,
                   bool storeindexpos,GtError *err)
{
  gt_error_check(err);
  if (writefmascii(indexname, fm, specialcharinfo,storeindexpos,err) != 0)
  {
    return -1;
  }
  if (writefmdata (indexname, fm, storeindexpos,err) != 0)
  {
    return -2;
  }
  return 0;
}
