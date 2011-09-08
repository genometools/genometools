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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include "core/endianess_api.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/str_array.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/defined-types.h"
#include "core/format64.h"
#include "core/readmode.h"
#include "core/encseq.h"

#include "esa-fileend.h"
#include "stamp.h"

#define PRJSPECIALOUT(VAL)\
        fprintf(outprj,"%s=%lu\n",#VAL,gt_encseq_##VAL(encseq))

static void showprjinfo(FILE *outprj,
                        GtReadmode readmode,
                        const GtEncseq *encseq,
                        unsigned long numberofallsortedsuffixes,
                        unsigned int prefixlength,
                        unsigned long numoflargelcpvalues,
                        double averagelcp,
                        unsigned long maxbranchdepth,
                        const Definedunsignedlong *longest)
{
  unsigned long totallength;
  unsigned long numofsequences;

  totallength = gt_encseq_total_length(encseq);
  fprintf(outprj,"totallength=%lu\n",totallength);
  PRJSPECIALOUT(specialcharacters);
  PRJSPECIALOUT(specialranges);
  PRJSPECIALOUT(realspecialranges);
  PRJSPECIALOUT(lengthofspecialprefix);
  PRJSPECIALOUT(lengthofspecialsuffix);
  PRJSPECIALOUT(wildcards);
  PRJSPECIALOUT(wildcardranges);
  PRJSPECIALOUT(realwildcardranges);
  PRJSPECIALOUT(lengthofwildcardprefix);
  PRJSPECIALOUT(lengthofwildcardsuffix);
  numofsequences = gt_encseq_num_of_sequences(encseq);
  fprintf(outprj,"numofsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofdbsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofquerysequences=0\n");
  fprintf(outprj,"numberofallsortedsuffixes=%lu\n",numberofallsortedsuffixes);
  if (longest->defined)
  {
    fprintf(outprj,"longest=%lu\n",longest->valueunsignedlong);
  }
  fprintf(outprj,"prefixlength=%u\n",prefixlength);
  fprintf(outprj,"largelcpvalues=%lu\n",numoflargelcpvalues);
  fprintf(outprj,"averagelcp=%.2f\n",averagelcp);
  fprintf(outprj,"maxbranchdepth=%lu\n",maxbranchdepth);
  fprintf(outprj,"integersize=%u\n",
                  (unsigned int) (sizeof (unsigned long) * CHAR_BIT));
  fprintf(outprj,"littleendian=%c\n",gt_is_little_endian() ? '1' : '0');
  fprintf(outprj,"readmode=%u\n",(unsigned int) readmode);
  fprintf(outprj,"mirrored=%c\n", gt_encseq_is_mirrored(encseq) ? '1' : '0');
}

int gt_outprjfile(const char *indexname,
                  GtReadmode readmode,
                  const GtEncseq *encseq,
                  unsigned long numberofallsortedsuffixes,
                  unsigned int prefixlength,
                  unsigned long numoflargelcpvalues,
                  double averagelcp,
                  unsigned long maxbranchdepth,
                  const Definedunsignedlong *longest,
                  GtError *err)
{
  FILE *prjfp;
  bool haserr = false;

  gt_error_check(err);
  prjfp = gt_fa_fopen_with_suffix(indexname,PROJECTFILESUFFIX,"wb",err);
  if (prjfp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    showprjinfo(prjfp,
                readmode,
                encseq,
                numberofallsortedsuffixes,
                prefixlength,
                numoflargelcpvalues,
                averagelcp,
                maxbranchdepth,
                longest);
    gt_fa_xfclose(prjfp);
  }
  return haserr ? -1 : 0;
}
