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

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/symboldef.h"
#include "spacedef.h"
#include "esafileend.h"
#include "sfx-outlcp.h"

#include "opensfxfile.pr"

 struct Outlcpinfo
{
  FILE *outfplcptab,
       *outfpllvtab;
  Seqpos numoflargelcpvalues,
         maxbranchdepth;
};

Outlcpinfo *newlcpoutfileinfo(const Str *indexname,Error *err,bool origin)
{
  bool haserr = false;
  Outlcpinfo *outlcpinfo;

  ALLOCASSIGNSPACE(outlcpinfo,NULL,Outlcpinfo,1);
  if (indexname == NULL)
  {
    outlcpinfo->outfplcptab = NULL;
    outlcpinfo->outfpllvtab = NULL;
  } else
  {
    outlcpinfo->outfplcptab = opensfxfile(indexname,
                                          origin ? LCPTABSUFFIX
                                                 : LCPTABSUFFIX "2",
                                          "wb",err);
    if (outlcpinfo->outfplcptab == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      outlcpinfo->outfpllvtab
        = opensfxfile(indexname,origin ? LARGELCPTABSUFFIX
                                       : LARGELCPTABSUFFIX "2",
                                       "wb",err);
      if (outlcpinfo->outfpllvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  outlcpinfo->numoflargelcpvalues = 0;
  outlcpinfo->maxbranchdepth = 0;
  if (haserr)
  {
    FREESPACE(outlcpinfo);
    return NULL;
  }
  return outlcpinfo;
}

int outlcpvalue(Seqpos lcpvalue,Seqpos pos,Seqpos pageoffset,
                Outlcpinfo *outlcpinfo,Error *err)
{
  Uchar outvalue;
  bool haserr = false;

  if (lcpvalue >= (Seqpos) UCHAR_MAX)
  {
    Largelcpvalue largelcpvalue;

    outlcpinfo->numoflargelcpvalues++;
    largelcpvalue.position = pageoffset + pos;
    largelcpvalue.value = lcpvalue;
    if (fwrite(&largelcpvalue,sizeof (Largelcpvalue),(size_t) 1,
               outlcpinfo->outfpllvtab) != (size_t) 1)
    {
      error_set(err,"cannot write 1 item of size %lu: "
                        "errormsg=\"%s\"",
                        (unsigned long) sizeof (Largelcpvalue),
                        strerror(errno));
      haserr = true;
    }
    outvalue = (Uchar) UCHAR_MAX;
  } else
  {
    outvalue = (Uchar) lcpvalue;
  }
  if (!haserr && fwrite(&outvalue,sizeof (Uchar),(size_t) 1,
                        outlcpinfo->outfplcptab) != (size_t) 1)
  {
    error_set(err,"cannot write 1 item of size %lu: "
                      "errormsg=\"%s\"",
                      (unsigned long) sizeof (Uchar),
                      strerror(errno));
    haserr = true;
  }
  if (outlcpinfo->maxbranchdepth < lcpvalue)
  {
    outlcpinfo->maxbranchdepth = lcpvalue;
  }
  return haserr ? -1 : 0;
}

void freeoutlcptab(Outlcpinfo **outlcpinfo,Error *err)
{
  fa_fclose((*outlcpinfo)->outfplcptab);
  fa_fclose((*outlcpinfo)->outfpllvtab);
  FREESPACE(*outlcpinfo);
}

Seqpos getnumoflargelcpvalues(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->numoflargelcpvalues;
}

Seqpos getmaxbranchdepth(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->maxbranchdepth;
}
