#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/symboldef.h"
#include "esafileend.h"
#include "sfx-outlcp.h"

#include "opensfxfile.pr"

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)\
        if (!haserr && (FLAG))\
        {\
          PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",env);\
          if ((PTR) == NULL)\
          {\
            haserr = true;\
          }\

int initlcpoutfileinfo(Outlcpinfo *outlcpinfo,const Str *indexname,Env *env,
                       bool origin)
{
  bool haserr = false;

  if (indexname == NULL)
  {
    outlcpinfo->outfplcptab = NULL;
    outlcpinfo->outfpllvtab = NULL;
  } else
  {
    outlcpinfo->outfplcptab = opensfxfile(indexname,
                                          origin ? LCPTABSUFFIX
                                                 : LCPTABSUFFIX "2",
                                          "wb",env);
    if (outlcpinfo->outfplcptab == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      outlcpinfo->outfpllvtab 
        = opensfxfile(indexname,origin ? LARGELCPTABSUFFIX
                                       : LARGELCPTABSUFFIX "2",
                                       "wb",env);
      if (outlcpinfo->outfpllvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  outlcpinfo->numoflargelcpvalues = 0;
  outlcpinfo->maxbranchdepth = 0;
  return haserr ? -1 : 0;
}

int outlcpvalue(Seqpos lcpvalue,Seqpos pos,Seqpos pageoffset,
                Outlcpinfo *outlcpinfo,Env *env)
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
      env_error_set(env,"cannot write 1 item of size %lu: "
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
    env_error_set(env,"cannot write 1 item of size %lu: "
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

void freeoutlcptab(Outlcpinfo *outlcpinfo,Env *env)
{
  fa_fclose(outlcpinfo->outfplcptab);
  fa_fclose(outlcpinfo->outfpllvtab);
}
