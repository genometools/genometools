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

#include <errno.h>
#include <limits.h>
#include <string.h>
#include "sarr-def.h"
#include "spacedef.h"
#include "emimergeesa.h"
#include "esafileend.h"
#include "verbose-def.h"

#include "esa-merge.pr"
#include "encseq2offset.pr"

typedef struct
{
  Str *outfilename;
  FILE *fp;
} NameandFILE;

typedef struct
{
  NameandFILE outsuf,
              outlcp,
              outllv;
  Seqpos currentlcpindex,
         absstartpostable[SIZEOFMERGERESULTBUFFER];
} Mergeoutinfo;

static int initNameandFILE(NameandFILE *nf,
                            const Str *outindex,
                            const char *suffix,
                            Env *env)
{
  env_error_check(env);
  nf->outfilename = str_clone(outindex,env);
  str_append_cstr(nf->outfilename,suffix,env);
  nf->fp = env_fa_fopen(env,str_get(nf->outfilename),"wb");
  if (nf->fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",str_get(nf->outfilename),
                                                    strerror(errno));
    return -1;
  }
  return 0;
}

static void freeNameandFILE(NameandFILE *nf,Env *env)
{
  env_fa_xfclose(nf->fp,env);
  str_delete(nf->outfilename,env);
}

static int outputsuflcpllv(void *processinfo,
                           const Seqpos *sequenceoffsettable,
                           const Suflcpbuffer *buf,
                           Env *env)
{
  Mergeoutinfo *mergeoutinfo = (Mergeoutinfo *) processinfo;

  unsigned int i, lastindex;
  Seqpos lcpvalue;
  Largelcpvalue currentexception;
  Uchar smallvalue;
  bool haserr = false;

  env_error_check(env);
  for (i=0; i<buf->nextstoreidx; i++)
  {
    mergeoutinfo->absstartpostable[i]
      = sequenceoffsettable[buf->suftabstore[i].idx] +
        buf->suftabstore[i].startpos;
  }
  if (fwrite(mergeoutinfo->absstartpostable,
            sizeof (Seqpos),
            (size_t) buf->nextstoreidx,
            mergeoutinfo->outsuf.fp)
         != (size_t) buf->nextstoreidx)
  {
    env_error_set(env,"fwrite(%s) of %u Seqpos-value failed: %s",
                  str_get(mergeoutinfo->outsuf.outfilename),
                  buf->nextstoreidx,strerror(errno));
    haserr = true;
  }
  if (!haserr)
  {
    if (buf->lastpage)
    {
      lastindex = buf->nextstoreidx - 1;
    } else
    {
      lastindex = buf->nextstoreidx;
    }
    for (i=0; i<lastindex; i++)
    {
      lcpvalue = buf->lcptabstore[i];
      if (lcpvalue < (Seqpos) UCHAR_MAX)
      {
        smallvalue = (Uchar) lcpvalue;
      } else
      {
        currentexception.position = mergeoutinfo->currentlcpindex;
        currentexception.value = lcpvalue;
        if (fwrite(&currentexception,sizeof (Largelcpvalue),
                 (size_t) 1,mergeoutinfo->outllv.fp) != (size_t) 1)
        {
          env_error_set(env,"fwrite(%s) of Largelcpvalue failed: %s",
                        str_get(mergeoutinfo->outllv.outfilename),
                        strerror(errno));
          haserr = true;
          break;
        }
        smallvalue = (Uchar) UCHAR_MAX;
      }
      if (fwrite(&smallvalue,sizeof (Uchar),(size_t) 1,
                mergeoutinfo->outlcp.fp) != (size_t) 1)
      {
        env_error_set(env,"fwrite(%s) of Uchar failed: %s",
                       str_get(mergeoutinfo->outlcp.outfilename),
                       strerror(errno));
        haserr = true;
        break;
      }
      mergeoutinfo->currentlcpindex++;
    }
  }
  return haserr ? -1 : 0;
}

static int mergeandstoreindex(const Str *storeindex,
                              Emissionmergedesa *emmesa,
                              Env *env)
{
  Mergeoutinfo mergeoutinfo;
  Uchar smalllcpvalue;
  Specialcharinfo specialcharinfo;
  Seqpos *sequenceoffsettable, totallength;
  bool haserr = false;

  env_error_check(env);
  if (initNameandFILE(&mergeoutinfo.outsuf,storeindex,SUFTABSUFFIX,env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (initNameandFILE(&mergeoutinfo.outlcp,storeindex,LCPTABSUFFIX,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (initNameandFILE(&mergeoutinfo.outllv,storeindex,LARGELCPTABSUFFIX,
                       env) != 0)
    {
      haserr = true;
    }
  }
  smalllcpvalue = 0;
  if (!haserr && fwrite(&smalllcpvalue,sizeof (Uchar),(size_t) 1,
                mergeoutinfo.outlcp.fp) != (size_t) 1)
  {
    env_error_set(env,"fwrite(%s) failed: %s",
                  str_get(mergeoutinfo.outlcp.outfilename),
                  strerror(errno));
    haserr = true;
  }
  if (!haserr)
  {
    mergeoutinfo.currentlcpindex = (Seqpos) 1;
    sequenceoffsettable = encseqtable2seqoffsets(&totallength,
                                                 &specialcharinfo,
                                                 emmesa->suffixarraytable,
                                                 emmesa->numofindexes,
                                                 env);
    assert(sequenceoffsettable != NULL);
    while (emmesa->numofentries > 0)
    {
      if (stepdeleteandinsertothersuffixes(emmesa,env) != 0)
      {
        haserr = true;
        break;
      }
      if (outputsuflcpllv(&mergeoutinfo,
                         sequenceoffsettable,
                         &emmesa->buf,
                         env) != 0)
      {
        haserr = true;
        break;
      }
    }
    FREESPACE(sequenceoffsettable);
  }
  freeNameandFILE(&mergeoutinfo.outsuf,env);
  freeNameandFILE(&mergeoutinfo.outlcp,env);
  freeNameandFILE(&mergeoutinfo.outllv,env);
  return haserr ? -1 : 0;
}

int performtheindexmerging(const Str *storeindex,
                           const StrArray *indexnametab,
                           Verboseinfo *verboseinfo,
                           Env *env)
{
  Emissionmergedesa emmesa;
  unsigned int demand = SARR_ESQTAB | SARR_SUFTAB | SARR_LCPTAB;
  bool haserr = false;

  env_error_check(env);
  if (initEmissionmergedesa(&emmesa,
                            indexnametab,
                            demand,
                            verboseinfo,
                            env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (strarray_size(indexnametab) > (unsigned long) 1)
    {
      if (mergeandstoreindex(storeindex,&emmesa,env) != 0)
      {
        haserr = true;
      }
    } else
    {
      env_error_set(env,"merging requires more than one index");
      haserr = true;
    }
  }
  wraptEmissionmergedesa(&emmesa,env);
  return haserr ? -1 : 0;
}
