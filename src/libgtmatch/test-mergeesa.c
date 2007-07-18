/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <errno.h>
#include "sarr-def.h"
#include "emimergeesa.h"
#include "esafileend.h"

#include "mergeesa.pr"
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

static void initNameandFILE(NameandFILE *nf,
                            const Str *outindex,
                            const char *suffix,
                            Env *env)
{
  nf->outfilename = str_clone(outindex,env);
  str_append_cstr(nf->outfilename,suffix,env);
  nf->fp = env_fa_fopen(env,str_get(nf->outfilename),"wb");
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
  
  uint32_t i, lastindex;
  Seqpos lcpvalue;
  Largelcpvalue currentexception;
  Uchar smallvalue;
  bool haserr = false;

  for(i=0; i<buf->nextstoreidx; i++)
  {
    mergeoutinfo->absstartpostable[i] 
      = sequenceoffsettable[buf->suftabstore[i].idx] +
        buf->suftabstore[i].startpos;
  }
  if(fwrite(mergeoutinfo->absstartpostable,
            sizeof(Seqpos),
            (size_t) buf->nextstoreidx,
            mergeoutinfo->outsuf.fp) 
         != (size_t) buf->nextstoreidx)
  {
    env_error_set(env,"fwrite(%s) of %u Seqpos-value failed: %s",
                  str_get(mergeoutinfo->outsuf.outfilename),
                  (unsigned int) buf->nextstoreidx,strerror(errno));
    haserr = true;
  }
  if(!haserr)
  {
    if(buf->lastpage)
    {
      lastindex = buf->nextstoreidx - 1;
    } else
    {
      lastindex = buf->nextstoreidx;
    }
    for(i=0; i<lastindex; i++)
    {
      lcpvalue = buf->lcptabstore[i];
      if(lcpvalue < UCHAR_MAX)
      {
        smallvalue = (Uchar) lcpvalue;
      } else
      {
        currentexception.position = mergeoutinfo->currentlcpindex;
        currentexception.value = lcpvalue;
        if(fwrite(&currentexception,sizeof(Largelcpvalue),
                 (size_t) 1,mergeoutinfo->outllv.fp) != (size_t) 1)
        {
          env_error_set(env,"fwrite(%s) of Largelcpvalue failed: %s",
                        str_get(mergeoutinfo->outllv.outfilename),
                        strerror(errno));
          haserr = true;
          break;
        }
        smallvalue = UCHAR_MAX;
      }
      if(fwrite(&smallvalue,sizeof(Uchar),(size_t) 1,
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

  initNameandFILE(&mergeoutinfo.outsuf,storeindex,SUFTABSUFFIX,env);
  initNameandFILE(&mergeoutinfo.outlcp,storeindex,LCPTABSUFFIX,env);
  initNameandFILE(&mergeoutinfo.outllv,storeindex,LARGELCPTABSUFFIX,env);
  smalllcpvalue = 0;
  if(fwrite(&smalllcpvalue,sizeof(Uchar),(size_t) 1,
            mergeoutinfo.outlcp.fp) != (size_t) 1)
  {
    env_error_set(env,"fwrite(%s) failed: %s",
                  str_get(mergeoutinfo.outlcp.outfilename),
                  strerror(errno));
    haserr = true;
  }
  if(!haserr)
  {
    mergeoutinfo.currentlcpindex = (Seqpos) 1;
    sequenceoffsettable = encseqtable2seqoffsets(&totallength,
                                                 &specialcharinfo,
                                                 emmesa->suffixarraytable,
                                                 emmesa->numofindexes,
                                                 env);
    assert(sequenceoffsettable != NULL);
    while(emmesa->numofentries > 0)
    {
      if(stepdeleteandinsertothersuffixes(emmesa,env) != 0)
      {
        haserr = true;
        break;
      }
      if(outputsuflcpllv(&mergeoutinfo,
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

/*

static int dothemerging(Mergeesacallinfo *mergeesacallinfo)
{
  Emissionmergedesa emmesa;
  Uint demand;

  if(mergeesacallinfo->numofindexes == UintConst(1) &&
     mergeesacallinfo->userdefinedminlength.defined)
  {
    demand = BWTTABSTREAM | SUFTABSTREAM | LCPTABSTREAM;
  } else
  {
    demand = TISTABSTREAM | SUFTABSTREAM | LCPTABSTREAM;
  }
  if(initEmissionmergedesa(&emmesa,
                           (const char **) mergeesacallinfo->indexnamelist,
                           mergeesacallinfo->numofindexes,
                           demand) != 0)
  {
    return -1;
  }
  if(mergeesacallinfo->storeindex != NULL)
  {
    if(mergeesacallinfo->numofindexes > UintConst(1))
    {
      if(mergeandstoreindex(mergeesacallinfo->is32bit,
                            mergeesacallinfo->storeindex,&emmesa) != 0)
      {
        return -2;
      }
    } else
    {
      ERROR0("merging requires more than one index");
      return -3;
    }
  } else
  {
    if(mergeesacallinfo->maxuniquelength > 0)
    {
      if(mergeesacallinfo->numofindexes > UintConst(1))
      {
        if(multienumuniquedist(&emmesa,
                               mergeesacallinfo->minuniquelength,
                               mergeesacallinfo->maxuniquelength) != 0)
        {
          return -4;
        }
      } else
      {
        if(singleenumuniquedist(&emmesa.encseqtable[0],
                                &emmesa.esastreamtable[0],
                                mergeesacallinfo->minuniquelength,
                                mergeesacallinfo->maxuniquelength) != 0)
        {
          return -5;
        }
      }
    } else
    {
      if(mergeesacallinfo->userdefinedminlength.defined)
      {
        if(mergeesacallinfo->numofindexes > UintConst(1))
        {
          Specialcharinfo specialcharinfo;
          Uint64 *sequenceoffsettable, 
                 totallength; 

          sequenceoffsettable 
            = encseqtable2seqoffsets(mergeesacallinfo->is32bit,
                                     &totallength,
                                     &specialcharinfo,
                                     emmesa.encseqtable,
                                     emmesa.numofindexes);
          if(sequenceoffsettable == NULL)
          {
            return -2;
          }
          if(mergevmatmaxoutdynamic(
                      &emmesa,
                      UintConst(1),
                      mergeesacallinfo->userdefinedminlength.uintvalue,
                      NULL,
                      (void *) sequenceoffsettable,
                      (void *) mergesimpleexactselfmatchoutput) != 0)
          {
            return -6;
          }
          FREESPACE(sequenceoffsettable);
        } else
        {
          if(strmvmatmaxoutdynamic(
                  &emmesa.esastreamtable[0],
                  UintConst(1),
                  mergeesacallinfo->userdefinedminlength.uintvalue,
                  NULL,
                  NULL,
                  (void *) simpleexactselfmatchoutput) != 0)
          {
            return -7;
          }
        }
      }
      {
         printf("# construct the fmindex \"%s\"\n",
                    mergeesacallinfo->storefmindex);
      }
    }
  }
  wraptEmissionmergedesa(&emmesa);
  return 0;
}
*/
