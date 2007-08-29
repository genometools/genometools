/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "libgtcore/array.h"
#include "libgtcore/str.h"
#include "arraydef.h"
#include "sfx-ri-def.h"
#include "esafileend.h"
#include "fmindex.h"
#include "sarr-def.h"
#include "stamp.h"

#include "readnextline.pr"
#include "opensfxfile.pr"
#include "sfx-readint.pr"
#include "sfx-map.pr"
#include "fmi-keyval.pr"
#include "fmi-mapspec.pr"

bool fmindexexists(const Str *indexname,Env *env)
{
  env_error_check(env);
  if (!indexfilealreadyexists(indexname,FMASCIIFILESUFFIX,env))
  {
    return false;
  }
  if (!indexfilealreadyexists(indexname,FMDATAFILESUFFIX,env))
  {
    return false;
  }
  return true;
}

static int scanfmafileviafileptr(Fmindex *fmindex,
                                 bool *storeindexpos,
                                 const Str *indexname,
                                 FILE *fpin,
                                 Env *env)
{
  ArrayUchar linebuffer;
  bool haserr = false;
  Array *riktab;
  uint32_t linenum,
           intstoreindexpos;

  env_error_check(env);
  riktab = array_new(sizeofReadintkeys(),env);
  SETREADINTKEYS("bwtlength",&fmindex->bwtlength,NULL);
  SETREADINTKEYS("longest",&fmindex->longestsuffixpos,NULL);
  SETREADINTKEYS("storeindexpos",&intstoreindexpos,NULL);
  SETREADINTKEYS("log2blocksize",&fmindex->log2bsize,NULL);
  SETREADINTKEYS("log2markdist",&fmindex->log2markdist,NULL);
  SETREADINTKEYS("specialcharacters",
                 &fmindex->specialcharinfo.specialcharacters,NULL);
  SETREADINTKEYS("specialranges",
                 &fmindex->specialcharinfo.specialranges,NULL);
  SETREADINTKEYS("lengthofspecialprefix",
                 &fmindex->specialcharinfo.lengthofspecialprefix,NULL);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &fmindex->specialcharinfo.lengthofspecialsuffix,NULL);
  SETREADINTKEYS("suffixlength",&fmindex->suffixlength,NULL);
  if (!haserr)
  {
    INITARRAY(&linebuffer,Uchar);
    for (linenum = 0; /* Nothing */; linenum++)
    {
      linebuffer.nextfreeUchar = 0;
      if (readnextline(fpin,&linebuffer,env) == EOF)
      {
        break;
      }
      if (analyzeuintline(indexname,
                         FMASCIIFILESUFFIX,
                         linenum,
                         linebuffer.spaceUchar,
                         linebuffer.nextfreeUchar,
                         riktab,
                         env) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (!haserr && allkeysdefined(indexname,FMASCIIFILESUFFIX,riktab,false,
                                env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (intstoreindexpos == (uint32_t) 1)
    {
      *storeindexpos = true;
    } else
    {
      if (intstoreindexpos == 0)
      {
        *storeindexpos = false;
      } else
      {
        env_error_set(env,"illegal value in line matching \"storeindexpos=\"");
        haserr = true;
      }
    }
  }
  FREEARRAY(&linebuffer,Uchar);
  array_delete(riktab,env);
  return haserr ? -1 : 0;
}

void freefmindex(Fmindex *fmindex,Env *env)
{
  if (fmindex->mappedptr != NULL)
  {
    env_fa_xmunmap(fmindex->mappedptr,env);
  }
  if (fmindex->bwtformatching != NULL)
  {
    freeEncodedsequence(&fmindex->bwtformatching,env);
  }
  if (fmindex->alphabet != NULL)
  {
    freeAlphabet(&fmindex->alphabet,env);
  }
}

static Encodedsequence *mapbwtencoding(const Str *indexname,Env *env)
{
  Suffixarray suffixarray;
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  if (mapsuffixarray(&suffixarray,&totallength,SARR_ESQTAB,indexname,
                    false,env) != 0)
  {
    haserr = true;
  }
  freeAlphabet(&suffixarray.alpha,env);
  strarray_delete(suffixarray.filenametab,env);
  FREESPACE(suffixarray.filelengthtab);
  if (haserr)
  {
    freeEncodedsequence(&suffixarray.encseq,env);
    return NULL;
  }
  return suffixarray.encseq;
}

int mapfmindex (Fmindex *fmindex,const Str *indexname,Env *env)
{
  FILE *fpin = NULL;
  bool haserr = false, storeindexpos = true;

  env_error_check(env);
  fmindex->mappedptr = NULL;
  fmindex->bwtformatching = NULL;
  fpin = opensfxfile(indexname,FMASCIIFILESUFFIX,"rb",env);
  if (fpin == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (scanfmafileviafileptr(fmindex,
                             &storeindexpos,
                             indexname,
                             fpin,
                             env) != 0)
    {
      haserr = true;
    }
  }
  env_fa_xfclose(fpin, env);
  if (!haserr)
  {
    Str *tmpfilename;

    fmindex->specpos.nextfreePairBwtidx
      = (unsigned long) determinenumberofspecialstostore(
                                          &fmindex->specialcharinfo);
    fmindex->specpos.spacePairBwtidx = NULL;
    fmindex->specpos.allocatedPairBwtidx = 0;
    tmpfilename = str_clone(indexname,env);
    str_append_cstr(tmpfilename,ALPHABETFILESUFFIX,env);
    fmindex->alphabet = assigninputalphabet(false,
                                            false,
                                            tmpfilename,
                                            NULL,
                                            env);
    if (fmindex->alphabet == NULL)
    {
      haserr = true;
    }
    str_delete(tmpfilename,env);
  }
  if (!haserr)
  {
    fmindex->bwtformatching = mapbwtencoding(indexname,env);
    if (fmindex->bwtformatching == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    Str *tmpfilename;

    computefmkeyvalues (fmindex,
                        fmindex->bwtlength,
                        fmindex->log2bsize,
                        fmindex->log2markdist,
                        getmapsizeAlphabet(fmindex->alphabet),
                        fmindex->suffixlength,
                        storeindexpos,
                        &fmindex->specialcharinfo);
    tmpfilename = str_clone(indexname,env);
    str_append_cstr(tmpfilename,FMDATAFILESUFFIX,env);
    if (fillfmmapspecstartptr(fmindex,storeindexpos,tmpfilename,env) != 0)
    {
      haserr = true;
    }
    str_delete(tmpfilename,env);
  }
  if (haserr)
  {
    freefmindex(fmindex,env);
  }
  return haserr ? -1 : 0;
}
