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

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "libgtcore/fa.h"
#include "libgtcore/array.h"
#include "libgtcore/str.h"
#include "sfx-ri-def.h"
#include "esafileend.h"
#include "fmindex.h"
#include "sarr-def.h"
#include "verbose-def.h"
#include "spacedef.h"
#include "stamp.h"

#include "opensfxfile.pr"
#include "esa-map.pr"
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
                                 Verboseinfo *verboseinfo,
                                 Env *env)
{
  bool haserr = false;
  Array *riktab;
  unsigned int intstoreindexpos;

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
    Str *currentline;
    unsigned int linenum;

    currentline = str_new();
    for (linenum = 0; str_read_next_line(currentline, fpin, env) != EOF;
         linenum++)
    {
      if (analyzeuintline(indexname,
                         FMASCIIFILESUFFIX,
                         linenum,
                         str_get(currentline),
                         str_length(currentline),
                         riktab,
                         env) != 0)
      {
        haserr = true;
        break;
      }
      str_reset(currentline);
    }
    str_delete(currentline);
  }
  if (!haserr && allkeysdefined(indexname,FMASCIIFILESUFFIX,riktab,
                                verboseinfo,env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (intstoreindexpos == (unsigned int) 1)
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
  array_delete(riktab,env);
  return haserr ? -1 : 0;
}

void freefmindex(Fmindex *fmindex,Env *env)
{
  if (fmindex->mappedptr != NULL)
  {
    fa_xmunmap(fmindex->mappedptr);
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

static Encodedsequence *mapbwtencoding(const Str *indexname,
                                       Verboseinfo *verboseinfo,
                                       Env *env)
{
  Suffixarray suffixarray;
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  if (mapsuffixarray(&suffixarray,&totallength,SARR_ESQTAB,indexname,
                     verboseinfo,env) != 0)
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

int mapfmindex (Fmindex *fmindex,const Str *indexname,
                Verboseinfo *verboseinfo,Env *env)
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
                             verboseinfo,
                             env) != 0)
    {
      haserr = true;
    }
  }
  fa_xfclose(fpin);
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
    str_delete(tmpfilename);
  }
  if (!haserr)
  {
    fmindex->bwtformatching = mapbwtencoding(indexname,verboseinfo,env);
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
    str_delete(tmpfilename);
  }
  if (haserr)
  {
    freefmindex(fmindex,env);
  }
  return haserr ? -1 : 0;
}
