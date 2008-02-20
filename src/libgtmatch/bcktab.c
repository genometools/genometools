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

#include <string.h>
#include <errno.h>
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/fa.h"
#include "esafileend.h"
#include "spacedef.h"
#include "bckbound.h"

#include "initbasepower.pr"

/* DELETE THE FOLLOWING THREE FUNCTIONS */

static void *genericmaponlytable(const Str *indexname,const char *suffix,
                                 size_t *numofbytes,Error *err)
{
  Str *tmpfilename;
  void *ptr;
  bool haserr = false;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,suffix);
  ptr = fa_mmap_read(str_get(tmpfilename),numofbytes);
  if (ptr == NULL)
  {
    error_set(err,"cannot map file \"%s\": %s",str_get(tmpfilename),
                  strerror(errno));
    haserr = true;
  }
  str_delete(tmpfilename);
  return haserr ? NULL : ptr;
}

static int checkmappedfilesize(size_t numofbytes,Seqpos expectedunits,
                               size_t sizeofunit,Error *err)
{
  error_check(err);
  if (expectedunits != (Seqpos) (numofbytes/sizeofunit))
  {
    error_set(err,"number of mapped units = %lu != " FormatSeqpos
                      " = expected number of integers",
                      (unsigned long) (numofbytes/sizeofunit),
                      PRINTSeqposcast(expectedunits));
    return -1;
  }
  return 0;
}

static void *genericmaptable(const Str *indexname,
                             const char *suffix,
                             Seqpos expectedunits,size_t sizeofunit,
                             Error *err)
{
  size_t numofbytes;

  void *ptr = genericmaponlytable(indexname,suffix,&numofbytes,err);
  if (ptr == NULL)
  {
    return NULL;
  }
  if (checkmappedfilesize(numofbytes,expectedunits,sizeofunit,err) != 0)
  {
    fa_xmunmap(ptr);
    return NULL;
  }
  return ptr;
}

int mapbcktab(Bcktab *bcktab,
              const Str *indexname,
              unsigned int numofchars,
              unsigned int prefixlength,
              Error *err)
{
  Codetype numofspecialcodes;

  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->filltable = initfilltable(bcktab->basepower,prefixlength);
  bcktab->multimappower = initmultimappower(numofchars,prefixlength);
  bcktab->leftborder
    = genericmaptable(indexname,
                      BCKTABSUFFIX,
                      (Seqpos) (bcktab->numofallcodes + 1 + numofspecialcodes),
                      sizeof (Seqpos),
                      err);
  if (bcktab->leftborder == NULL)
  {
    bcktab->countspecialcodes = NULL;
    return -1;
  }
  bcktab->countspecialcodes = bcktab->leftborder + bcktab->numofallcodes + 1;
  return 0;
}

void freebcktab(Bcktab *bcktab)
{
  fa_xmunmap((void *) bcktab->leftborder);
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  multimappowerfree(&bcktab->multimappower);
  FREESPACE(bcktab->filltable);
  FREESPACE(bcktab->basepower);
}

void initbcktabwithNULL(Bcktab *bcktab)
{
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->multimappower = NULL;
  bcktab->filltable = NULL;
  bcktab->basepower = NULL;
}

static unsigned long **initdistprefixindexcounts(const Codetype *basepower,
                                                 unsigned int prefixlength)
{
  if (prefixlength > 2U)
  {
    unsigned int idx;
    unsigned long *counters, numofcounters, **distpfxidx;

    for (numofcounters = 0, idx=1U; idx <= prefixlength-2; idx++)
    {
      numofcounters += basepower[idx];
    }
    assert(numofcounters > 0);
    ALLOCASSIGNSPACE(distpfxidx,NULL,unsigned long *,prefixlength-2);
    ALLOCASSIGNSPACE(counters,NULL,unsigned long,numofcounters);
    memset(counters,0,(size_t) sizeof (*counters) * numofcounters);
    distpfxidx[0] = counters;
    for (idx=1U; idx<prefixlength-2; idx++)
    {
      distpfxidx[idx] = distpfxidx[idx-1] + basepower[idx];
    }
    return distpfxidx;
  } else
  {
    return NULL;
  }
}

int allocBcktab(Bcktab *bcktab,
                unsigned int numofchars,
                unsigned int prefixlength,
                unsigned int codebits,
                unsigned int maxcodevalue,
                Error *err)
{
  bool haserr = false;

  bcktab->distpfxidx = NULL;
  bcktab->filltable = NULL;
  bcktab->basepower = NULL;
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->filltable = initfilltable(bcktab->basepower,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  if (bcktab->numofallcodes-1 > maxcodevalue)
  {
    error_set(err,"alphasize^prefixlength-1 = %u does not fit into "
                  " %u bits: choose smaller value for prefixlength",
                  bcktab->numofallcodes-1,
                  codebits);
    haserr = true;
  } else
  {
    bcktab->distpfxidx = initdistprefixindexcounts(bcktab->basepower,
                                                   prefixlength);
  }
  if (!haserr)
  {
    ALLOCASSIGNSPACE(bcktab->leftborder,NULL,Seqpos,
                     bcktab->numofallcodes+1);
    memset(bcktab->leftborder,0,
           sizeof (*bcktab->leftborder) *
           (size_t) bcktab->numofallcodes);
    ALLOCASSIGNSPACE(bcktab->countspecialcodes,NULL,Seqpos,
                     bcktab->numofspecialcodes);
    memset(bcktab->countspecialcodes,0,
           sizeof (*bcktab->countspecialcodes) *
                  (size_t) bcktab->numofspecialcodes);
  }
  return haserr ? -1 : 0;
}

void updatebckspecials(Bcktab *bcktab,
                       Codetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex,
                       unsigned int prefixlength)
{
  if (prefixindex < prefixlength-1)
  {
    Codetype ordercode = (code - bcktab->filltable[prefixindex])/
                         (bcktab->filltable[prefixindex]+1);
    bcktab->distpfxidx[prefixindex-1][ordercode]++;
  }
  bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
}

int bcktab2file(FILE *fp,
                const Bcktab *bcktab,
                Error *err)
{
  if (fwrite(bcktab->leftborder,
             sizeof (*bcktab->leftborder),
             (size_t) (bcktab->numofallcodes+1),
             fp)
             != (size_t) (bcktab->numofallcodes+1))
  {
    error_set(err,"cannot write %u items of size %u: errormsg=\"%s\"",
              bcktab->numofallcodes+1,
              (unsigned int) sizeof (*bcktab->leftborder),
              strerror(errno));
    return -1;
  }
  if (fwrite(bcktab->countspecialcodes,
             sizeof (*bcktab->countspecialcodes),
             (size_t) bcktab->numofspecialcodes,
             fp)
             != (size_t) bcktab->numofspecialcodes)
  {
    error_set(err,"cannot write %u items of size %u: errormsg=\"%s\"",
              bcktab->numofspecialcodes,
              (unsigned int) sizeof (*bcktab->countspecialcodes),
              strerror(errno));
    return -2;
  }
  return 0;
}
