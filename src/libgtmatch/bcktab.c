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
#include "libgtcore/symboldef.h"
#include "esafileend.h"
#include "spacedef.h"
#include "bcktab.h"

#include "initbasepower.pr"

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS))

/*
  ADD prefixlength
*/

struct Bcktab
{
  Seqpos totallength,
         *leftborder,
         *countspecialcodes;
  Codetype numofallcodes,
           numofspecialcodes,
           **multimappower,
           *basepower,
           *filltable;
  unsigned long **distpfxidx;
};

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

void freebcktab(Bcktab **bcktab,bool mapped)
{
  Bcktab *bcktabptr = *bcktab;
  if (mapped)
  {
    fa_xmunmap((void *) bcktabptr->leftborder);
    bcktabptr->leftborder = NULL;
    bcktabptr->countspecialcodes = NULL;
  } else
  {
    FREESPACE(bcktabptr->leftborder);
    FREESPACE(bcktabptr->countspecialcodes);
  }
  if (bcktabptr->multimappower != NULL)
  {
    multimappowerfree(&bcktabptr->multimappower);
  }
  if (bcktabptr->distpfxidx != NULL)
  {
    FREESPACE(bcktabptr->distpfxidx[0]);
    FREESPACE(bcktabptr->distpfxidx);
  }
  FREESPACE(bcktabptr->filltable);
  FREESPACE(bcktabptr->basepower);
  FREESPACE(*bcktab);
}

static void initbcktabwithNULL(Bcktab *bcktab)
{
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->multimappower = NULL;
  bcktab->filltable = NULL;
  bcktab->basepower = NULL;
  bcktab->distpfxidx = NULL;
}

Bcktab *mapbcktab(const Str *indexname,
                  Seqpos totallength,
                  unsigned int numofchars,
                  unsigned int prefixlength,
                  Error *err)
{
  Bcktab *bcktab;
  bool haserr = false;

  ALLOCASSIGNSPACE(bcktab,NULL,Bcktab,1);
  bcktab->totallength = totallength;
  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->filltable = initfilltable(bcktab->basepower,prefixlength);
  bcktab->multimappower = initmultimappower(numofchars,prefixlength);
  bcktab->distpfxidx = NULL;
  bcktab->leftborder
    = genericmaptable(indexname,
                      BCKTABSUFFIX,
                      (Seqpos) (bcktab->numofallcodes + 1 + 
                                bcktab->numofspecialcodes),
                      sizeof (Seqpos),
                      err);
  if (bcktab->leftborder == NULL)
  {
    bcktab->countspecialcodes = NULL;
    haserr = true;
  } else
  {
    bcktab->countspecialcodes = bcktab->leftborder + bcktab->numofallcodes + 1;
  }
  if (haserr)
  {
    freebcktab(&bcktab,false);
    return NULL;
  }
  return bcktab;
}

static unsigned long **initdistprefixindexcounts(const Codetype *basepower,
                                                 unsigned int prefixlength)
{
  if (prefixlength > 1U)
  {
    unsigned int idx;
    unsigned long *counters, numofcounters, **distpfxidx;

    for (numofcounters = 0, idx=1U; idx <= prefixlength-1; idx++)
    {
      numofcounters += basepower[idx];
    }
    assert(numofcounters > 0);
    ALLOCASSIGNSPACE(distpfxidx,NULL,unsigned long *,prefixlength-1);
    ALLOCASSIGNSPACE(counters,NULL,unsigned long,numofcounters);
    memset(counters,0,(size_t) sizeof (*counters) * numofcounters);
    distpfxidx[0] = counters;
    for (idx=1U; idx<prefixlength-1; idx++)
    {
      distpfxidx[idx] = distpfxidx[idx-1] + basepower[idx];
    }
    return distpfxidx;
  } else
  {
    return NULL;
  }
}

Bcktab *allocBcktab(Seqpos totallength,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    unsigned int codebits,
                    unsigned int maxcodevalue,
                    Error *err)
{
  Bcktab *bcktab;
  bool haserr = false;

  ALLOCASSIGNSPACE(bcktab,NULL,Bcktab,1);
  initbcktabwithNULL(bcktab);
  bcktab->totallength = totallength;
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
  if (haserr)
  {
    freebcktab(&bcktab,false);
    return NULL;
  }
  return bcktab;
}

void updatebckspecials(Bcktab *bcktab,
                       Codetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex,
                       unsigned int prefixlength)
{
  Codetype ordercode = (code - bcktab->filltable[prefixindex])/
                       (bcktab->filltable[prefixindex]+1);
  assert(prefixindex > 0);
  if (code == 0)
  {
    printf("increment [prefixindex=%u,pos=%u]\n",prefixindex,ordercode);
  }
  bcktab->distpfxidx[prefixindex-1][ordercode]++;
  if (prefixindex == prefixlength-1)
  {
    assert(FROMCODE2SPECIALCODE(code,numofchars) == ordercode);
  }
  bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
}

void addfinalbckspecials(Bcktab *bcktab,unsigned int numofchars,
                         Seqpos specialcharacters)
{
  Codetype specialcode;

  specialcode = FROMCODE2SPECIALCODE(bcktab->filltable[0],numofchars);
  bcktab->countspecialcodes[specialcode] += specialcharacters + 1;
}

static long fromcode2countspecialcodes(Codetype code,
                                       unsigned int prefixlength,
                                       const Bcktab *bcktab)
{
  if (code >= bcktab->filltable[prefixlength-1])
  {
    Codetype ordercode = code - bcktab->filltable[prefixlength-1];
    Codetype divisor = bcktab->filltable[prefixlength-1] + 1;
    if (ordercode % divisor == 0)
    {
      ordercode /= divisor;
      return bcktab->countspecialcodes[ordercode];
    }
  }
  return 0;
}

static void pfxidxpartialsums(unsigned long *count,
                              Codetype code,
                              const Bcktab *bcktab,
                              unsigned int prefixlength)
{
  unsigned int prefixindex;
  unsigned long sum = 0; 
  Codetype ordercode, divisor;

  count[0] = 0;
  for (prefixindex=1U; prefixindex<=prefixlength-1; prefixindex++)
  {
    count[prefixindex] = 0;
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        count[prefixindex] = bcktab->distpfxidx[prefixindex-1][ordercode];
        sum += count[prefixindex];
        if (prefixindex == prefixlength -1)
        {
          if (sum != bcktab->countspecialcodes[ordercode])
          {
            fprintf(stderr,"code %u: sum = %lu != %u = specialsinbucket\n",
                    code,sum,bcktab->countspecialcodes[ordercode]);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }
  assert(sum == fromcode2countspecialcodes(code,prefixlength,bcktab));
  if (prefixlength > 2U)
  {
    for (prefixindex = prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      count[prefixindex] += count[prefixindex+1];
    }
    if (sum != count[1])
    {
      fprintf(stderr,"code %u: sum = %lu != %lu = count[1]\n",
              code,sum,count[1]);
      exit(EXIT_FAILURE);
    }
  }
}

void checkcountspecialcodes(const Bcktab *bcktab,unsigned int prefixlength)
{
  Codetype code;
  unsigned long *count;

  if (prefixlength >= 2U)
  {
    ALLOCASSIGNSPACE(count,NULL,unsigned long,prefixlength);
    for (code=0; code<bcktab->numofallcodes; code++)
    {
      pfxidxpartialsums(count,
                        code,
                        bcktab,
                        prefixlength);
    }
    FREESPACE(count);
  }
}

Codetype codedownscale(const Bcktab *bcktab,
                       Codetype code,
                       unsigned int prefixindex,
                       unsigned int maxprefixlen)
{
  unsigned int remain;

  code -= bcktab->filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (bcktab->filltable[remain]+1);
  code *= bcktab->basepower[remain];
  code += bcktab->filltable[prefixindex];
  return code;
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

unsigned int calcbucketboundsparts(Bucketspecification *bucketspec,
                                   const Bcktab *bcktab,
                                   Codetype code,
                                   Codetype maxcode,
                                   Seqpos totalwidth,
                                   unsigned int rightchar,
                                   unsigned int numofchars)
{
  bucketspec->left = bcktab->leftborder[code];
  if (code == maxcode)
  {
    assert(totalwidth >= bucketspec->left);
    bucketspec->nonspecialsinbucket
      = (unsigned long) (totalwidth - bucketspec->left);
  } else
  {
    if (bcktab->leftborder[code+1] > 0)
    {
      bucketspec->nonspecialsinbucket
        = (unsigned long) (bcktab->leftborder[code+1] - bucketspec->left);
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
  }
  assert(rightchar == code % numofchars);
  if (rightchar == numofchars - 1)
  {
    bucketspec->specialsinbucket
      = (unsigned long)
        bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
    if (bucketspec->nonspecialsinbucket >= bucketspec->specialsinbucket)
    {
      bucketspec->nonspecialsinbucket -= bucketspec->specialsinbucket;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
    rightchar = 0;
  } else
  {
    bucketspec->specialsinbucket = 0;
    rightchar++;
  }
  return rightchar;
}

void calcbucketboundaries(Bucketspecification *bucketspec,
                          const Bcktab *bcktab,
                          Codetype code)
{
  unsigned int numofchars = bcktab->basepower[1];
  assert(code != bcktab->numofallcodes);
  (void) calcbucketboundsparts(bucketspec,
                               bcktab,
                               code,
                               bcktab->numofallcodes, /* code < numofallcodes */
                               0,                     /* not necessary */
                               code % numofchars,
                               numofchars);
}

unsigned int pfxidx2lcpvalues(Uchar *lcpsubtab,
                              unsigned long specialsinbucket,
                              const Bcktab *bcktab,
                              Codetype code,
                              unsigned int prefixlength)
{
  unsigned int prefixindex, maxvalue = 0;
  Codetype ordercode, divisor;
  unsigned long idx, insertindex = specialsinbucket-1;

  for (prefixindex=1U; prefixindex<prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0 &&
            insertindex > 0)
        {
          for (idx=0;
               insertindex > 0 &&
               idx < bcktab->distpfxidx[prefixindex-1][ordercode];
               idx++, insertindex--)
          {
            lcpsubtab[insertindex] = (Uchar) prefixindex;
          }
          if (maxvalue < prefixindex)
          {
            maxvalue = prefixindex;
          }
        }
      }
    }
  }
  while (insertindex > 0)
  {
    lcpsubtab[insertindex] = (Uchar) (prefixlength-1);
    if (maxvalue < prefixindex)
    {
      maxvalue = prefixindex;
    }
    insertindex--;
  }
  return maxvalue;
}

const Codetype **bcktab_multimappower(const Bcktab *bcktab)
{
  return (const Codetype **) bcktab->multimappower;
}

Codetype bcktab_filltable(const Bcktab *bcktab,unsigned int idx)
{
  return bcktab->filltable[idx];
}

Seqpos *bcktab_leftborder(Bcktab *bcktab)
{
  return bcktab->leftborder;
}

Codetype bcktab_numofallcodes(Bcktab *bcktab)
{
  return bcktab->numofallcodes;
}
