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
#include "libgtcore/chardef.h"
#include "esafileend.h"
#include "mapspec-def.h"
#include "spacedef.h"
#include "bcktab.h"

#include "initbasepower.pr"

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS))

struct Bcktab
{
  Seqpos totallength,
         *leftborder;
  Codetype numofallcodes,
           numofspecialcodes,
           **multimappower,
           *basepower,
           *filltable;
  unsigned int prefixlength;
  unsigned long sizeofrep,
                *countspecialcodes,
                **distpfxidx;
  Uchar *qgrambuffer;
  void *mappedptr;
};

static unsigned long numofdistpfxidxcounters(const Codetype *basepower,
                                             unsigned int prefixlength)
{
  if (prefixlength > 2U)
  {
    unsigned long numofcounters = 0;
    unsigned int idx;

    for (idx=1U; idx < prefixlength-1; idx++)
    {
      numofcounters += basepower[idx];
    }
    return numofcounters;
  }
  return 0;
}

unsigned long sizeofbuckettable(unsigned int numofchars,
                                unsigned int prefixlength)
{
  unsigned long sizeofrep;
  Codetype *basepower, numofallcodes, numofspecialcodes;

  basepower = initbasepower(numofchars,prefixlength);
  numofallcodes = basepower[prefixlength];
  numofspecialcodes = basepower[prefixlength-1];
  sizeofrep
    = (unsigned long)
      sizeof (Seqpos) * (numofallcodes + 1) +
      sizeof (Seqpos) * numofspecialcodes +
      sizeof (unsigned long) * numofdistpfxidxcounters(basepower,
                                                       prefixlength);
  FREESPACE(basepower);
  return sizeofrep;
}

static void setdistpfxidxptrs(unsigned long **distpfxidx,
                              unsigned long *ptr,
                              const Codetype *basepower,
                              unsigned int prefixlength)
{
  unsigned int idx;

  distpfxidx[0] = ptr;
  for (idx=1U; idx<prefixlength-1; idx++)
  {
    distpfxidx[idx] = distpfxidx[idx-1] + basepower[idx];
  }
}

static unsigned long **allocdistprefixindexcounts(const Codetype *basepower,
                                                  unsigned int prefixlength,
                                                  Verboseinfo *verboseinfo)
{
  if (prefixlength > 2U)
  {
    unsigned long numofcounters;

    numofcounters = numofdistpfxidxcounters(basepower,prefixlength);
    if (numofcounters > 0)
    {
      unsigned long *counters, **distpfxidx;

      ALLOCASSIGNSPACE(distpfxidx,NULL,unsigned long *,prefixlength-1);
      ALLOCASSIGNSPACE(counters,NULL,unsigned long,numofcounters);
      showverbose(verboseinfo,"sizeof (distpfxidx)=%lu",
                  (unsigned long) sizeof (*distpfxidx) * (prefixlength-1) +
                                  sizeof (*counters) * numofcounters);
      memset(counters,0,(size_t) sizeof (*counters) * numofcounters);
      setdistpfxidxptrs(distpfxidx,counters,basepower,prefixlength);
      return distpfxidx;
    }
  }
  return NULL;
}

static Bcktab *newBcktab(unsigned int numofchars,
                         unsigned int prefixlength,
                         Seqpos totallength)
{
  Bcktab *bcktab;

  ALLOCASSIGNSPACE(bcktab,NULL,Bcktab,1);
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->distpfxidx = NULL;
  bcktab->totallength = totallength;
  bcktab->mappedptr = NULL;
  bcktab->prefixlength = prefixlength;
  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->filltable = initfilltable(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->multimappower = initmultimappower(numofchars,prefixlength);
  ALLOCASSIGNSPACE(bcktab->qgrambuffer,NULL,Uchar,prefixlength);
  bcktab->sizeofrep
    = (unsigned long)
      sizeof (*bcktab->leftborder) * (bcktab->numofallcodes + 1) +
      sizeof (*bcktab->countspecialcodes) * bcktab->numofspecialcodes +
      sizeof (unsigned long) * numofdistpfxidxcounters(bcktab->basepower,
                                                       bcktab->prefixlength);
  return bcktab;
}

Bcktab *allocBcktab(Seqpos totallength,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    unsigned int codebits,
                    Codetype maxcodevalue,
                    Verboseinfo *verboseinfo,
                    Error *err)
{
  Bcktab *bcktab;
  bool haserr = false;

  bcktab = newBcktab(numofchars,prefixlength,totallength);
  if (maxcodevalue > 0 && bcktab->numofallcodes-1 > maxcodevalue)
  {
    error_set(err,"alphasize^prefixlength-1 = " FormatCodetype
                  " does not fit into %u"
                  " bits: choose smaller value for prefixlength",
                  bcktab->numofallcodes-1,
                  codebits);
    haserr = true;
  }
  if (!haserr)
  {
    ALLOCASSIGNSPACE(bcktab->leftborder,NULL,Seqpos,
                     bcktab->numofallcodes+1);
    showverbose(verboseinfo,"sizeof (leftborder)=%lu",
              (unsigned long) sizeof (*bcktab->leftborder) *
                              (bcktab->numofallcodes+1));
    memset(bcktab->leftborder,0,
           sizeof (*bcktab->leftborder) *
           (size_t) (bcktab->numofallcodes+1));
    ALLOCASSIGNSPACE(bcktab->countspecialcodes,NULL,unsigned long,
                     bcktab->numofspecialcodes);
    showverbose(verboseinfo,
                "sizeof (countspecialcodes)=%lu",
                (unsigned long) sizeof (*bcktab->countspecialcodes) *
                bcktab->numofspecialcodes);
    memset(bcktab->countspecialcodes,0,
           sizeof (*bcktab->countspecialcodes) *
                  (size_t) bcktab->numofspecialcodes);
    bcktab->distpfxidx = allocdistprefixindexcounts(bcktab->basepower,
                                                    prefixlength,
                                                    verboseinfo);
  }
  if (haserr)
  {
    freebcktab(&bcktab);
    return NULL;
  }
  return bcktab;
}

static void assignbcktabmapspecification(ArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode)
{
  Bcktab *bcktab = (Bcktab *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long numofcounters;

  NEWMAPSPEC(bcktab->leftborder,Seqpos,
             (unsigned long) (bcktab->numofallcodes+1));
  NEWMAPSPEC(bcktab->countspecialcodes,Unsignedlong,
             (unsigned long) bcktab->numofspecialcodes);
  numofcounters = numofdistpfxidxcounters((const Codetype *) bcktab->basepower,
                                          bcktab->prefixlength);
  if (numofcounters > 0)
  {
    if (!writemode)
    {
      ALLOCASSIGNSPACE(bcktab->distpfxidx,NULL,unsigned long *,
                       bcktab->prefixlength-1);
    }
    NEWMAPSPEC(bcktab->distpfxidx[0],Unsignedlong,numofcounters);
  }
}

int bcktab2file(FILE *fp,const Bcktab *bcktab,Error *err)
{
  error_check(err);
  return flushtheindex2file(fp,
                            assignbcktabmapspecification,
                            (Bcktab *) bcktab,
                            bcktab->sizeofrep,
                            err);
}

static int fillbcktabmapspecstartptr(Bcktab *bcktab,
                                     const Str *indexname,
                                     Error *err)
{
  bool haserr = false;
  Str *tmpfilename;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,BCKTABSUFFIX);
  if (fillmapspecstartptr(assignbcktabmapspecification,
                          &bcktab->mappedptr,
                          bcktab,
                          tmpfilename,
                          bcktab->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

Bcktab *mapbcktab(const Str *indexname,
                  Seqpos totallength,
                  unsigned int numofchars,
                  unsigned int prefixlength,
                  Error *err)
{
  Bcktab *bcktab;

  bcktab = newBcktab(numofchars,prefixlength,totallength);
  if (fillbcktabmapspecstartptr(bcktab,
                                indexname,
                                err) != 0)
  {
    freebcktab(&bcktab);
    return NULL;
  }
  if (bcktab->distpfxidx != NULL)
  {
    setdistpfxidxptrs(bcktab->distpfxidx,bcktab->distpfxidx[0],
                      bcktab->basepower,bcktab->prefixlength);
  }
#ifdef SKDEBUG
  checkcountspecialcodes(bcktab);
#endif
  return bcktab;
}

void freebcktab(Bcktab **bcktab)
{
  Bcktab *bcktabptr = *bcktab;

  if (bcktabptr->mappedptr != NULL) /* use mapped file */
  {
    fa_xmunmap(bcktabptr->mappedptr);
    bcktabptr->mappedptr = NULL;
    bcktabptr->leftborder = NULL;
    bcktabptr->countspecialcodes = NULL;
    if (bcktabptr->distpfxidx != NULL)
    {
      bcktabptr->distpfxidx[0] = NULL;
    }
  } else
  {
    FREESPACE(bcktabptr->leftborder);
    FREESPACE(bcktabptr->countspecialcodes);
    if (bcktabptr->distpfxidx != NULL)
    {
      FREESPACE(bcktabptr->distpfxidx[0]);
    }
  }
  FREESPACE(bcktabptr->distpfxidx);
  if (bcktabptr->multimappower != NULL)
  {
    multimappowerfree(&bcktabptr->multimappower);
  }
  FREESPACE(bcktabptr->filltable);
  FREESPACE(bcktabptr->basepower);
  FREESPACE(bcktabptr->qgrambuffer);
  FREESPACE(*bcktab);
}

void updatebckspecials(Bcktab *bcktab,
                       Codetype code,
                       unsigned int numofchars,
                       unsigned int prefixindex)
{
  assert(prefixindex > 0);
  if (prefixindex < bcktab->prefixlength-1)
  {
    Codetype ordercode = (code - bcktab->filltable[prefixindex])/
                         (bcktab->filltable[prefixindex]+1);
    bcktab->distpfxidx[prefixindex-1][ordercode]++;
  }
  bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)]++;
}

void addfinalbckspecials(Bcktab *bcktab,unsigned int numofchars,
                         Seqpos specialcharacters)
{
  Codetype specialcode;

  specialcode = FROMCODE2SPECIALCODE(bcktab->filltable[0],numofchars);
  bcktab->countspecialcodes[specialcode]
    += (unsigned long) (specialcharacters + 1);
}

#ifdef SKDEBUG
static unsigned long fromcode2countspecialcodes(Codetype code,
                                                const Bcktab *bcktab)
{
  if (code >= bcktab->filltable[bcktab->prefixlength-1])
  {
    Codetype ordercode = code - bcktab->filltable[bcktab->prefixlength-1];
    Codetype divisor = bcktab->filltable[bcktab->prefixlength-1] + 1;
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
                              const Bcktab *bcktab)
{
  unsigned int prefixindex;
  unsigned long sum = 0, specialsinbucket;
  Codetype ordercode, divisor;

  memset(count,0,sizeof (*count) * (size_t) bcktab->prefixlength);
  for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        count[prefixindex] = bcktab->distpfxidx[prefixindex-1][ordercode];
        sum += count[prefixindex];
      }
    } else
    {
      break;
    }
  }
  specialsinbucket = fromcode2countspecialcodes(code,bcktab);
  assert(sum <= specialsinbucket);
  count[bcktab->prefixlength-1] = specialsinbucket - sum;
  if (bcktab->prefixlength > 2U)
  {
    for (prefixindex = bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      count[prefixindex] += count[prefixindex+1];
    }
    if (specialsinbucket != count[1])
    {
      fprintf(stderr,"code " FormatCodetype ": sum = %lu != %lu = count[1]\n",
              code,sum,count[1]);
      exit(EXIT_FAILURE);
    }
  }
}

void checkcountspecialcodes(const Bcktab *bcktab)
{
  Codetype code;
  unsigned long *count;

  if (bcktab->prefixlength >= 2U)
  {
    ALLOCASSIGNSPACE(count,NULL,unsigned long,bcktab->prefixlength);
    for (code=0; code<bcktab->numofallcodes; code++)
    {
      pfxidxpartialsums(count,
                        code,
                        bcktab);
    }
    FREESPACE(count);
  }
}
#endif

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
  assert(rightchar == (unsigned int) (code % numofchars));
  if (rightchar == numofchars - 1)
  {
    bucketspec->specialsinbucket
      = bcktab->countspecialcodes[FROMCODE2SPECIALCODE(code,numofchars)];
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
  unsigned int numofchars = (unsigned int) bcktab->basepower[1];
  assert(code != bcktab->numofallcodes);
  (void) calcbucketboundsparts(bucketspec,
                               bcktab,
                               code,
                               bcktab->numofallcodes, /* code < numofallcodes */
                               0,                     /* not necessary */
                               (unsigned int) (code % numofchars),
                               numofchars);
}

unsigned int singletonmaxprefixindex(const Bcktab *bcktab,Codetype code)
{
  if (bcktab->prefixlength > 2U)
  {
    Codetype ordercode, divisor;
    unsigned int prefixindex;

    for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      if (code >= bcktab->filltable[prefixindex])
      {
        ordercode = code - bcktab->filltable[prefixindex];
        divisor = bcktab->filltable[prefixindex] + 1;
        if (ordercode % divisor == 0)
        {
          ordercode /= divisor;
          if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0)
          {
            return prefixindex;
          }
        }
      } else
      {
        break;
      }
    }
  }
  return bcktab->prefixlength-1;
}

unsigned int pfxidx2lcpvalues(unsigned int *minprefixindex,
                              Uchar *lcpsubtab,
                              unsigned long specialsinbucket,
                              const Bcktab *bcktab,
                              Codetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  Codetype ordercode, divisor;
  unsigned long idx;
  Uchar *insertptr;

  *minprefixindex = bcktab->prefixlength;
  insertptr = lcpsubtab + specialsinbucket - 1;
  for (prefixindex=1U; prefixindex<bcktab->prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        if (bcktab->distpfxidx[prefixindex-1][ordercode] > 0)
        {
          maxprefixindex = prefixindex;
          if (*minprefixindex > prefixindex)
          {
            *minprefixindex = prefixindex;
          }
          for (idx=0; idx < bcktab->distpfxidx[prefixindex-1][ordercode]; idx++)
          {
            assert(insertptr >= lcpsubtab);
            *insertptr-- = (Uchar) prefixindex;
          }
        }
      }
    }
  }
  if (insertptr >= lcpsubtab)
  {
    maxprefixindex = bcktab->prefixlength-1;
    if (*minprefixindex == bcktab->prefixlength)
    {
      *minprefixindex = bcktab->prefixlength-1;
    }
    while (insertptr >= lcpsubtab)
    {
      *insertptr-- = (Uchar) (bcktab->prefixlength-1);
    }
  }
  return maxprefixindex;
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

Codetype bcktab_numofallcodes(const Bcktab *bcktab)
{
  return bcktab->numofallcodes;
}

#ifdef SKDEBUG
#include "qgram2code.h"
void consistencyofsuffix(int line,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         const Bcktab *bcktab,
                         unsigned int numofchars,
                         const Suffixwithcode *suffix)
{
  unsigned int idx, firstspecial = bcktab->prefixlength, gramfirstspecial;
  Codetype qgramcode = 0;
  Uchar cc = 0;

  for (idx=0; idx<bcktab->prefixlength; idx++)
  {
    if ((Seqpos) (suffix->startpos + idx) >= bcktab->totallength)
    {
      firstspecial = idx;
      break;
    }
    cc = getencodedchar(encseq,suffix->startpos + idx, readmode);
    if (ISSPECIAL(cc))
    {
      firstspecial = idx;
      break;
    }
    bcktab->qgrambuffer[idx] = cc;
  }
  for (idx=firstspecial; idx<bcktab->prefixlength; idx++)
  {
    bcktab->qgrambuffer[idx] = (Uchar) (numofchars-1);
  }
  gramfirstspecial = qgram2code(&qgramcode,
                                (const Codetype **) bcktab->multimappower,
                                bcktab->prefixlength,
                                bcktab->qgrambuffer);
  assert(gramfirstspecial == bcktab->prefixlength);
  assert(qgramcode == suffix->code);
  if (firstspecial != suffix->prefixindex)
  {
    fprintf(stderr,"line %d: code=%u: ",line,suffix->code);
    fprintf(stderr,"firstspecial = %u != %u = suffix->prefixindex\n",
                    firstspecial,suffix->prefixindex);
    exit(EXIT_FAILURE);
  }
}
#endif
