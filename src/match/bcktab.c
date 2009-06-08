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
#include <math.h>
#include "core/error.h"
#include "core/str.h"
#include "core/fa.h"
#include "core/symboldef.h"
#include "core/chardef.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "esa-fileend.h"
#include "mapspec-def.h"
#include "spacedef.h"
#include "bcktab.h"
#include "format64.h"
#include "initbasepower.h"

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS))

typedef struct
{
  unsigned long nonspecialsmaxbucketsize,
                specialsmaxbucketsize,
                maxbucketsize,
                log2nonspecialbucketsizedist[GT_MAXLOG2VALUE+1],
                log2specialbucketsizedist[GT_MAXLOG2VALUE+1];
} Maxbucketinfo;

struct Bcktab
{
  Seqpos *leftborder;
  Codetype numofallcodes,
           numofspecialcodes,
           **multimappower,
           *basepower,
           *filltable;
  unsigned int prefixlength;
  unsigned long sizeofrep,
                *countspecialcodes,
                **distpfxidx;
  GtUchar *qgrambuffer;
  Maxbucketinfo maxbucketinfo;
  unsigned int optimalnumofbits;
  unsigned short logofremaining;
  bool allocated;
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

uint64_t sizeofbuckettable(unsigned int numofchars,
                           unsigned int prefixlength)
{
  uint64_t sizeofrep;
  double numofallcodes, numofspecialcodes;
  Codetype *basepower;

  numofallcodes = pow((double) numofchars,(double) prefixlength);
  numofspecialcodes = pow((double) numofchars,(double) (prefixlength-1));
  if (prefixlength >= 2U)
  {
    basepower = initbasepower(numofchars,prefixlength-2);
  } else
  {
    basepower = NULL;
  }
  sizeofrep
    = (uint64_t)
      sizeof (Seqpos) * (numofallcodes + 1.0) +
      sizeof (Seqpos) * numofspecialcodes +
      sizeof (unsigned long) * numofdistpfxidxcounters(basepower,prefixlength);
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
                                                  unsigned int prefixlength)
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
      /*
      showverbose(verboseinfo,"sizeof (distpfxidx)=%lu",
                  (unsigned long) sizeof (*distpfxidx) * (prefixlength-1) +
                                  sizeof (*counters) * numofcounters);
      */
      memset(counters,0,(size_t) sizeof (*counters) * numofcounters);
      setdistpfxidxptrs(distpfxidx,counters,basepower,prefixlength);
      return distpfxidx;
    }
  }
  return NULL;
}

static Bcktab *newBcktab(unsigned int numofchars,
                         unsigned int prefixlength)
{
  Bcktab *bcktab;

  ALLOCASSIGNSPACE(bcktab,NULL,Bcktab,1);
  bcktab->leftborder = NULL;
  bcktab->countspecialcodes = NULL;
  bcktab->distpfxidx = NULL;
  bcktab->mappedptr = NULL;
  bcktab->prefixlength = prefixlength;
  bcktab->basepower = initbasepower(numofchars,prefixlength);
  bcktab->filltable = initfilltable(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->multimappower = initmultimappower(numofchars,prefixlength);
  bcktab->allocated = false;
  ALLOCASSIGNSPACE(bcktab->qgrambuffer,NULL,GtUchar,prefixlength);
  bcktab->sizeofrep
    = (unsigned long)
      sizeof (*bcktab->leftborder) * (bcktab->numofallcodes + 1) +
      sizeof (*bcktab->countspecialcodes) * bcktab->numofspecialcodes +
      sizeof (unsigned long) * numofdistpfxidxcounters(bcktab->basepower,
                                                       bcktab->prefixlength);
  return bcktab;
}

Bcktab *allocBcktab(unsigned int numofchars,
                    unsigned int prefixlength,
                    bool storespecialcodes,
                    Verboseinfo *verboseinfo,
                    GtError *err)
{
  Bcktab *bcktab;
  bool haserr = false;

  bcktab = newBcktab(numofchars,prefixlength);
  bcktab->allocated = true;
  if (storespecialcodes && bcktab->numofallcodes > 0 &&
      bcktab->numofallcodes-1 > (unsigned long) MAXCODEVALUE)
  {
    gt_error_set(err,"alphasize^prefixlength-1 = " FormatCodetype
                  " does not fit into %u"
                  " bits: choose smaller value for prefixlength",
                  bcktab->numofallcodes-1,
                  CODEBITS);
    haserr = true;
  }
  if (!haserr)
  {
    ALLOCASSIGNSPACE(bcktab->leftborder,NULL,Seqpos,
                     bcktab->numofallcodes+1);
    /*
    showverbose(verboseinfo,"sizeof (leftborder)=%lu",
              (unsigned long) sizeof (*bcktab->leftborder) *
                              (bcktab->numofallcodes+1));
    */
    memset(bcktab->leftborder,0,
           sizeof (*bcktab->leftborder) *
           (size_t) (bcktab->numofallcodes+1));
    ALLOCASSIGNSPACE(bcktab->countspecialcodes,NULL,unsigned long,
                     bcktab->numofspecialcodes);
    /*
    showverbose(verboseinfo,
                "sizeof (countspecialcodes)=%lu",
                (unsigned long) sizeof (*bcktab->countspecialcodes) *
                bcktab->numofspecialcodes);
    */
    memset(bcktab->countspecialcodes,0,
           sizeof (*bcktab->countspecialcodes) *
                  (size_t) bcktab->numofspecialcodes);
    bcktab->distpfxidx = allocdistprefixindexcounts(bcktab->basepower,
                                                    prefixlength);
    showverbose(verboseinfo,"sizeof (bcktab)=" Formatuint64_t " bytes",
                PRINTuint64_tcast(sizeofbuckettable(numofchars,prefixlength)));
  }
  if (haserr)
  {
    bcktab_delete(&bcktab);
    return NULL;
  }
  return bcktab;
}

static void assignbcktabmapspecification(GtArrayMapspecification *mapspectable,
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

int bcktab2file(FILE *fp,const Bcktab *bcktab,GtError *err)
{
  gt_error_check(err);
  return flushtheindex2file(fp,
                            assignbcktabmapspecification,
                            (Bcktab *) bcktab,
                            bcktab->sizeofrep,
                            err);
}

static int fillbcktabmapspecstartptr(Bcktab *bcktab,
                                     const GtStr *indexname,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,BCKTABSUFFIX);
  if (fillmapspecstartptr(assignbcktabmapspecification,
                          &bcktab->mappedptr,
                          bcktab,
                          tmpfilename,
                          bcktab->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

Bcktab *mapbcktab(const GtStr *indexname,
                  unsigned int numofchars,
                  unsigned int prefixlength,
                  GtError *err)
{
  Bcktab *bcktab;

  bcktab = newBcktab(numofchars,prefixlength);
  bcktab->allocated = false;
  if (fillbcktabmapspecstartptr(bcktab,
                                indexname,
                                err) != 0)
  {
    bcktab_delete(&bcktab);
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

void showbcktab(const Bcktab *bcktab)
{
  unsigned int prefixindex;

  for (prefixindex=1U; prefixindex < bcktab->prefixlength-1; prefixindex++)
  {
    Codetype code;
    unsigned long sum = 0;
    for (code = 0; code < bcktab->basepower[prefixindex]; code++)
    {
      sum += bcktab->distpfxidx[prefixindex-1][code];
      printf("distpfxidx[%u][%lu]=%lu\n",
              prefixindex,code,bcktab->distpfxidx[prefixindex-1][code]);
    }
    printf("sum %lu\n",sum);
  }
}

void bcktab_delete(Bcktab **bcktab)
{
  Bcktab *bcktabptr = *bcktab;

  /*showbcktab(bcktabptr);*/
  if (bcktabptr->allocated)
  {
    FREESPACE(bcktabptr->leftborder);
    FREESPACE(bcktabptr->countspecialcodes);
    if (bcktabptr->distpfxidx != NULL)
    {
      FREESPACE(bcktabptr->distpfxidx[0]);
    }
  } else
  {
    if (bcktabptr->mappedptr != NULL)
    {
      gt_fa_xmunmap(bcktabptr->mappedptr);
    }
    bcktabptr->mappedptr = NULL;
    bcktabptr->leftborder = NULL;
    bcktabptr->countspecialcodes = NULL;
    if (bcktabptr->distpfxidx != NULL)
    {
      bcktabptr->distpfxidx[0] = NULL;
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
  gt_assert(prefixindex > 0);
  gt_assert(prefixindex <= bcktab->prefixlength);
  gt_assert(code < bcktab->numofallcodes);
  if (prefixindex < bcktab->prefixlength-1)
  {
    Codetype ordercode = (code - bcktab->filltable[prefixindex])/
                         (bcktab->filltable[prefixindex]+1);
    /*
    printf("prefixindex=%u\n",(unsigned int) prefixindex);
    printf("code=%u\n",(unsigned int) code);
    printf("ordercode=%u\n",(unsigned int) ordercode);
    */
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
  gt_assert(sum <= specialsinbucket);
  count[bcktab->prefixlength-1] = specialsinbucket - sum;
  if (bcktab->prefixlength > 2U)
  {
    for (prefixindex = bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      count[prefixindex] += count[prefixindex+1];
    }
#ifndef NDEBUG
    if (specialsinbucket != count[1])
    {
      fprintf(stderr,"code " FormatCodetype ": sum = %lu != %lu = count[1]\n",
              code,sum,count[1]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
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
    gt_assert(totalwidth >= bucketspec->left);
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
  gt_assert(rightchar == (unsigned int) (code % numofchars));
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
  gt_assert(code != bcktab->numofallcodes);
  (void) calcbucketboundsparts(bucketspec,
                               bcktab,
                               code,
                               bcktab->numofallcodes, /* code < numofallcodes */
                               0,                     /* not necessary */
                               (unsigned int) (code % numofchars),
                               numofchars);
}

static void updatelog2values(unsigned long *tab,unsigned long maxvalue)
{
  unsigned long multi = 1UL, idx, sum = 1UL;

  for (idx = 0; /* Nothing */; idx++)
  {
    tab[idx] += multi;
    multi *= 2;
    if (sum+multi > maxvalue)
    {
      tab[idx+1] += maxvalue - sum;
      break;
    }
    sum += multi;
  }
}

static unsigned int calc_optimalnumofbits(unsigned short *logofremaining,
                                          const unsigned long *log2tab,
                                          Seqpos totallength)
{
  unsigned int lastbitset = 0, maxbits, optbits = 0;
  unsigned long currentsum = 0, total = 0, optcurrentsum = 0;
  unsigned short tmplogofremaining;
  const size_t size_entry = sizeof (uint32_t) +
                            sizeof (Seqpos) +
                            sizeof (unsigned long);
  size_t savedbitsinbytes,
         hashtablesize, saved, maxsaved = 0;

  for (maxbits = 0; maxbits <= (unsigned int) GT_MAXLOG2VALUE; maxbits++)
  {
    if (log2tab[maxbits] > 0)
    {
      total += log2tab[maxbits];
      lastbitset = maxbits;
    }
  }
#ifdef WITHsave
  printf("lastbitset=%u\n",lastbitset);
#endif
  for (maxbits = 0; maxbits <= lastbitset; maxbits++)
  {
    if (log2tab[maxbits] > 0)
    {
      currentsum += log2tab[maxbits];
      tmplogofremaining = (unsigned short) gt_determinebitspervalue(
                                           (uint64_t) (total - currentsum));
      hashtablesize = ((1 << tmplogofremaining)-1) * size_entry;
      savedbitsinbytes = (totallength/CHAR_BIT) * (lastbitset - maxbits + 1);
#ifdef WITHsave
      printf("savedbitsintbytes=%lu,hashtablesize=%lu\n",
              (unsigned long) savedbitsinbytes,(unsigned long) hashtablesize);
#endif
      if (savedbitsinbytes > hashtablesize)
      {
        saved = savedbitsinbytes - hashtablesize;
#ifdef WITHsave
        printf("saved=%lu\n",(unsigned long) saved);
#endif
        if (saved > maxsaved)
        {
          maxsaved = saved;
#ifdef WITHsave
          printf("maxsaved=%lu\n",(unsigned long) maxsaved);
#endif
          optcurrentsum = currentsum;
          optbits = maxbits;
        }
      }
    }
  }
  *logofremaining = (unsigned short)
                    gt_determinebitspervalue((uint64_t)
                                             (total - optcurrentsum));
#ifdef WITHsave
  printf("store %lu values in hashtable (%lu>=%lu bytes)\n",
         (unsigned long) (total - optcurrentsum),
         (unsigned long) ((1 << (*logofremaining))-1) * size_entry,
         (total - optcurrentsum) * size_entry);
#endif
  return optbits;
}

void determinemaxbucketsize(Bcktab *bcktab,
                            const Codetype mincode,
                            const Codetype maxcode,
                            Seqpos partwidth,
                            unsigned int numofchars,
                            bool hashexceptions,
                            Seqpos totallength) /* relevant for hashexception */
{
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  Bucketspecification bucketspec;
  Codetype code;

  bcktab->maxbucketinfo.specialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.nonspecialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.maxbucketsize = 1UL;
  if (hashexceptions)
  {
    memset(bcktab->maxbucketinfo.log2nonspecialbucketsizedist,0,
           sizeof (*bcktab->maxbucketinfo.log2nonspecialbucketsizedist) *
           (GT_MAXLOG2VALUE+1));
    memset(bcktab->maxbucketinfo.log2specialbucketsizedist,0,
           sizeof (*bcktab->maxbucketinfo.log2specialbucketsizedist) *
           (GT_MAXLOG2VALUE+1));
  }
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      partwidth,
                                      rightchar,
                                      numofchars);
    if (bucketspec.nonspecialsinbucket >
        bcktab->maxbucketinfo.nonspecialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.nonspecialsmaxbucketsize
        = bucketspec.nonspecialsinbucket;
    }
    if (bucketspec.specialsinbucket >
        bcktab->maxbucketinfo.specialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.specialsmaxbucketsize = bucketspec.specialsinbucket;
    }
    if (bucketspec.nonspecialsinbucket + bucketspec.specialsinbucket
        > bcktab->maxbucketinfo.maxbucketsize)
    {
      bcktab->maxbucketinfo.maxbucketsize = bucketspec.nonspecialsinbucket +
                                            bucketspec.specialsinbucket;
    }
    if (hashexceptions)
    {
      if (bucketspec.nonspecialsinbucket >= 1UL)
      {
        updatelog2values(bcktab->maxbucketinfo.log2nonspecialbucketsizedist,
                         bucketspec.nonspecialsinbucket);
      }
      if (bucketspec.specialsinbucket > 1UL)
      {
        updatelog2values(bcktab->maxbucketinfo.log2specialbucketsizedist,
                         bucketspec.specialsinbucket);
      }
    }
  }
  if (hashexceptions)
  {
    bcktab->optimalnumofbits
      = calc_optimalnumofbits(&bcktab->logofremaining,
                              bcktab->maxbucketinfo.
                                      log2nonspecialbucketsizedist,
                              totallength);
  } else
  {
    bcktab->optimalnumofbits = 0;
    bcktab->logofremaining = 0;
  }
  /*
  showverbose(verboseinfo,"maxbucket (specials)=%lu",
              bcktab->maxbucketinfo.specialsmaxbucketsize);
  showverbose(verboseinfo,"maxbucket (nonspecials)=%lu",
              bcktab->maxbucketinfo.nonspecialsmaxbucketsize);
  showverbose(verboseinfo,"maxbucket (all)=%lu",
              bcktab->maxbucketinfo.maxbucketsize);
  */
}

static void showlog2info(const char *tag,const unsigned long *log2tab,
                         Verboseinfo *verboseinfo)
{
  if (verboseinfo != NULL)
  {
    int maxbits;
    unsigned long currentsum = 0, total = 0;
    for (maxbits = 0; maxbits <= GT_MAXLOG2VALUE; maxbits++)
    {
      total += log2tab[maxbits];
    }
    for (maxbits = 0; maxbits <= GT_MAXLOG2VALUE; maxbits++)
    {
      if (log2tab[maxbits] > 0)
      {
        currentsum += log2tab[maxbits];
        showverbose(verboseinfo,"%s[%d]=%lu (%.4f)",tag,maxbits,
                                 log2tab[maxbits],(double) currentsum/total);
      }
    }
    showverbose(verboseinfo,"total=%lu",total);
  }
}

void bcktab_showlog2info(const Bcktab *bcktab,Verboseinfo *verboseinfo)
{
  showlog2info("log2nonspecialbucketsizedist",
                bcktab->maxbucketinfo.log2nonspecialbucketsizedist,
                verboseinfo);
  showlog2info("log2specialbucketsizedist",
                bcktab->maxbucketinfo.log2specialbucketsizedist,
                verboseinfo);
}

unsigned long bcktab_specialsmaxbucketsize(const Bcktab *bcktab)
{
  return bcktab->maxbucketinfo.specialsmaxbucketsize;
}

unsigned long bcktab_nonspecialsmaxbucketsize(const Bcktab *bcktab)
{
  return bcktab->maxbucketinfo.nonspecialsmaxbucketsize;
}

unsigned int bcktab_optimalnumofbits(unsigned short *logofremaining,
                                     const Bcktab *bcktab)
{
  *logofremaining = bcktab->logofremaining;
  return bcktab->optimalnumofbits;
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

unsigned long distpfxidxpartialsums(const Bcktab *bcktab,Codetype code,
                                    unsigned int lowerbound)
{
  Codetype ordercode, divisor;
  unsigned long sum = 0;
  unsigned int prefixindex;

  for (prefixindex=bcktab->prefixlength-2; prefixindex>lowerbound;
       prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        sum += bcktab->distpfxidx[prefixindex-1][ordercode];
      }
    } else
    {
      break;
    }
  }
  return sum;
}

unsigned int pfxidx2lcpvalues(unsigned int *minprefixindex,
                              uint8_t *lcpsubtab,
                              unsigned long specialsinbucket,
                              const Bcktab *bcktab,
                              Codetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  Codetype ordercode, divisor;
  unsigned long idx;
  uint8_t *insertptr;

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
            gt_assert(insertptr >= lcpsubtab);
            *insertptr-- = (uint8_t) prefixindex;
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
      *insertptr-- = (uint8_t) (bcktab->prefixlength-1);
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

void bcktab_leftborderpartialsums(Bcktab *bcktab,Seqpos numofsuffixestosort)
{
  Seqpos *optr;

  for (optr = bcktab->leftborder + 1;
       optr < bcktab->leftborder + bcktab->numofallcodes; optr++)
  {
    *optr += *(optr-1);
  }
  bcktab->leftborder[bcktab->numofallcodes] = numofsuffixestosort;
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
  Seqpos totallength;
  GtUchar cc = 0;

  totallength = getencseqtotallength(encseq);
  for (idx=0; idx<bcktab->prefixlength; idx++)
  {
    if ((Seqpos) (suffix->startpos + idx) >= totallength)
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
    bcktab->qgrambuffer[idx] = (GtUchar) (numofchars-1);
  }
  gramfirstspecial = qgram2code(&qgramcode,
                                (const Codetype **) bcktab->multimappower,
                                bcktab->prefixlength,
                                bcktab->qgrambuffer);
  gt_assert(gramfirstspecial == bcktab->prefixlength);
  gt_assert(qgramcode == suffix->code);
  if (firstspecial != suffix->prefixindex)
  {
    fprintf(stderr,"line %d: code=%u: ",line,suffix->code);
    fprintf(stderr,"firstspecial = %u != %u = suffix->prefixindex\n",
                    firstspecial,suffix->prefixindex);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif
