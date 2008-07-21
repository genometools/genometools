/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/ma.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "defined-types.h"
#include "absdfstrans-imp.h"
#include "initeqsvec.h"

typedef struct
{
  unsigned long prefixofsuffix;
} Parallelmstats;

typedef struct
{
  unsigned long patternlength,
                mstatlength[INTWORDSIZE],
                *eqsvector;
} Matchtaskinfo;

#ifdef SKDEBUG

static void pms_showParallelmstats(const DECLAREPTRDFSSTATE(aliascol),
                                   unsigned long depth,
                                   const void *dfsconstinfo)
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  const Myerscolumn *col = (const Myerscolumn *) aliascol;
  bool first = true;

  unsigned long idx, backmask;

  printf("[");
  for (idx=0, backmask = 1UL; idx<mti->patternlength; idx++, backmask <<= 1)
  {
    if (col->prefixofsuffix & backmask)
    {
      if (first)
      {
        printf("%lu",idx);
        first = false;
      } else
      {
        printf(",%lu",idx);
      }
    }
  }
  printf("\n");
}

#endif

static void *pms_allocatedfsconstinfo(unsigned int alphasize)
{
  Matchtaskinfo *mti = ma_malloc(sizeof(Matchtaskinfo));
  mti->eqsvector = ma_malloc(sizeof(*mti->eqsvector) * alphasize);
  return mti;
}

static void pms_initdfsconstinfo(void *dfsconstinfo,
                                 unsigned int alphasize,
                                 const Uchar *pattern,
                                 unsigned long patternlength,
                                 UNUSED unsigned long maxdistance,
                                 UNUSED Seqpos maxintervalwidth)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;
#ifdef SKDEBUG
  int a;
#endif

  initeqsvector(mti->eqsvector,(unsigned long) alphasize,pattern,patternlength);
#ifdef SKDEBUG
  for (a=0; a<4; a++)
  {
    char buffer[32+1];
    uint32_t2string(buffer,mti->eqsvector[a]);
    printf("# %d->%s\n",a,buffer);
  }
#endif
  mti->patternlength = patternlength;
}

static void pms_extractdfsconstinfo(UNUSED void *dfsconstinfo)
{
  unsigned long idx;
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

  for (idx=0; idx<mti->patternlength; idx++)
  {
    if (mti->mstatlength[idx] > 0)
    {
      printf("%lu->%lu ",idx,mti->mstatlength[idx]);
    }
  }
  printf("\n");
}

static void pms_freedfsconstinfo(void **dfsconstinfo)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) *dfsconstinfo;

  ma_free(mti->eqsvector);
  ma_free(mti);
  *dfsconstinfo = NULL;
}

static unsigned long zerosontheright(unsigned long v)
{
  unsigned long c;     /* c will be the number of zero bits on the right,
                         so if v is 1101000 (base 2), then c will be 3 */
  assert(v > 0);
  assert(sizeof (unsigned long) == (size_t) 4);
  if (v & 0x1)
  {
    c = 0; /* special case for odd v (assumed to happen half of the time) */
  } else
  {
    c = 1UL;
    if ((v & 0xffff) == 0)
    {
      v >>= 16;
      c += 16UL;
    }
    if ((v & 0xff) == 0)
    {
      v >>= 8;
      c += 8UL;
    }
    if ((v & 0xf) == 0)
    {
      v >>= 4;
      c += 4UL;
    }
    if ((v & 0x3) == 0)
    {
      v >>= 2;
      c += 2UL;
    }
    c -= v & 0x1;
  }
  return c;
}

static void pms_initParallelmstats(DECLAREPTRDFSSTATE(aliascolumn),
                                   void *dfsconstinfo)
{
  Parallelmstats *column = (Parallelmstats *) aliascolumn;
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;
  unsigned long idx;

  column->prefixofsuffix = ~0UL;
  assert(mti->patternlength <= (unsigned long) INTWORDSIZE);
  for (idx = 0; idx<mti->patternlength; idx++)
  {
    mti->mstatlength[idx] = 0;
  }
}

static unsigned long pms_nextstepfullmatches(
                              DECLAREPTRDFSSTATE(aliascolumn),
                              UNUSED Seqpos width,
                              unsigned long currentdepth,
                              void *dfsconstinfo)
{
  Parallelmstats *limdfsstate = (Parallelmstats *) aliascolumn;

  if (limdfsstate->prefixofsuffix == 0)
  {
    return 0; /* stop depth first traversal */
  } else
  {
    Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;
    unsigned long tmp, bitindex, first1;
    tmp = limdfsstate->prefixofsuffix;
    bitindex = 0;
    do
    {
      first1 = zerosontheright(tmp);
      assert(bitindex + first1 < mti->patternlength);
      if (mti->mstatlength[bitindex+first1] < currentdepth)
      {
        mti->mstatlength[bitindex+first1] = currentdepth;
      }
      tmp >>= (first1+1);
      bitindex += (first1+1);
    } while (tmp != 0);
  }
  return 1UL; /* continue with depth first traversal */
}

static void pms_nextParallelmstats(const void *dfsconstinfo,
                                   DECLAREPTRDFSSTATE(aliasoutcol),
                                   unsigned long currentdepth,
                                   Uchar currentchar,
                                   const DECLAREPTRDFSSTATE(aliasincol))
{
#ifdef SKDEBUG
  char buffer1[32+1], buffer2[32+1];
#endif
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Parallelmstats *outcol = (Parallelmstats *) aliasoutcol;
  const Parallelmstats *incol = (const Parallelmstats *) aliasincol;

  assert(ISNOTSPECIAL(currentchar));
  assert(currentdepth > 0);

  outcol->prefixofsuffix = incol->prefixofsuffix &
                           (mti->eqsvector[currentchar] << (currentdepth-1));
#ifdef SKDEBUG
  uint32_t2string(buffer1,(uint32_t) incol->prefixofsuffix);
  uint32_t2string(buffer2,(uint32_t) outcol->prefixofsuffix);
  printf("next(%s,%u,depth=%lu)->%s\n",buffer1,currentchar,
                                       currentdepth,buffer2);
#endif
}

static void pms_inplacenextParallelmstats(const void *dfsconstinfo,
                                          DECLAREPTRDFSSTATE(aliascol),
                                          unsigned long currentdepth,
                                          Uchar currentchar)
{
#ifdef SKDEBUG
  char buffer1[32+1], buffer2[32+1];
#endif
  unsigned long tmp;
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Parallelmstats *col = (Parallelmstats *) aliascol;

  assert(ISNOTSPECIAL(currentchar));
  tmp = col->prefixofsuffix;
  col->prefixofsuffix &= (mti->eqsvector[currentchar] << (currentdepth-1));
#ifdef DEBUG
  uint32_t2string(buffer1,(uint32_t) tmp);
  uint32_t2string(buffer2,(uint32_t) col->prefixofsuffix);
  printf("inplacenext(%s,%u,%lu)->%s\n",buffer1,currentchar,currentdepth,
                                        buffer2);
#endif
}

const AbstractDfstransformer *pms_AbstractDfstransformer(void)
{
  static const AbstractDfstransformer pms_adfst =
  {
    sizeof (Parallelmstats),
    pms_allocatedfsconstinfo,
    pms_initdfsconstinfo,
    pms_extractdfsconstinfo,
    pms_freedfsconstinfo,
    pms_initParallelmstats,
    pms_nextstepfullmatches,
    pms_nextParallelmstats,
    pms_inplacenextParallelmstats,
#ifdef SKDEBUG
    pms_showParallelmstats,
#endif
  };
  return &pms_adfst;
}

/*
  define bitvector prefixofsuffix_{d} such that after processing a sequence v
  of length d we have: for all i\in[0,m-1]
  prefixofsuffix_{d}[i] is 1 iff P[i..i+d-1] = v[0..d-1]

  Let eqsvector_{a} be a vector of size m such that
  eqsvector_{a}[i]=1 if P[i]=a

  Let d=0 (i.e. at the root). Then
  P[i..i+d-1]=P[i..i-1]=\varepsilon=v[0..-1]=v[0..d-1] for all i \in[0..m-1]
  and hence prefixofsuffix_{d}[i]=1. In other words
  prefixofsuffix_{d} = 1^{m}.

  Now suppose d > 0 and assume we have computed
  prefixofsuffix_{d-1}. Then by definition
  prefixofsuffix_{d}[i]
    iff P[i..i+d-1] = v[0..d-1]
    iff P[i..i+d-2] = v[0..d-2] && P[i+d-1]=v[d-1]
    iff prefixofsuffix_{d-1][i]=1 && eqsvector_{v[d-1]}[i+d-1]=1
    iff prefixofsuffix_{d-1][i] & eqsvector_{v[d-1]}[i+d-1]

  All values in prefixofsuffix_{d} are independent and can be computed
  in parallel by

  prefixofsuffix_{d} = prefixofsuffix_{d-1} & (eqsvector_{a} << (d-1))
  where a=v[d-1]

  prefixofsuffix_{d] = 0 and
  prefixofsuffix_{d-1} != 0 then  for all i satisfying
  prefixofsuffix_{d-1][i] = 1 do:
    if mstats[i]<d then mstats[i]=d and store first suffixposition of current
    interval.
*/
