#include "libgtcore/symboldef.h"
#include "libgtcore/ma.h"

typedef struct
{
  unsigned long *endindex,
                *positions;
} Charatpos;

Charatpos *newCharatpos(unsigned long patternlength,unsigned int alphasize)
{
  Charatpos *catpos;

  catpos = ma_malloc(sizeof(*catpos));
  catpos->endindex = ma_malloc(sizeof(unsigned long) * alphasize);
  catpos->positions = ma_malloc(sizeof(unsigned long) * patternlength);
  return catpos;
}

Charatpos *reinitCharatpos(Charatpos *catpos,
                           const Uchar *pattern,unsigned long patternlength,
                           unsigned int alphasize)
{
  const Uchar *pptr;
  unsigned long partialsum, tmp;
  unsigned int idx;

  for (idx=0; idx<alphasize; idx++)
  {
    catpos->endindex[idx] = 0;
  }
  for (pptr=pattern; pptr<pattern+ patternlength; pptr++)
  {
    catpos->endindex[(int) *pptr]++;
  }
  partialsum = catpos->endindex[0];
  catpos->endindex[0] = 0;
  for (idx=1U; idx<alphasize; idx++)
  {
    tmp = catpos->endindex[idx];
    catpos->endindex[idx] = partialsum;
    partialsum += tmp;
  }
  for (pptr=pattern; pptr<pattern+patternlength; pptr++)
  {
    catpos->positions[catpos->endindex[*pptr]++]
      = (unsigned long) (pptr - pattern);
  }
  return catpos;
}

void wrapCharatpos(Charatpos **catposptr)
{
  Charatpos *catpos = *catposptr;
  ma_free(catpos->endindex);
  ma_free(catpos->positions);
  ma_free(catpos);
  *catposptr = NULL;
}

void maintainnullcols(const Charatpos *catpos,
                      unsigned long *front0,Uchar cc,unsigned long depth)
{
  unsigned long idx;

  for (idx = (cc == 0) ? 0 : catpos->endindex[cc-1];
       idx < catpos->endindex[cc]; idx++)
  {
    unsigned long pos = catpos->positions[idx];
    if (front0[pos] == depth)
    {
      front0[pos]++;
    }
  }
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
