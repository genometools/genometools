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
  define bitvector prefixofsuffix such that after processing a sequence v 
  of length d we have: for all i\in[0,m]
  prefixofsuffix[i] is 1 iff P[i..i+d-1] = v[0..d-1]
*/
