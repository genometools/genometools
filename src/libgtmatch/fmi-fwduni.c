#include <math.h>
#include "chardef.h"
#include "divmodmul.h"
#include "fmindex.h"

#include "fmi-occ.gen"

Seqpos skfmuniqueforward (const Fmindex *fm,
                          const Uchar *qstart,
                          const Uchar *qend)
{
  Uchar cc;
  const Uchar *qptr;
  Bwtbound bwtbound;

  qptr = qstart;
  cc = *qptr++;
  if(ISSPECIAL(cc))
  {
    return 0;
  }
  bwtbound.lbound = fm->tfreq[cc];
  bwtbound.ubound = fm->tfreq[cc+1];
  while (qptr < qend && bwtbound.lbound + 1 < bwtbound.ubound)
  {
    cc = *qptr;
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    bwtbound.lbound = fm->tfreq[cc] + skfmocc (fm, cc, bwtbound.lbound);
    bwtbound.ubound = fm->tfreq[cc] + skfmocc (fm, cc, bwtbound.ubound);
    qptr++;
  }
  if(bwtbound.lbound + 1 == bwtbound.ubound)
  {
    return (Seqpos) (qptr - qstart);
  }
  return 0;
}
