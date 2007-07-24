#include <math.h>
#include "chardef.h"
#include "divmodmul.h"
#include "fmindex.h"

#include "fmi-occ.gen"

unsigned long skfmuniqueforward (const Fmindex *fmindex,
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
  bwtbound.lbound = fmindex->tfreq[cc];
  bwtbound.ubound = fmindex->tfreq[cc+1];
  while (qptr < qend && bwtbound.lbound + 1 < bwtbound.ubound)
  {
    cc = *qptr;
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    bwtbound.lbound = fmindex->tfreq[cc] + 
                      fmoccurrence (fmindex, cc, bwtbound.lbound);
    bwtbound.ubound = fmindex->tfreq[cc] + 
                      fmoccurrence (fmindex, cc, bwtbound.ubound);
    qptr++;
  }
  if(bwtbound.lbound + 1 == bwtbound.ubound)
  {
    return (unsigned long) (qptr - qstart);
  }
  return 0;
}
