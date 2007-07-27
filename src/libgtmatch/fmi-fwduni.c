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
#undef mydebug
#ifdef mydebug
  printf("# start cc=%u\n",cc);
#endif
  if(ISSPECIAL(cc))
  {
    return 0;
  }
  bwtbound.lbound = fmindex->tfreq[cc];
  bwtbound.ubound = fmindex->tfreq[cc+1];
#ifdef mydebug
  printf("# bounds=%u,%u = %u occurrences\n",
              (unsigned int) bwtbound.lbound,
              (unsigned int) bwtbound.ubound,
              bwtbound.ubound - bwtbound.lbound);
#endif
  while (qptr < qend && bwtbound.lbound + 1 < bwtbound.ubound)
  {
    cc = *qptr;
#ifdef mydebug
    printf("# cc=%u\n",cc);
#endif
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    bwtbound.lbound = fmindex->tfreq[cc] + 
                      fmoccurrence (fmindex, cc, bwtbound.lbound);
    bwtbound.ubound = fmindex->tfreq[cc] + 
                      fmoccurrence (fmindex, cc, bwtbound.ubound);
#ifdef mydebug
    printf("# bounds=%u,%u = %u occurrences\n",
                (unsigned int) bwtbound.lbound,
                (unsigned int) bwtbound.ubound,
                bwtbound.ubound - bwtbound.lbound);
#endif
    qptr++;
  }
  if(bwtbound.lbound + 1 == bwtbound.ubound)
  {
    return (unsigned long) (qptr - qstart);
  }
  return 0;
}
