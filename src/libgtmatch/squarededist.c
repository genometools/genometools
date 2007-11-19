#include "libgtcore/env.h"
#include "libgtcore/symboldef.h"
#include "spacedef.h"

static unsigned long squarededistunit2 (const Uchar *u, unsigned long m,
                                        const Uchar *v, unsigned long n,
                                        Env *env)
{
  unsigned long val, we, nw, *ecol, *ecolptr;
  const Uchar *uptr, *vptr;

  ALLOCASSIGNSPACE(ecol,NULL,unsigned long,m+1);
  for (*ecol = 0, ecolptr = ecol+1, uptr = u; uptr < u + m; ecolptr++, uptr++)
  {
    *ecolptr = *(ecolptr-1) + 1;
  }
  for (vptr = v; vptr < v + n; vptr++)
  {
    nw = *ecol;
    *ecol = nw + 1;
    for (ecolptr = ecol+1, uptr = u; uptr < u + m; ecolptr++, uptr++)
    {
      we = *ecolptr;
      *ecolptr = *(ecolptr-1) + 1;
      if (*uptr == *vptr)
      {
        val = nw;
      } else
      {
        val = nw + 1;
      }
      if (val < *ecolptr)
      {
        *ecolptr = val;
      }
      if ((val = we + 1) < *ecolptr)
      {
        *ecolptr = val;
      }
      nw = we;
    }
  }
  val = *(ecolptr-1);
  FREESPACE(ecol);
  return val;
}

unsigned long squarededistunit (const Uchar *u, unsigned long m,
                                const Uchar *v, unsigned long n,
                                Env *env)
{
  if (m < n)
  {
    return squarededistunit2(u,m,v,n,env);
  }
  return squarededistunit2(v,n,u,m,env);
}
