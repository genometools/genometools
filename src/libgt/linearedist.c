/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "linearedist.h"
#include "minmax.h"
#include "xansi.h"

static void fillDPtable(unsigned long *dptable,
                        const char *u, unsigned long n,
                        const char *v, unsigned long m)
{
  unsigned long i, j , nw, we;
  for (i = 0; i <= n; i++)
    dptable[i] = i;
  for (j = 1; j <= m; j++) {
    nw = dptable[0];
    dptable[0] = j;
    for (i = 1; i <= n; i++) {
      we = dptable[i];
      dptable[i] = nw + (u[i-1] == v[j-1] ? 0 : 1); /* replacement */
      if (dptable[i-1] + 1 < dptable[i]) /* deletion */
        dptable[i] = dptable[i-1] + 1;
      if (we + 1 < dptable[i]) /* insertion */
        dptable[i] = we + 1;
      nw = we;
    }
  }
}

unsigned long linearedist(const char *u, unsigned long n,
                          const char *v, unsigned long m, Env *env)
{
  unsigned long *dptable, edist;
  dptable = env_ma_malloc(env, sizeof (unsigned long) * (MIN(n,m) + 1));
  fillDPtable(dptable, n <= m ? u : v, MIN(n,m), n <= m ? v : u, MAX(n,m));
  edist = dptable[MIN(n,m)];
  env_ma_free(dptable, env);
  return edist;
}
