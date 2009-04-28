/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#include "core/ma.h"
#include "core/assert_api.h"
#include "divmodmul.h"
#include "diff-cover.h"

struct Differencecover
{
  unsigned int vparam, logmod, size, *coverrank;
  Diffvalue *diffvalues;
};

static Diffvalue differencecovertab[] =
{
  0UL,
  0UL,1UL,
  0UL,1UL,2UL,
  0UL,1UL,2UL,4UL,
  0UL,1UL,2UL,5UL,8UL,
  0UL,1UL,2UL,3UL,7UL,11UL,19UL,
  0UL,1UL,2UL,5UL,14UL,16UL,34UL,42UL,59UL,
  0UL,1UL,3UL,7UL,17UL,40UL,55UL,64UL,75UL,85UL,104UL,109UL,117UL
};

static unsigned int differencecoversizes[] = {1U,2U,3U,4U,5U,7U,9U,13U};

static void fillcoverrank(Differencecover *dcov)
{
  unsigned int i, j;

  dcov->coverrank = gt_malloc(sizeof(*dcov->coverrank) * dcov->vparam);
  for (i=0, j=0; i<dcov->vparam; i++)
  {
    dcov->coverrank[i] = j;
    if (j<dcov->size && dcov->diffvalues[j] <= (Diffvalue) i)
    {
      j++;
    }
  }
}

Differencecover *differencecover_new(unsigned int vparam)
{
  size_t logmod;
  unsigned int offset = 0, v = 1U;
  Differencecover *dcov;
  bool found = false;

  dcov = gt_malloc(sizeof (*dcov));
  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    if (v == vparam)
    {
      dcov->size = differencecoversizes[logmod];
      dcov->diffvalues = differencecovertab + offset;
      found = true;
      break;
    }
    offset += differencecoversizes[logmod];
    v = MULT2(v);
  }
  if (!found)
  {
    gt_free(dcov);
    return NULL;
  }
  dcov->logmod = (unsigned int) logmod;
  dcov->vparam = 1U << logmod;
  fillcoverrank(dcov);
  return dcov;
}

void differencecover_delete(Differencecover *dcov)
{
  gt_free(dcov->coverrank);
  gt_free(dcov);
}

void differencecovers_check(void)
{
  Differencecover *dcov;
  size_t logmod, next = 0;
  unsigned int j, vparam;

  for (logmod = 0;
       logmod < sizeof (differencecoversizes)/sizeof (differencecoversizes[0]);
       logmod++)
  {
    vparam = 1U << logmod;
    dcov = differencecover_new(vparam);
    if (dcov == NULL)
    {
      fprintf(stderr,"no difference cover for v=%u\n",vparam);
      exit(EXIT_FAILURE);
    }
    for (j = 0; j<dcov->size; j++)
    {
      gt_assert(dcov->diffvalues[j] == differencecovertab[next]);
      next++;
    }
    differencecover_delete(dcov);
  }
  gt_assert(next == sizeof (differencecovertab)/sizeof (differencecovertab[0]));
  printf("# %u difference covers checked\n",(unsigned int) logmod);
}
