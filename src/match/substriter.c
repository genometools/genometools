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

#include <inttypes.h>
#include "core/str_array.h"
#include "core/seq_iterator_api.h"
#include "core/chardef.h"
#include "core/codetype.h"
#include "core/ma_api.h"
#include "substriter.h"
#include "qgram2code.h"
#include "initbasepower.h"

struct Substriter
{
  const GtUchar *currentptr, *start;
  unsigned long remaining;
  GtCodetype currentcode;
  unsigned int qvalue,
               numofchars;
  GtCodetype **multimappower;
};

Substriter *gt_substriter_new(const GtAlphabet *alphabet,unsigned int qvalue)
{
  Substriter *substriter;

  substriter = gt_malloc(sizeof *substriter);
  substriter->qvalue = qvalue;
  substriter->numofchars = gt_alphabet_num_of_chars(alphabet);
  substriter->multimappower = gt_initmultimappower(substriter->numofchars,
                                                   qvalue);
  return substriter;
}

void gt_substriter_init(Substriter *substriter,const GtUchar *start,
                    unsigned long len)
{
  substriter->start = substriter->currentptr = start;
  substriter->remaining = len;
  gt_assert(substriter->remaining > 0);
}

int gt_substriter_next(Substriter *substriter)
{
  unsigned int firstspecial;

  while (true)
  {
    if (substriter->remaining >= (unsigned long) substriter->qvalue)
    {
      firstspecial = qgram2code(&substriter->currentcode,
                                (const GtCodetype **) substriter->multimappower,
                                substriter->qvalue,
                                substriter->currentptr);
      if (firstspecial == substriter->qvalue)
      {
        substriter->remaining--;
        substriter->currentptr++;
        break;
      }
      substriter->remaining -= firstspecial + 1;
      substriter->currentptr += firstspecial + 1;
    } else
    {
      return 0;
    }
  }
  return 1;
}

void gt_substriter_delete(Substriter *substriter)
{
  if (substriter != NULL)
  {
    gt_multimappower_delete(substriter->multimappower);
    gt_free(substriter);
  }
}
