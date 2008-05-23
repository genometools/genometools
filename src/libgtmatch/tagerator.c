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

#include "libgtcore/unused.h"
#include "libgtcore/strarray.h"
#include "libgtcore/ma.h"
#include "libgtcore/error.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/seqiterator.h"
#include "libgtmatch/tagerator.h"
#include "libgtmatch/sarr-def.h"

#include "libgtmatch/esa-map.pr"

int runtagerator(const TageratorOptions *tageratoroptions,Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  int retval;
  unsigned long len;
  unsigned int demand = SARR_SUFTAB | SARR_ESQTAB;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     demand,
                     tageratoroptions->indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  seqit = seqiterator_new(tageratoroptions->tagfiles, NULL, true);
  while (true)
  {
    retval = seqiterator_next(seqit, &sequence, &len, &desc, err);
    if (retval != 1)
    {
      break;
    }
    ma_free(desc);
  }
  seqiterator_delete(seqit);
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
