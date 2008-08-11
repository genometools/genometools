/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "libgtcore/str.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/eis-voiditf.h"
#include "libgtmatch/esa-map.pr"

typedef struct
{
  Str *indexname;
  unsigned long maxdepth;
} Prebwtoptions;

int runprebwt(const Prebwtoptions *prebwtoptions,Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  void *packedindex = NULL;
  bool haserr = false;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     0,
                     prebwtoptions->indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    packedindex = loadvoidBWTSeqForSA(prebwtoptions->indexname,
                                      &suffixarray,
                                      totallength, err);
    if (packedindex == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    Matchbound *boundsarray;
    unsigned int alphasize = getnumofcharsAlphabet(suffixarray.alpha);
    unsigned long numofbounds;

    numofbounds = (unsigned long) pow((double) alphasize,
                                      (double) prebwtoptions->maxdepth);
    boundsarray = ma_malloc(sizeof(*boundsarray) * numofbounds);
    pck_precomputebounds(boundsarray,
                         numofbounds,
                         (const void *) packedindex,
                         alphasize,
                         totallength,
                         prebwtoptions->maxdepth);
  }
  freesuffixarray(&suffixarray);
  if (packedindex != NULL)
  {
    deletevoidBWTSeq(packedindex);
  }
  return haserr ? -1 : 0;
}
