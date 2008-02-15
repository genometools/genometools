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

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include "intcode-def.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "bckbound.h"

typedef struct
{
  const Seqpos *bcktab,
               *countspecialcodes;
  Codetype numofallcodes;
  unsigned int prefixlength;
  Bucketspecification bucketspec;
} Bucketenumerator;

Bucketenumerator *newbucketenumerator(const Seqpos *bcktab,
                                      const Seqpos *countspecialcodes,
                                      Codetype numofallcodes,
                                      unsigned int prefixlength)
{
  Bucketenumerator *bucketenumerator;

  ALLOCASSIGNSPACE(bucketenumerator,NULL,Bucketenumerator,1);
  bucketenumerator->bcktab = bcktab;
  bucketenumerator->countspecialcodes = countspecialcodes;
  bucketenumerator->numofallcodes = numofallcodes;
  bucketenumerator->prefixlength = prefixlength;
  return bucketenumerator;
}
