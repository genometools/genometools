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

#ifndef TAGERATOR_H
#define TAGERATOR_H
#include <stdbool.h>
#include "core/strarray.h"
#include "core/error.h"

typedef struct
{
  GtStrArray *tagfiles;
  GtStr *esaindexname,
      *pckindexname;
  bool online,  /* perform online search, for testing */
       docompare, /* compare results with online search */
       replacewildcard, /* replace wildcards by random symbol */
       nofwdmatch, /* do not perform matching on forward strand */
       norcmatch, /* do not perform matching on reverse complemented strand */
       nowildcards, /* ignore matches containing wildcards */
       skpp; /* Skip prefix of pattern without counting errors */
  long maxdistance; /* maximal number of allowed differences >= 0 */
  int userdefinedmaxdepth;   /* use pckbuckets only up to this depth */
  unsigned long maxintervalwidth; /* max width of interval */
} TageratorOptions;

int runtagerator(const TageratorOptions *tageratoroptions,GT_Error *err);

#endif
