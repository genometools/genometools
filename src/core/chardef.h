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

#ifndef CHARDEF_H
#define CHARDEF_H

#include "core/types_api.h"
#include "core/assert_api.h"
#include <limits.h>

/*
  This file defines some character values used when storing
  multiple sequences.
*/

/*
  separator symbol in multiple seq
*/

#define SEPARATOR       UCHAR_MAX

/*
  wildcard symbol in multiple seq
*/

#define WILDCARD        (SEPARATOR-1)

/*
  undefined character, only to be used in conjunction with symbol maps
*/

#define UNDEFCHAR       (SEPARATOR-2)

/*
  either WILDCARD or SEPARATOR
*/

#define ISSPECIAL(C)    ((C) >= (GtUchar) WILDCARD)

/*
  neither WILDCARD nor SEPARATOR
*/

#define ISNOTSPECIAL(C) ((C) < (GtUchar) WILDCARD)

/*
  undefined character, only to be used in conjunction with the Burrows-Wheeler
  transform
*/

#define UNDEFBWTCHAR    WILDCARD

/*
  Either special character or UNDEFBWTCHAR
*/

#define ISBWTSPECIAL(C) ((C) >= (GtUchar) UNDEFBWTCHAR)

/*@unused@*/ static inline unsigned long
                            containsspecialbytestring(const GtUchar *seq,
                                                      unsigned long offset,
                                                      unsigned long len)
{
  const GtUchar *sptr;

  gt_assert(offset < len);
  for (sptr=seq+offset; sptr < seq + len; sptr++)
  {
    if (ISSPECIAL(*sptr))
    {
      return (unsigned long) (sptr - seq);
    }
  }
  return len;
}

typedef struct
{
  unsigned long specialcharacters,      /* total number of special syms */
                specialranges,          /* number of stored ranges (of maximal
                                           length according to the chosen
                                           representation) with special syms */
                realspecialranges,      /* number of ranges with special syms */
                lengthofspecialprefix,  /* number of specials at start of
                                           sequence */
                lengthofspecialsuffix,  /* number of specials at end of
                                           sequence */
                wildcards,              /* total number of wildcards */
                wildcardranges,         /* number of stored ranges (of maximal
                                           length according to the chosen
                                           representation) with wildcards */
                realwildcardranges,     /* number of ranges with wildcards */
                lengthofwildcardprefix, /* number of wildcards at start of
                                           sequence */
                lengthofwildcardsuffix, /* number of wildcards at end of
                                           sequence */
                lengthoflongestnonspecial, /* length of longest non-special
                                              stretch */
                exceptioncharacters,    /* total number of exception chars */
                exceptionranges,        /* number of stored exception ranges */
                realexceptionranges;    /* number of exception ranges */
} GtSpecialcharinfo;

#endif
