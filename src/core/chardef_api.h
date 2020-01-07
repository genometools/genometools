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

#ifndef CHARDEF_API_H
#define CHARDEF_API_H

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

#define GT_SEPARATOR       UCHAR_MAX

/*
  wildcard symbol in multiple seq
*/

#define GT_WILDCARD        (GT_SEPARATOR-1)

/*
  undefined character, only to be used in conjunction with symbol maps
*/

#define GT_UNDEFCHAR       (GT_SEPARATOR-2)

/*
  either GT_WILDCARD or GT_SEPARATOR
*/

#define GT_ISSPECIAL(C)    ((C) >= (GtUchar) GT_WILDCARD)

/*
  neither GT_WILDCARD nor GT_SEPARATOR
*/

#define GT_ISNOTSPECIAL(C) ((C) < (GtUchar) GT_WILDCARD)

/*
  undefined character, only to be used in conjunction with the Burrows-Wheeler
  transform
*/

#define GT_UNDEFBWTCHAR    GT_WILDCARD

/*
  Either special character or GT_UNDEFBWTCHAR
*/

#define GT_ISBWTSPECIAL(C) ((C) >= (GtUchar) GT_UNDEFBWTCHAR)

/*@unused@*/ static inline GtUword
                            gt_containsspecialbytestring(const GtUchar *seq,
                                                      GtUword offset,
                                                      GtUword len)
{
  const GtUchar *sptr;

  gt_assert(offset < len);
  for (sptr=seq+offset; sptr < seq + len; sptr++)
  {
    if (GT_ISSPECIAL(*sptr))
    {
      return (GtUword) (sptr - seq);
    }
  }
  return len;
}

typedef struct
{
  GtUword specialcharacters,      /* total number of special syms */
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
