/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>
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

#ifndef TYPES_API_H
#define TYPES_API_H

/* Define the conversion string for '%lld' in platform independent fashion. */
#ifndef _WIN32
#define GT_LLD "%lld"
#else
#define GT_LLD "%I64d"
#endif

/* Define the conversion string for '%llu' in platform independent fashion. */
#ifndef _WIN32
#define GT_LLU "%llu"
#else
#define GT_LLU "%I64u"
#endif

/* Define the conversion string for 'ld' in platform independent fashion. */
#ifndef _WIN64
#define GT_LDS "ld"
#else
#define GT_LDS "I64d"
#endif

/* Define the conversion string for 'ld' in platform independent fashion. */
#define GT_LD "%"GT_LDS

/* Define the conversion string for 'lu' in platform independent fashion. */
#ifndef _WIN64
#define GT_LUS "lu"
#else
#define GT_LUS "I64u"
#endif

/* Define the conversion string for '%lu' in platform independent fashion. */
#define GT_LU "%"GT_LUS

/* Define the conversion string for '%zu' in platform independent fashion. */
#if !defined(_WIN32)
#define GT_ZU "%zu"
#elif defined(_WIN64)
#define GT_ZU "%I64u"
#else
#define GT_ZU "%u"
#endif

/* Define GtUword as an unsigned integer with the machine word size (4 byte on
   32-bit systems and 8 byte on 64-bit systems). */
#ifdef _WIN64
typedef unsigned long long GtUword;
#else
typedef unsigned long GtUword;
#endif

/* Define GtWord as a signed integer with the machine word size (4 byte on
   32-bit systems and 8 byte on 64-bit systems). */
#ifdef _WIN64
typedef long long GtWord;
#else
typedef long GtWord;
#endif

typedef long long          GtInt64;
typedef unsigned long long GtUint64;

typedef unsigned char GtUchar;

/* deprecated */
typedef GtUword GtUlong;

/* deprecated */
typedef struct {
  GtUlong a, b;
} GtUlongPair;

#endif
