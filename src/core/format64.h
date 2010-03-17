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

#ifndef FORMAT64_H
#define FORMAT64_H
#include <inttypes.h>

/*
  This file contains some basic type definitions for 64 bit integer types.
*/

#ifdef S_SPLINT_S
#define Formatuint64_t "%lu"
#define SCANint64_tcast(X) ((long *) (X))
#define PRINTuint64_tcast(X) ((unsigned long) (X))
#else
#define Formatuint64_t "%" PRIu64
#define SCANint64_tcast(X) (X)
#define PRINTuint64_tcast(X) (X)
#endif

#ifdef S_SPLINT_S
#define FormatScanint64_t "%ld"
#else
#define FormatScanint64_t "%" SCNd64
#endif

#endif
