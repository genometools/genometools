/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FORMAT64_H
#define FORMAT64_H
#include <inttypes.h>

/*
  This file contains some basic type definitions for 64 bit integer types.
*/

#ifdef S_SPLINT_S
#define Formatuint64_t "%lu"
#define Scanuint64_tcast(X) ((unsigned long *) (X))
#define PRINTuint64_tcast(X) ((unsigned long) (X))
#else
#define Formatuint64_t "%" PRIu64
#define Scanuint64_tcast(X) (X)
#define PRINTuint64_tcast(X) (X)
#endif

#ifdef S_SPLINT_S
#define FormatScanint64_t "%lu"
#else
#define FormatScanint64_t "%" SCNd64
#endif

#endif
