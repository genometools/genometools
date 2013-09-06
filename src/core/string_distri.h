/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STRING_DISTRI_H
#define STRING_DISTRI_H

#include "core/error.h"
#include "core/file.h"

/* A discrete distribution */
typedef struct GtStringDistri GtStringDistri;

typedef void (*GtStringDistriIterFunc)(const char *string,
                                       GtUword occurrences,
                                       double probability, void *data);

GtStringDistri* gt_string_distri_new(void);
void            gt_string_distri_add(GtStringDistri*, const char*);
/* <string_distri> must contain at least one element with given <key>. */
void            gt_string_distri_sub(GtStringDistri *string_distri,
                                     const char *key);
GtUword   gt_string_distri_get(const GtStringDistri*, const char*);
/* return probability */
double          gt_string_distri_get_prob(const GtStringDistri*, const char*);
void            gt_string_distri_foreach(const GtStringDistri*,
                                         GtStringDistriIterFunc, void *data);
void            gt_string_distri_delete(GtStringDistri*);

#endif
