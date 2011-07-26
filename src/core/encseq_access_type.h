/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCSEQ_ACCESS_TYPE_H
#define ENCSEQ_ACCESS_TYPE_H

#include "core/defined-types.h"
#include "core/error_api.h"

typedef enum
{
  GT_ACCESS_TYPE_DIRECTACCESS,
  GT_ACCESS_TYPE_BYTECOMPRESS,
  GT_ACCESS_TYPE_EQUALLENGTH,
  GT_ACCESS_TYPE_BITACCESS,
  GT_ACCESS_TYPE_UCHARTABLES,
  GT_ACCESS_TYPE_USHORTTABLES,
  GT_ACCESS_TYPE_UINT32TABLES,
  GT_ACCESS_TYPE_UNDEFINED
} GtEncseqAccessType;

GtEncseqAccessType gt_encseq_access_type_get(const char *str);
bool               gt_encseq_access_type_isviautables(GtEncseqAccessType sat);
const char*        gt_encseq_access_type_str(GtEncseqAccessType at);
const char*        gt_encseq_access_type_list(void);
int                gt_encseq_access_type_determine(unsigned long *specialranges,
                                         unsigned long *wildcardranges,
                                         unsigned long totallength,
                                         unsigned long numofsequences,
                                         unsigned long numofdbfiles,
                                         unsigned long lengthofalphadef,
                                         unsigned long lengthofdbfilenames,
                                         const unsigned long *specialrangestab,
                                         const unsigned long *wildcardrangestab,
                                         const Definedunsignedlong *equallength,
                                         unsigned int numofchars,
                                         const char *str_sat,
                                         GtError *err);
#endif
