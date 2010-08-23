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

#include <string.h>
#include "core/assert_api.h"
#include "core/encseq_access_type.h"

typedef struct
{
  GtEncseqAccessType sat;
  char *name;
} GtWrittenAccessType;

static GtWrittenAccessType wpa[] = {
  {GT_ACCESS_TYPE_DIRECTACCESS, "direct"},
  {GT_ACCESS_TYPE_BYTECOMPRESS, "bytecompress"},
  {GT_ACCESS_TYPE_BITACCESS, "bit"},
  {GT_ACCESS_TYPE_UCHARTABLES, "uchar"},
  {GT_ACCESS_TYPE_USHORTTABLES, "ushort"},
  {GT_ACCESS_TYPE_UINT32TABLES, "uint32"}
};

const char* gt_encseq_access_type_list(void)
{
  return "direct, bytecompress, bit, uchar, ushort, uint32";
}

const char* gt_encseq_access_type_str(GtEncseqAccessType at)
{
  gt_assert((int) at < (int) GT_ACCESS_TYPE_UNDEFINED);
  return wpa[at].name;
}

GtEncseqAccessType gt_encseq_access_type_get(const char *str)
{
  size_t i;

  for (i=0; i<sizeof (wpa)/sizeof (wpa[0]); i++)
  {
    if (strcmp(str,wpa[i].name) == 0)
    {
      return wpa[i].sat;
    }
  }
  return GT_ACCESS_TYPE_UNDEFINED;
}
