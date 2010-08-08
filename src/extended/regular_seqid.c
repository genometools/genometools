/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/types_api.h"
#include "extended/regular_seqid.h"

void gt_regular_seqid_save(GtStr *seqid, GtStr *description)
{
  unsigned long i, len;
  unsigned char *desc, cc;

  gt_assert(seqid && description);

  len  = gt_str_length(description);
  desc = (GtUchar*) gt_str_get(description);

  i = 0;

  if ((len >= 2) && (desc[0] == 'g') && (desc[1] == 'i') && (desc[2] == '|')) {
    /* skip 'gi|' */
    i = 3;
  }
  else if ((len >= 2) && (desc[0] == 'S') && (desc[1] == 'Q') &&
           (desc[2] == ';')) {
    /* skip 'SQ;' */
    i = 3;
  }
  else if ((len >= 3) && (desc[0] == '(') && (desc[1] == 'g') &&
           (desc[2] == 'i') && (desc[3] == '|')) {
    /* skip '(gi|' */
    i = 4;
  }
  else if ((len >= 3) && (desc[0] == 'r') && (desc[1] == 'e') &&
           (desc[2] == 'f') && (desc[3] == '|')) {
    /* skip 'ref|' */
    i = 4;
  }

  for (/* init already done */ ; i < len; i++) {
    cc = desc[i];
    if (cc == ':' || cc == '|' || cc == '\t' || cc == ' ')
      break;
    gt_str_append_char(seqid, desc[i]);
  }
}
