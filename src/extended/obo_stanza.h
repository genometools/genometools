/*
  Copyright (c) 2008-2009, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008            Center for Bioinformatics, University of Hamburg

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

#ifndef OBO_STANZA_H
#define OBO_STANZA_H

#include "core/str_api.h"

typedef struct GtOBOStanza GtOBOStanza;

GtOBOStanza*  gt_obo_stanza_new(const char *type, unsigned long line,
                                GtStr *filename);
void          gt_obo_stanza_delete(GtOBOStanza *obo_stanza);
void          gt_obo_stanza_add(GtOBOStanza *obo_stanza,
                                const char *tag, const char *value);
const char*   gt_obo_stanza_get_type(const GtOBOStanza *obo_stanza);
const char*   gt_obo_stanza_get_value(const GtOBOStanza *obo_stanza,
                                      const char *stanza_key,
                                      unsigned long num);
unsigned long gt_obo_stanza_size(const GtOBOStanza *obo_stanza,
                                 const char *stanza_key);
const char*   gt_obo_stanza_filename(const GtOBOStanza *obo_stanza);
unsigned long gt_obo_stanza_line(const GtOBOStanza *obo_stanza);

#endif
