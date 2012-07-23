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

#include "core/cstr_api.h"
#include "core/fptr_api.h"
#include "core/hashmap_api.h"
#include "core/ma_api.h"
#include "core/str_array_api.h"
#include "extended/obo_stanza.h"

struct GtOBOStanza {
  char *type;
  GtHashmap *content;
  unsigned long line;
  GtStr *filename;
};

GtOBOStanza* gt_obo_stanza_new(const char *type, unsigned long line,
                               GtStr *filename)
{
  GtOBOStanza *obo_stanza = gt_malloc(sizeof *obo_stanza);
  obo_stanza->type = gt_cstr_dup(type);
  obo_stanza->content = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                       (GtFree) gt_str_array_delete);
  obo_stanza->line = line;
  obo_stanza->filename = gt_str_ref(filename);
  return obo_stanza;
}

void gt_obo_stanza_delete(GtOBOStanza *obo_stanza)
{
  if (!obo_stanza) return;
  gt_str_delete(obo_stanza->filename);
  gt_hashmap_delete(obo_stanza->content);
  gt_free(obo_stanza->type);
  gt_free(obo_stanza);
}

void gt_obo_stanza_add(GtOBOStanza *obo_stanza, const char *tag,
                       const char *value)
{
  GtStrArray *sa;
  gt_assert(obo_stanza && tag && value);
  if (!(sa = gt_hashmap_get(obo_stanza->content, tag))) {
    sa = gt_str_array_new();
    gt_str_array_add_cstr(sa, value);
    gt_hashmap_add(obo_stanza->content, gt_cstr_dup(tag), sa);
  }
  else
    gt_str_array_add_cstr(sa, value);
}

const char* gt_obo_stanza_get_type(const GtOBOStanza *obo_stanza)
{
  gt_assert(obo_stanza);
  return obo_stanza->type;
}

const char* gt_obo_stanza_get_value(const GtOBOStanza *obo_stanza,
                                    const char *stanza_key,
                                    unsigned long num)
{
  GtStrArray *sa;
  gt_assert(obo_stanza);
  if ((sa = gt_hashmap_get(obo_stanza->content, stanza_key)))
    return gt_str_array_get(sa, num);
  return NULL;
}

unsigned long gt_obo_stanza_size(const GtOBOStanza *obo_stanza,
                                 const char *stanza_key)
{
  GtStrArray *sa;
  gt_assert(obo_stanza);
  if ((sa = gt_hashmap_get(obo_stanza->content, stanza_key)))
    return gt_str_array_size(sa);
  return 0;
}

const char* gt_obo_stanza_filename(const GtOBOStanza *obo_stanza)
{
  gt_assert(obo_stanza);
  return gt_str_get(obo_stanza->filename);
}

unsigned long gt_obo_stanza_line(const GtOBOStanza *obo_stanza)
{
  gt_assert(obo_stanza);
  return obo_stanza->line;
}
