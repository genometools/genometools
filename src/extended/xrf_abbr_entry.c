/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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
#include "core/str_api.h"
#include "extended/xrf_abbr_entry.h"

struct GtXRFAbbrEntry {
  GtHashmap *content;
  GtUword line;
  GtStr *filename;
};

GtXRFAbbrEntry* gt_xrf_abbr_entry_new(GtUword line,
                                      GtStr *filename)
{
  GtXRFAbbrEntry *abbr_entry = gt_malloc(sizeof *abbr_entry);
  abbr_entry->content = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                       (GtFree) gt_str_delete);
  abbr_entry->line = line;
  abbr_entry->filename = gt_str_ref(filename);
  return abbr_entry;
}

void gt_xrf_abbr_entry_delete(GtXRFAbbrEntry *abbr_entry)
{
  if (!abbr_entry) return;
  gt_str_delete(abbr_entry->filename);
  gt_hashmap_delete(abbr_entry->content);
  gt_free(abbr_entry);
}

void gt_xrf_abbr_entry_add(GtXRFAbbrEntry *abbr_entry, const char *tag,
                           const char *value)
{
  GtStr *s;
  gt_assert(abbr_entry && tag && value);
  if (!(s = gt_hashmap_get(abbr_entry->content, tag))) {
    s = gt_str_new_cstr(value);
    gt_hashmap_add(abbr_entry->content, gt_cstr_dup(tag), s);
  }
}

const char* gt_xrf_abbr_entry_get_value(const GtXRFAbbrEntry *abbr_entry,
                                        const char *entry_key)
{
  GtStr *s;
  gt_assert(abbr_entry);
  if ((s = gt_hashmap_get(abbr_entry->content, entry_key)))
    return gt_str_get(s);
  return NULL;
}

const char* gt_xrf_abbr_entry_filename(const GtXRFAbbrEntry *abbr_entry)
{
  gt_assert(abbr_entry);
  return gt_str_get(abbr_entry->filename);
}

GtUword gt_xrf_abbr_entry_line(const GtXRFAbbrEntry *abbr_entry)
{
  gt_assert(abbr_entry);
  return abbr_entry->line;
}
