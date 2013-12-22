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

#include <string.h>
#include "core/cstr_api.h"
#include "core/grep_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/splitter_api.h"
#include "core/unused_api.h"
#include "extended/xrf_abbr_entry.h"
#include "extended/xrf_abbr_parse_tree.h"
#include "extended/xrf_checker_api.h"

struct GtXRFChecker {
  GtHashmap *abbrvs;
  GtXRFAbbrParseTree *xpt;
  GtSplitter *splitter;
  GtUword reference_count;
};

GtXRFChecker* gt_xrf_checker_ref(GtXRFChecker *xrc)
{
  if (!xrc) return NULL;
  xrc->reference_count++;
  return xrc;
}

bool gt_xrf_checker_is_valid(GtXRFChecker *xrc, const char *value, GtError *err)
{
  bool valid = true;
  char *myvalue = gt_cstr_dup(value),
       *dbid = NULL,
       *localid = NULL;
  GtXRFAbbrEntry *e;
  GtUword nof_tokens, i;
  gt_assert(xrc && value);
  gt_error_check(err);

  /* XXX: Thread safety! */
  gt_splitter_reset(xrc->splitter);
  gt_splitter_split(xrc->splitter, myvalue, strlen(myvalue), ',');
  nof_tokens = gt_splitter_size(xrc->splitter);

  for (i = 0; valid && i < nof_tokens; i++) {
    dbid = gt_splitter_get_token(xrc->splitter, i);

    if (!(localid = strchr(dbid, ':'))) {
      gt_error_set(err, "xref \"%s\": separator colon missing", value);
      valid = false;
    }
    if (valid) {
      *localid = '\0';
      if (*(++localid) == '\0') {
        gt_error_set(err, "xref \"%s\": local ID (part after colon) missing",
                     value);
        valid = false;
      }
    }
    if (valid) {
      gt_assert(dbid && localid);
      if (!(e = gt_hashmap_get(xrc->abbrvs, dbid))) {
        gt_error_set(err, "xref \"%s\": unknown database abbreviation \"%s\"",
                     value, dbid);
        valid = false;
      }
    }
    if (valid) {
      /* TODO: use #defines here. */
      const char *regex = NULL;
      gt_assert(e);
      if ((regex = gt_xrf_abbr_entry_get_value(e, "local_id_syntax"))) {
        bool match = false;
        GT_UNUSED int rval;
        rval = gt_grep(&match, regex, localid, NULL);
        gt_assert(rval == 0); /* regex format has been checked before */
        if (!match) {
          gt_error_set(err, "xref \"%s\": local ID \"%s\" does not "
                            "conform to syntax \"%s\" for the "
                            "%s database",
                       value, localid, regex, dbid);
          valid = false;
        }
      }
    }
  }

  gt_free(myvalue);
  return valid;
}

GtXRFChecker* gt_xrf_checker_new(const char *file_path, GtError *err)
{
  GtXRFChecker *xrc;
  GtUword i;
  gt_error_check(err);
  gt_assert(file_path);

  xrc = gt_calloc(1UL, sizeof (GtXRFChecker));
  xrc->xpt = gt_xrf_abbr_parse_tree_new(file_path, err);
  if (!xrc->xpt) {
    gt_xrf_checker_delete(xrc);
    return NULL;
  }
  xrc->abbrvs = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  for (i = 0; i < gt_xrf_abbr_parse_tree_num_of_entries(xrc->xpt); i++) {
    const GtXRFAbbrEntry *e = gt_xrf_abbr_parse_tree_get_entry(xrc->xpt, i);
    const char *synonym;
    gt_hashmap_add(xrc->abbrvs,
                   (void*) gt_xrf_abbr_entry_get_value(e, "abbreviation"),
                   (void*) e);
    if ((synonym = gt_xrf_abbr_entry_get_value(e, "synonym"))) {
      gt_hashmap_add(xrc->abbrvs, (void*) synonym, (void*) e);
    }
  }
  xrc->splitter = gt_splitter_new();
  return xrc;
}

void gt_xrf_checker_delete(GtXRFChecker *xrc)
{
  if (!xrc) return;
  if (xrc->reference_count) {
    xrc->reference_count--;
    return;
  }
  gt_xrf_abbr_parse_tree_delete(xrc->xpt);
  gt_hashmap_delete(xrc->abbrvs);
  gt_splitter_delete(xrc->splitter);
  gt_free(xrc);
}
