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

#ifndef XRF_ABBR_PARSE_TREE_H
#define XRF_ABBR_PARSE_TREE_H

#include "core/error.h"
#include "extended/xrf_abbr_entry.h"

typedef struct GtXRFAbbrParseTree GtXRFAbbrParseTree;

/* Parse the XRF abbreviation file given with <file_path> and return the
   result as a GtXRFAbbrParseTree. If an error occurs during parsing,
   NULL is returned and <err> is set. */
GtXRFAbbrParseTree*   gt_xrf_abbr_parse_tree_new(const char *file_path,
                                                 GtError *err);
void                  gt_xrf_abbr_parse_tree_delete(GtXRFAbbrParseTree*);
/* Return the type of entry number <entry_number>. */
const char*           gt_xrf_abbr_parse_tree_get_entry_type(const
                                                     GtXRFAbbrParseTree*,
                                                     GtUword entry_num);
/* Return the value of entry <entry_key> in entry number <entry_number>. */
const char*           gt_xrf_abbr_parse_tree_get_entry_value(const
                                                      GtXRFAbbrParseTree*,
                                                      GtUword entry_num,
                                                      const char *entry_key);
/* Return OBO entry with number <entry_number. */
const GtXRFAbbrEntry* gt_xrf_abbr_parse_tree_get_entry(const
                                                       GtXRFAbbrParseTree*,
                                                       GtUword entry_num);
/* Return the number of entries. */
GtUword               gt_xrf_abbr_parse_tree_num_of_entries(const
                                                       GtXRFAbbrParseTree*);

#endif
