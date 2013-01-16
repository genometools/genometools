/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef OBO_PARSE_TREE_H
#define OBO_PARSE_TREE_H

#include "core/error.h"
#include "extended/obo_stanza.h"

typedef struct GtOBOParseTree GtOBOParseTree;

/* Parse the OBO file given with <obo_file_path> and return the result as an
   GtOBOParseTree.
   If an error occurs during parsing, NULL is returned and <err> is set. */
GtOBOParseTree*    gt_obo_parse_tree_new(const char *obo_file_path,
                                         GtError *err);
void               gt_obo_parse_tree_delete(GtOBOParseTree*);
/* Return the type of stanza number <stanza_number>. */
const char*        gt_obo_parse_tree_get_stanza_type(const GtOBOParseTree*,
                                                     unsigned long stanza_num);
/* Return the value of entry <stanza_key> in stanza number <stanza_number>. */
const char*        gt_obo_parse_tree_get_stanza_value(const GtOBOParseTree*,
                                                      unsigned long stanza_num,
                                                      const char *stanza_key);
/* Return OBO stanza with number <stanza_number. */
const GtOBOStanza* gt_obo_parse_tree_get_stanza(const GtOBOParseTree*,
                                                unsigned long stanza_num);
/* Return the number of stanzas. */
unsigned long      gt_obo_parse_tree_num_of_stanzas(const GtOBOParseTree*);

#endif
