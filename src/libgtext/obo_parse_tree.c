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

#include "libgtcore/array.h"
#include "libgtcore/fa.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/obo_parse_tree.h"

typedef struct {
  Hashtable *content;
} OBOHeader;

static OBOHeader* obo_header_new(void)
{
  OBOHeader *obo_header = ma_malloc(sizeof *obo_header);
  obo_header->content = hashtable_new(HASH_STRING, ma_free_func, ma_free_func);
  return obo_header;
}

static void obo_header_delete(OBOHeader *obo_header)
{
  if (!obo_header) return;
  hashtable_delete(obo_header->content);
  ma_free(obo_header);
}

typedef struct {
  Hashtable *content;
} OBOStanza;

#if 0
static OBOStanza* obo_stanza_new(void)
{
  OBOStanza *obo_stanza = ma_malloc(sizeof *obo_stanza);
  obo_stanza->content = hashtable_new(HASH_STRING, ma_free_func, ma_free_func);
  return obo_stanza;
}
#endif

static void obo_stanza_delete(OBOStanza *obo_stanza)
{
  if (!obo_stanza) return;
  hashtable_delete(obo_stanza->content);
  ma_free(obo_stanza);
}

static const char* obo_stanza_get_value(const OBOStanza *obo_stanza,
                                        const char *stanza_key)
{
  assert(obo_stanza);
  return hashtable_get(obo_stanza->content, stanza_key);
}

struct OBOParseTree {
  OBOHeader *obo_header;
  Array *stanzas;
};

static int parse_obo_file(UNUSED OBOParseTree *obo_parse_tree,
                          UNUSED FILE *obo_file, UNUSED Error *err)
{
  int had_err = 0;
  error_check(err);
  assert(obo_parse_tree && obo_file);
  /* XXX */
  return had_err;
}

OBOParseTree* obo_parse_tree_new(const char *obo_file_path, Error *err)
{
  OBOParseTree *obo_parse_tree;
  FILE *obo_file;
  error_check(err);
  assert(obo_file_path);
  obo_file = fa_xfopen(obo_file_path, "r");
  obo_parse_tree = ma_malloc(sizeof *obo_parse_tree);
  obo_parse_tree->obo_header = obo_header_new();
  obo_parse_tree->stanzas = array_new(sizeof (OBOStanza*));
  if (parse_obo_file(obo_parse_tree, obo_file, err)) {
    obo_parse_tree_delete(obo_parse_tree);
    fa_fclose(obo_file);
    return NULL;
  }
  fa_xfclose(obo_file);
  return obo_parse_tree;
}

void obo_parse_tree_delete(OBOParseTree *obo_parse_tree)
{
  unsigned long i;
  if (!obo_parse_tree) return;
  for (i = 0; i < array_size(obo_parse_tree->stanzas); i++)
    obo_stanza_delete(*(OBOStanza**) array_get(obo_parse_tree->stanzas, i));
  array_delete(obo_parse_tree->stanzas);
  obo_header_delete(obo_parse_tree->obo_header);
  ma_free(obo_parse_tree);
}

const char* obo_parse_tree_get_stanza_value(OBOParseTree *obo_parse_tree,
                                            unsigned long stanza_num,
                                            const char *stanza_key)
{
  assert(obo_parse_tree);
  return obo_stanza_get_value(*(OBOStanza**)
                              array_get(obo_parse_tree->stanzas, stanza_num),
                              stanza_key);
}

unsigned long obo_parse_tree_num_of_stanzas(const OBOParseTree *obo_parse_tree)
{
  assert(obo_parse_tree);
  return array_size(obo_parse_tree->stanzas);
}
