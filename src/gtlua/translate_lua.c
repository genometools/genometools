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

#include <string.h>
#include "lauxlib.h"
#include "core/codon_iterator_simple_api.h"
#include "core/translator.h"
#include "gtlua/translate_lua.h"

static int translate_dna_lua(lua_State *L)
{
  GtStr *protein;
  GtTranslator *tr;
  int rval;
  char translated;
  unsigned int frame;
  const char *dna = luaL_checkstring(L, 1);
  protein = gt_str_new();

  GtCodonIterator *ci = gt_codon_iterator_simple_new(dna,
                                                     strlen(dna),
                                                     NULL);
  tr = gt_translator_new(ci);
  rval = gt_translator_next(tr, &translated, &frame, NULL);
  while (!rval && translated) {
    gt_str_append_char(protein, translated);
    rval = gt_translator_next(tr, &translated, &frame, NULL);
  }
  lua_pushstring(L, gt_str_get(protein));
  gt_str_delete(protein);
  gt_translator_delete(tr);
  gt_codon_iterator_delete(ci);
  return 1;
}

static const struct luaL_Reg translate_lib_f [] = {
  { "translate_dna", translate_dna_lua },
  { NULL, NULL }
};

int gt_lua_open_translate(lua_State *L)
{
#ifndef NDEBUG
  int stack_size;
#endif
  gt_assert(L);
#ifndef NDEBUG
  stack_size = lua_gettop(L);
#endif
  luaL_register(L, "gt", translate_lib_f);
  lua_pop(L, 1);
  gt_assert(lua_gettop(L) == stack_size);
  return 1;
}
