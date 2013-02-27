/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_table.h"
#include "core/mathsupport.h"
#include "core/symbol.h"
#include "core/thread_api.h"
#include "core/unused_api.h"

static GtCstrTable *symbols = NULL;
static GtMutex *symbol_mutex = NULL;

void gt_symbol_init(void)
{
  if (!symbols)
    symbols = gt_cstr_table_new();
  if (!symbol_mutex)
    symbol_mutex = gt_mutex_new();
}

const char* gt_symbol(const char *cstr)
{
  const char *symbol;
  if (!cstr)
    return NULL;
  gt_mutex_lock(symbol_mutex);
  if (!(symbol = gt_cstr_table_get(symbols, cstr))) {
    gt_cstr_table_add(symbols, cstr);
    symbol = gt_cstr_table_get(symbols, cstr);
  }
  gt_mutex_unlock(symbol_mutex);
  return symbol;
}

void gt_symbol_clean(void)
{
  gt_cstr_table_delete(symbols);
  gt_mutex_delete(symbol_mutex);
}

/* we use randomly generated numbers to test the symbol mechanism */
#define NUMBER_OF_SYMBOLS 10000
#define MAX_SYMBOL        1000

static void* test_symbol(GT_UNUSED void *data)
{
  GtStr *symbol;
  unsigned long i;
  symbol = gt_str_new();
  for (i = 0; i < NUMBER_OF_SYMBOLS; i++) {
    gt_str_reset(symbol);
    gt_str_append_ulong(symbol, gt_rand_max(MAX_SYMBOL));
    gt_symbol(gt_str_get(symbol));
    gt_assert(!strcmp(gt_symbol(gt_str_get(symbol)), gt_str_get(symbol)));
  }
  gt_str_delete(symbol);
  return NULL;
}

int gt_symbol_unit_test(GtError *err)
{
  gt_error_check(err);
  return gt_multithread(test_symbol, NULL, err);
}
