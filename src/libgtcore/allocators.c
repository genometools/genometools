/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/allocators.h"
#include "libgtcore/cstr.h"
#include "libgtcore/cstr_array.h"
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/splitter.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/warning.h"
#include "libgtcore/xansi.h"

static bool spacepeak = false;

static OPrval parse_env_options(int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  op = option_parser_new("GT_ENV_OPTIONS='[option ...]' ...",
                         "Parse the options contained in the "
                         "environment variable GT_ENV_OPTIONS.");
  o = option_new_bool("spacepeak", "show space peak on stdout upon deletion",
                      &spacepeak, false);
  option_parser_add_option(op, o);
  option_parser_set_max_args(op, 0);
  oprval = option_parser_parse(op, NULL, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void proc_gt_env_options(void)
{
  int argc;
  char *env_options, **argv;
  Splitter *splitter;
  Error *err;
  /* construct argument vector from $GT_ENV_OPTIONS */
  env_options = getenv("GT_ENV_OPTIONS");
  if (!env_options)
    return;
  env_options = cstr_dup(env_options); /* make writeable copy */
  splitter = splitter_new();
  splitter_split(splitter, env_options, strlen(env_options), ' ');
  argc = splitter_size(splitter);
  argv = cstr_array_preprend((const char**) splitter_get_tokens(splitter),
                             "env");
  argc++;
  /* parse options contained in $GT_ENV_OPTIONS */
  err = error_new();
  switch (parse_env_options(argc, (const char**) argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      fprintf(stderr, "error parsing $GT_ENV_OPTIONS: %s\n", error_get(err));
      error_unset(err);
      break;
    case OPTIONPARSER_REQUESTS_EXIT: break;
  }
  error_delete(err);
  ma_free(env_options);
  splitter_delete(splitter);
  cstr_array_delete(argv);
}

void allocators_init(void)
{
  const char *bookkeeping;
  bookkeeping = getenv("GT_MEM_BOOKKEEPING");
  ma_init(bookkeeping && !strcmp(bookkeeping, "on"));
  proc_gt_env_options();
  if (spacepeak && !(bookkeeping && !strcmp(bookkeeping, "on")))
    warning("GT_ENV_OPTIONS=-spacepeak used without GT_MEM_BOOKKEEPING=on");
}

static void allocators_atexit_func(void)
{
  (void) allocators_clean();
}

void allocators_reg_atexit_func(void)
{
  xatexit(allocators_atexit_func);
}

int allocators_clean(void)
{
  int fa_fptr_rval, fa_mmap_rval, ma_rval;
  if (spacepeak) {
    ma_show_space_peak(stdout);
    fa_show_space_peak(stdout);
  }
  fa_fptr_rval = fa_check_fptr_leak();
  fa_mmap_rval = fa_check_mmap_leak();
  fa_clean();
  ma_rval = ma_check_space_leak();
  ma_clean();
  return fa_fptr_rval || fa_mmap_rval || ma_rval;
}
