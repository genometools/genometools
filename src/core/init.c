/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifdef HAVE_MYSQL
#include <mysql/mysql.h>
#endif
#ifndef WITHOUT_CAIRO
/* #include <fontconfig.h> */
#endif
#include <string.h>
#include "core/class_alloc.h"
#include "core/class_alloc_lock.h"
#include "core/combinatorics.h"
#include "core/cstr_api.h"
#include "core/cstr_array.h"
#include "core/fa.h"
#include "core/init_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/showtime.h"
#include "core/spacepeak.h"
#include "core/splitter.h"
#include "core/symbol.h"
#include "core/versionfunc.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "core/yarandom.h"

static bool spacepeak = false;
static bool showtime = false;

static GtOPrval parse_env_options(int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *o;
  GtOPrval oprval;
  op = gt_option_parser_new("GT_ENV_OPTIONS='[option ...]' ...",
                         "Parse the options contained in the "
                         "environment variable GT_ENV_OPTIONS.");
  o = gt_option_new_bool("spacepeak", "show space peak on stdout upon deletion",
                         &spacepeak, false);
  gt_option_parser_add_option(op, o);
  o = gt_option_new_bool("showtime", "enable output for run-time statistics",
                         &showtime, false);
  gt_option_parser_add_option(op, o);
  gt_option_parser_set_max_args(op, 0);
  oprval = gt_option_parser_parse(op, NULL, argc, argv, gt_versionfunc, err);
  gt_option_parser_delete(op);
  return oprval;
}

static void proc_env_options(void)
{
  int argc;
  char *env_options, **argv;
  GtSplitter *splitter;
  GtError *err;
  /* construct argument vector from $GT_ENV_OPTIONS */
  env_options = getenv("GT_ENV_OPTIONS");
  if (!env_options)
    return;
  env_options = gt_cstr_dup(env_options); /* make writeable copy */
  splitter = gt_splitter_new();
  gt_splitter_split(splitter, env_options, strlen(env_options), ' ');
  argc = gt_splitter_size(splitter);
  argv = gt_cstr_array_preprend((const char**) gt_splitter_get_tokens(splitter),
                             "env");
  argc++;
  /* parse options contained in $GT_ENV_OPTIONS */
  err = gt_error_new();
  switch (parse_env_options(argc, (const char**) argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      fprintf(stderr, "error parsing $GT_ENV_OPTIONS: %s\n", gt_error_get(err));
      gt_error_unset(err);
      break;
    case GT_OPTION_PARSER_REQUESTS_EXIT: break;
  }
  gt_error_delete(err);
  gt_free(env_options);
  gt_splitter_delete(splitter);
  gt_cstr_array_delete(argv);
}

void gt_lib_init(void)
{
  const char *bookkeeping;
  bookkeeping = getenv("GT_MEM_BOOKKEEPING");
  gt_ma_init(bookkeeping && !strcmp(bookkeeping, "on"));
  proc_env_options();
  if (spacepeak && !(bookkeeping && !strcmp(bookkeeping, "on")))
    gt_warning("GT_ENV_OPTIONS=-spacepeak used without GT_MEM_BOOKKEEPING=on");
  gt_fa_init();
  if (spacepeak) {
    gt_spacepeak_init();
    gt_ma_enable_global_spacepeak();
    gt_fa_enable_global_spacepeak();
  }
  gt_log_init();
  if (showtime) gt_showtime_enable();
  gt_symbol_init();
  gt_class_alloc_lock_init();
  gt_ya_rand_init(0);
#ifdef HAVE_MYSQL
  mysql_library_init(0, NULL, NULL);
#endif
  gt_combinatorics_init();
}

static void gt_lib_atexit_func(void)
{
  (void) gt_lib_clean();
}

void gt_lib_reg_atexit_func(void)
{
  gt_xatexit(gt_lib_atexit_func);
}

int gt_lib_clean(void)
{
  int fa_fptr_rval, fa_mmap_rval, gt_rval;
  if (spacepeak) {
    gt_ma_show_space_peak(stdout);
    gt_fa_show_space_peak(stdout);
    gt_spacepeak_show_space_peak(stdout);
    gt_ma_disable_global_spacepeak();
  }
  fa_fptr_rval = gt_fa_check_fptr_leak();
  fa_mmap_rval = gt_fa_check_mmap_leak();
  gt_fa_clean();
  gt_symbol_clean();
  gt_class_alloc_clean();
  gt_class_alloc_lock_clean();
  gt_ya_rand_clean();
  gt_log_clean();
  gt_spacepeak_clean();
  gt_combinatorics_clean();
  gt_rval = gt_ma_check_space_leak();
  gt_ma_clean();
#ifdef HAVE_MYSQL
  mysql_library_end();
#endif
  return fa_fptr_rval || fa_mmap_rval || gt_rval;
}
