/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdarg.h>
#include "core/logger.h"

static GtLogger *gt_global_logger;

void gt_log_init()
{
  gt_global_logger = gt_logger_new(false, "debug: ", stderr);
}

void gt_log_enable(void)
{
  gt_logger_enable(gt_global_logger);
}

bool gt_log_enabled(void)
{
  return gt_logger_enabled(gt_global_logger);
}

void gt_log_log(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  gt_logger_log_va(gt_global_logger, format, ap);
  va_end(ap);
}

void gt_log_vlog(const char *format, va_list ap)
{
  gt_logger_log_va(gt_global_logger, format, ap);
}

FILE* gt_log_fp(void)
{
  return gt_logger_target(gt_global_logger);
}

void gt_log_set_fp(FILE *fp)
{
  gt_logger_set_target(gt_global_logger, fp);
}

void gt_log_clean()
{
  gt_logger_delete(gt_global_logger);
}
