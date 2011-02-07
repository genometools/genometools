/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/logger.h"
#include "core/ma_api.h"

struct GtLogger
{
  FILE *target;
  char *prefix;
  bool enabled;
};

GtLogger* gt_logger_new(bool enabled, const char *prefix, FILE *target)
{
  GtLogger *l;
  gt_assert(target);
  l = (GtLogger*) gt_calloc((size_t) 1, sizeof (GtLogger));
  l->target = target;
  l->prefix = gt_cstr_dup(prefix);
  l->enabled = enabled;
  return l;
}

void gt_logger_enable(GtLogger *logger)
{
  gt_assert(logger);
  logger->enabled = true;
}

bool gt_logger_enabled(GtLogger *logger)
{
  gt_assert(logger);
  return logger->enabled;
}

void gt_logger_disable(GtLogger *logger)
{
  gt_assert(logger);
  logger->enabled = false;
}

FILE* gt_logger_target(GtLogger *logger)
{
  gt_assert(logger);
  return logger->target;
}

void gt_logger_set_target(GtLogger *logger, FILE *fp)
{
  gt_assert(logger && fp);
  logger->target = fp;
}

void gt_logger_log_force(GtLogger *logger, const char *format, ...)
{
  va_list ap;
  if (!logger) return;
  gt_assert(format);
  va_start(ap, format);
  gt_logger_log_va_force(logger, format, ap);
  va_end(ap);
}

void gt_logger_log(GtLogger *logger, const char *format, ...)
{
  va_list ap;
  if (!logger) return;
  if (!logger->enabled) return;
  gt_assert(format);
  va_start(ap, format);
  gt_logger_log_va(logger, format, ap);
  va_end(ap);
}

void gt_logger_log_va_force(GtLogger *logger, const char *format, va_list ap)
{
  if (!logger) return;
  gt_assert(format && logger->target);
  if (logger->prefix != NULL)
    fprintf(logger->target, "%s", logger->prefix);
  (void) vfprintf(logger->target, format, ap);
  (void) putc('\n', logger->target);
}

void gt_logger_log_va(GtLogger *logger, const char *format, va_list ap)
{
  if (!logger) return;
  if (!logger->enabled) return;
  gt_assert(format && logger->target);
  if (logger->prefix != NULL)
    fprintf(logger->target, "%s", logger->prefix);
  (void) vfprintf(logger->target, format, ap);
  (void) putc('\n', logger->target);
}

void gt_logger_delete(GtLogger *logger)
{
  if (!logger) return;
  gt_free(logger->prefix);
  gt_free(logger);
}
