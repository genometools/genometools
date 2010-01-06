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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "core/unused_api.h"
#include "core/warning_api.h"

static GtWarningHandler warning_handler = gt_warning_default_handler;
static void *warning_data = NULL;

void gt_warning(const char *format, ...)
{
  va_list ap;
  if (warning_handler) {
    va_start(ap, format);
    warning_handler(warning_data, format, ap);
    va_end(ap);
  }
}

void gt_warning_disable(void)
{
  warning_handler = NULL;
  warning_data = NULL;
}

void gt_warning_set_handler(GtWarningHandler warn_handler, void *data)
{
  warning_handler = warn_handler;
  warning_data = data;
}

void gt_warning_default_handler(GT_UNUSED void *data, const char *format,
                                va_list ap)
{
  (void) fputs("warning: ", stderr);
  (void) vfprintf(stderr, format, ap);
  (void) putc('\n', stderr);
}

GtWarningHandler gt_warning_get_handler(void)
{
  return warning_handler;
}

void* gt_warning_get_data(void)
{
  return warning_data;
}
