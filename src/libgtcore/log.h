/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef LOG_H
#define LOG_H

#include <stdarg.h>
#include <stdio.h>

typedef struct Log Log;

/*
  Creates a new log object, the output is directed at stderr by default.
  @param memory allocator
  @return log object
*/
Log*  log_new(MA*);
/*
  Prints the log message obtained from format and following
  parameters according to Log object. The debug output is
  prefixed with the string "debug: " and finished by a newline.
*/
void  log_log(Log*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));
/*
  Prints the log message obtained from format and following
  parameter according to Log object like log_log. But in contrast to
  log_log log_vlog not accepts individual arguments but a single
  va_list argument instead.
*/
void  log_vlog(Log*, const char *format, va_list);
FILE* log_fp(Log*);
void  log_delete(Log*, MA*);

#endif
