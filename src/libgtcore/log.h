/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LOG_H
#define LOG_H

#include <stdarg.h>
#include <stdio.h>

typedef struct Log Log;

/**
 * Creates a new log object, the output is directed at stderr by default.
 * @param memory allocator
 * @return log object
 */
Log*  log_new(MA*);
/**
 * Prints the log message obtained from format and following
 * parameters according to Log object. The debug output is typically
 * prefixed with the string "debug: " and finished by a new-line
 */
void  log_log(Log*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));
/**
 * Prints the log message obtained from format and following
 * parameter according to Log object like log_log. But in contrast to
 * log_log log_vlog not accepts individual arguments but a single
 * va_list argument instead.
 */
void  log_vlog(Log*, const char *format, va_list);
FILE* log_fp(Log*);
void  log_delete(Log*, MA*);

#endif
