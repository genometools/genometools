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

Log*  log_new(MA*);
void  log_log(Log*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));
void  log_vlog(Log*, const char *format, va_list);
FILE* log_fp(Log*);
void  log_delete(Log*, MA*);

#endif
