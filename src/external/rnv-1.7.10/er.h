/* $Id: er.h 315 2004-01-20 00:22:29Z dvd $ */

#ifndef ER_H
#define ER_H 1

#include <stdarg.h>

extern int (*er_printf)(char *format,...);
extern int (*er_vprintf)(char *format,va_list ap);

extern int er_default_printf(char *format,...);
extern int er_default_vprintf(char *format,va_list ap);

#endif
