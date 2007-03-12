/* $Id: er.c,v 1.24 2004/01/20 00:25:09 dvd Exp $ */

#include <stdio.h>
#include "er.h"

int (*er_printf)(char *format,...)=&er_default_printf;
int (*er_vprintf)(char *format,va_list ap)=&er_default_vprintf;

int er_default_printf(char *format,...) {
  int ret;
  va_list ap; va_start(ap,format); ret=(*er_vprintf)(format,ap); va_end(ap);
  return ret;
}
int er_default_vprintf(char *format,va_list ap) {return vfprintf(stderr,format,ap);}
