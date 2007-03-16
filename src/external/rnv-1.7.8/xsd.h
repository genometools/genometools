/* $Id: xsd.h,v 1.14 2005/01/05 09:46:25 dvd Exp $ */

#include <stdarg.h>

#ifndef XSD_H
#define XSD_H 1

#define XSD_ER_TYP 0
#define XSD_ER_PAR 1
#define XSD_ER_PARVAL 2
#define XSD_ER_VAL 3
#define XSD_ER_NPAT 4
#define XSD_ER_WS 5
#define XSD_ER_ENUM 6

extern void (*xsd_verror_handler)(int erno,va_list ap);

extern void xsd_default_verror_handler(int erno,va_list ap);

extern void xsd_init(void);
extern void xsd_clear(void);

extern int xsd_allows(char *typ,char *ps,char *s,int n);
extern int xsd_equal(char *typ,char *val,char *s,int n);

extern void xsd_test(void);

#endif
