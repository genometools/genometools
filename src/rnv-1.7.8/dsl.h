/* $Id: dsl.h,v 1.3 2004/01/28 12:34:12 dvd Exp $ */

#ifndef DSL_H
#define DSL_H 1

#ifndef DSL_SCM
#define DSL_SCM 0
#endif

#define DSL_URL "http://davidashen.net/relaxng/scheme-datatypes"

extern void dsl_ld(char *dl);

extern int dsl_allows(char *typ,char *ps,char *s,int n);
extern int dsl_equal(char *typ,char *val,char *s,int n);

#endif
