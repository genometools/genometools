/* $Id: dxl.h,v 1.2 2004/01/28 12:34:12 dvd Exp $ */

#ifndef DXL_H
#define DXL_H 1

#ifndef DXL_EXC
#define DXL_EXC 0
#endif

#define DXL_URL "http://davidashen.net/relaxng/pluggable-datatypes"

extern char *dxl_cmd;

extern int dxl_allows(char *typ,char *ps,char *s,int n);
extern int dxl_equal(char *typ,char *val,char *s,int n);

#endif
