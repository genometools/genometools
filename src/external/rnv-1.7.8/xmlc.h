/* $Id: xmlc.h,v 1.2 2003/12/21 22:38:28 dvd Exp $ */

#ifndef XMLC_H
#define XMLC_H 1

/* character classes required for parsing XML */
extern int xmlc_white_space(int u);
extern int xmlc_base_char(int u);
extern int xmlc_ideographic(int u);
extern int xmlc_combining_char(int u);
extern int xmlc_digit(int u);
extern int xmlc_extender(int u);

extern int u_in_ranges(int u,int r[][2],int len);

#endif
