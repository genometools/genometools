/* $Id: u.h,v 1.15 2004/01/15 23:47:45 dvd Exp $ */

#ifndef U_H
#define U_H 1

#define U_MAXLEN 6

/* returns BOM length if the string starts with BOM */
extern int u_bom(char *s,int n);

/* computes a unicode character u off the head of s;
 returns number of bytes read. 0 means error.
 */
extern int u_get(int *up,char *s);

/* encodes u in utf-8, returns number of octets taken */
extern int u_put(char *s,int u);

/* number of unicode characters in the string; -1 means error */
extern int u_strlen(char *s);
extern int u_strnlen(char *s,int n);

/* checks whether a character falls within one of sorted ranges */
extern int u_in_ranges(int u,int r[][2],int len);

#endif
