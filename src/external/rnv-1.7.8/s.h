/* $Id: s.h,v 1.6 2005/01/05 09:46:25 dvd Exp $ */

#ifndef S_H
#define S_H 1

/* compares two strings, s1 is null terminated, s2 is n characters long */
extern int s_cmpn(char *s1,char *s2,int n2);

/* compares two tokens, s1 is null terminated, s2 is n characters long */
extern int s_tokcmpn(char *s1,char *s2,int n2);

/* hash value for a zero-terminated string */
extern int s_hval(char *s);

/* strdup is a non-standard function */
extern char *s_clone(char *s);

/* compute the absolute path from a relative path and a base path;
 the caller must ensure that there is enough space in r:
 size(r) > strlen(r)+strlen(b)
 returns a pointer to the string containing the relative path
 */
extern char *s_abspath(char *r,char *b);

/* find a string in a sorted array, return the index,
  or size on failure */
extern int s_tab(char *s,char *tab[],int size);
extern int s_ntab(char *s,int len,char *tab[],int size);

extern void s_test(void);

#endif
