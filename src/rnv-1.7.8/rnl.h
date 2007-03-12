/* $Id: rnl.h,v 1.1 2004/01/10 14:25:05 dvd Exp $ */

#ifndef RNL_H
#define RNL_H 1

extern void (*rnl_verror_handler)(int er_no,va_list ap);
extern void rnl_default_verror_handler(int erno,va_list ap);

extern void rnl_init(void);
extern void rnl_clear(void);

extern int rnl_fn(char *fn);
extern int rnl_fd(char *fn,int fd);
extern int rnl_s(char *fn,char *s,int len);

#endif
