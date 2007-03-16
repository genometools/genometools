/* $Id: rnd.h,v 1.11 2004/01/02 20:27:23 dvd Exp $ */

#include <stdarg.h>

#ifndef RND_H
#define RND_H 1

#define RND_ER_LOOPST 0
#define RND_ER_LOOPEL 1
#define RND_ER_CTYPE 2
#define RND_ER_BADSTART 3
#define RND_ER_BADMORE 4
#define RND_ER_BADEXPT 5
#define RND_ER_BADLIST 6
#define RND_ER_BADATTR 7

extern void (*rnd_verror_handler)(int er_no,va_list ap);

extern void rnd_default_verror_handler(int erno,va_list ap);

extern void rnd_init(void);
extern void rnd_clear(void);

extern int rnd_fixup(int start);

#endif
