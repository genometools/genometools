/* $Id: rnl.c,v 1.2 2004/01/15 23:47:45 dvd Exp $ */

#include <stdarg.h>
#include "erbit.h"
#include "rn.h"
#include "rnc.h"
#include "rnd.h"
#include "rnl.h"

void rnl_default_verror_handler(int erno,va_list ap) {
  if(erno&ERBIT_RNC) {
    rnc_default_verror_handler(erno&~ERBIT_RNC,ap);
  } else if(erno&ERBIT_RND) {
    rnd_default_verror_handler(erno&~ERBIT_RND,ap);
  }
}
void (*rnl_verror_handler)(int er_no,va_list ap)=&rnl_default_verror_handler;

static void verror_handler_rnc(int erno,va_list ap) {rnl_verror_handler(erno|ERBIT_RNC,ap);}
static void verror_handler_rnd(int erno,va_list ap) {rnl_verror_handler(erno|ERBIT_RND,ap);}

static int initialized=0;
void rnl_init(void) {
  if(!initialized) { initialized=1;
    rn_init();
    rnc_init(); rnc_verror_handler=&verror_handler_rnc;
    rnd_init(); rnd_verror_handler=&verror_handler_rnd;
  }
}

void rnl_clear(void) {}

static int load(struct rnc_source *sp) {
  int start=-1;
  if(!rnc_errors(sp)) start=rnc_parse(sp); rnc_close(sp);
  if(!rnc_errors(sp)&&(start=rnd_fixup(start))) {
    start=rn_compress_last(start);
  } else start=0;
  return start;
}

int rnl_fn(char *fn) {
  struct rnc_source src;
  rnc_open(&src,fn); return load(&src);
}

int rnl_fd(char *fn,int fd) {
  struct rnc_source src;
  rnc_bind(&src,fn,fd); return load(&src);
}

int rnl_s(char *fn,char *s,int len) {
  struct rnc_source src;
  rnc_stropen(&src,fn,s,len); return load(&src);
}
