/* $Id: m.c,v 1.9 2004/03/13 13:28:02 dvd Exp $ */

#include <stdlib.h>
#include <string.h>
#include "er.h"
#include "m.h"

#ifndef M_STATIC
#define M_STATIC 0
#endif

#if M_STATIC

#ifndef M_FILL
#define M_FILL '\0'
#endif

static char memory[M_STATIC];
static char *mp=memory,*pmp=memory;

void m_free(void *p) {
  if(p==pmp) {
    mp=pmp; pmp=(char*)-1;
  }
}

void *m_alloc(int length,int size) {
  char *p=mp, *q=mp; int n=length*size;
  pmp=mp; mp+=(n+sizeof(int)-1)/sizeof(int)*sizeof(int);
  if(mp>=memory+M_STATIC) {
    (*er_printf)("failed to allocate %i bytes of memory\n",length*size);
    exit(1);
  }
  if(M_FILL!=-1) while(q!=mp) *(q++)=M_FILL;
  return (char*)p;
}

#else

void m_free(void *p) {
  free(p);
}

void *m_alloc(int length,int size) {
  void *p=malloc(length*size);
  if(p==NULL) {
    (*er_printf)("failed to allocate %i bytes of memory\n",length*size);
    exit(1);
  }
  return p;
}

#endif

void *m_stretch(void *p,int newlen,int oldlen,int size) {
  void *newp=m_alloc(newlen,size);
  memcpy(newp,p,oldlen*size);
  m_free(p);
  return newp;
}
