/* $Id: m.h 309 2004-01-19 20:33:42Z dvd $ */

#ifndef M_H
#define M_H 1

extern void m_free(void *p);
extern void *m_alloc(int length,int size);
extern void *m_stretch(void *p,int newlen,int oldlen,int size);

#endif
