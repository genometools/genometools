/* $Id: sc.h,v 1.7 2004/01/15 23:47:45 dvd Exp $ */

#ifndef SC_H
#define SC_H 1

#define SC_RECSIZE 3 /* 0 - key, 1 - value, 2 - auxiliary */

struct sc_stack {
  int (*tab)[SC_RECSIZE];
  int len,base,top;
};

extern void sc_init(struct sc_stack *stp);
extern void sc_clear(struct sc_stack *stp);

extern void sc_open(struct sc_stack *stp);
extern void sc_lock(struct sc_stack *stp);
extern void sc_close(struct sc_stack *stp);

extern int sc_void(struct sc_stack *sp);
extern int sc_locked(struct sc_stack *stp);

extern int sc_find(struct sc_stack *stp,int key); /* returns 0 if not found, index in tab otherwise */
extern int sc_add(struct sc_stack *stp,int key,int val,int aux); /* returns index for the new record */

#endif
