/* $Id: rnx.h,v 1.7 2004/02/18 12:53:42 dvd Exp $ */

#ifndef RNX_H
#define RNX_H 1

extern void rnx_init(void);
extern void rnx_clear(void);

extern int rnx_n_exp,*rnx_exp;
extern void rnx_expected(int p,int req);

extern char *rnx_p2str(int p);
extern char *rnx_nc2str(int nc);

#endif
