/* $Id: xsd_tm.h,v 1.2 2003/12/25 03:05:34 dvd Exp $ */

#ifndef XSD_TM_H
#define XSD_TM_H 1

struct xsd_tm {int days,secs,mics,tz;};

/* fmt is a combination of ymdtz */
extern void xsd_mktm(struct xsd_tm *tmp,char *fmt,char *val);
extern void xsd_mktmn(struct xsd_tm *tmp,char *fmt,char *s,int n);

/* -1 - less, 0 - equal, 1 - greater, other - unknown */
extern int xsd_tmcmp(struct xsd_tm *tmp1, struct xsd_tm *tmp2);

#endif
