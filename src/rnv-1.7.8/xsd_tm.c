/* $Id: xsd_tm.c,v 1.9 2004/02/25 00:00:32 dvd Exp $ */

#include <stdlib.h> /*strtol*/
#include <limits.h>
#include <string.h> /*strlen*/
#include <assert.h>
#include "xsd_tm.h"

static int leap(int yr) {return !(yr%4)&&((yr%100)||!(yr%400));}
static int y2d(int yr) {return yr*365+yr/4-yr/100+yr/400;}
static int ymd2dn(int yr,int mo,int dy) {
  switch(mo) {
  case 12: dy+=30;
  case 11: dy+=31;
  case 10: dy+=30;
  case 9: dy+=31;
  case 8: dy+=31;
  case 7: dy+=30;
  case 6: dy+=31;
  case 5: dy+=30;
  case 4: dy+=31;
  case 3: dy+=28;
  case 2: dy+=31;
  case 1: break;
  }
  if(mo>2&&leap(yr)) ++dy;
  return dy;
}

static int ymd2ds(int yr,int mo,int dy) {
  return (yr>=0?y2d(yr-1):y2d(yr)-366)+ymd2dn(yr,mo,dy);
}

#define DAYSECS 86400
#define TZSECS 50400

static void addsecs(struct xsd_tm *tmp,int secs) {
  tmp->secs+=secs;
  if(tmp->secs<0) {
    --tmp->days;
    tmp->secs+=DAYSECS;
  } else if(tmp->secs>=DAYSECS) {
    ++tmp->days;
    tmp->secs-=DAYSECS;
  }
}

void xsd_mktmn(struct xsd_tm *tmp,char *fmt,char *s,int n) {
  char *end=s+n;
  int yr=2000,mo=1,dy=1,hr=0,mi=0,zh=15,zm=0;
  double se=0.0;
  for(;;) {
    if(s==end||!*fmt) break;
    switch(*s) {
    case '-':
      switch(*fmt) {
      case 'y': ++fmt; yr=strtol(s,&s,10); continue;
      case 'z': ++fmt; ++s; zh=strtol(s,&s,10); ++s; zm=strtol(s,&s,10); continue;
      }
      break;
    case '+': assert(*fmt=='z'); ++fmt;
	  zh=-strtol(s,&s,10); ++s; zm=-strtol(s,&s,10); continue;
    case 'Z': assert(*fmt=='z'); ++fmt; zh=0; zm=0; ++s; continue;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
      switch(*(fmt++)) {
      case 'y': yr=strtol(s,&s,10); continue;
      case 'm': mo=strtol(s,&s,10); continue;
      case 'd': dy=strtol(s,&s,10); continue;
      case 't': hr=strtol(s,&s,10); ++s; mi=strtol(s,&s,10); ++s; se=strtod(s,&s); continue;
      }
      break;
    }
    ++s;
  }
  tmp->mics=(int)((se-(int)se)*1000000+0.5);
  tmp->secs=(int)se+60*(mi+60*hr);
  tmp->days=ymd2ds(yr,mo,dy);
  if((tmp->tz=(zh!=15))) addsecs(tmp,60*(zm+60*zh));
}
void xsd_mktm(struct xsd_tm *tmp,char *fmt,char *val) {xsd_mktmn(tmp,fmt,val,strlen(val));}

static int tmcmp(struct xsd_tm *tmp1, struct xsd_tm *tmp2) {
  int dd=tmp1->days-tmp2->days, ds=tmp1->secs-tmp2->secs, dm=tmp1->mics-tmp2->mics;
  return dd<0?-1:dd>0?1:ds<0?-1:ds>0?1:dm<0?-1:dm>0?1:0;
}

extern int xsd_tmcmp(struct xsd_tm *tmp1, struct xsd_tm *tmp2) {
  if(tmp1->tz==tmp2->tz) {
    return tmcmp(tmp1,tmp2);
  } else if(tmp1->tz) {
    struct xsd_tm tm; tm.mics=tmp2->mics;
    tm.days=tmp2->days; tm.secs=tmp2->secs; addsecs(&tm,TZSECS);
    if(tmcmp(tmp1,&tm)==1) return 1;
    tm.days=tmp2->days; tm.secs=tmp2->secs; addsecs(&tm,-TZSECS);
    if(tmcmp(tmp1,&tm)==-1) return -1;
    return 2;
  } else return -xsd_tmcmp(tmp2,tmp1);
}

