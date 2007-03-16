/* $Id: u.c,v 1.24 2004/01/25 09:44:04 dvd Exp $ */

#include "u.h"

#define ux(u,c) (((u)<<6)|(c&0x3F))
#define u1(t) t[0]
#define u2(t) ux(t[0]&0x1F,t[1])
#define u3(t) ux(ux(t[0]&0xF,t[1]),t[2])
#define u4(t) ux(ux(ux(t[0]&0x7,t[1]),t[2]),t[3])
#define u5(t) ux(ux(ux(ux(t[0]&0x3,t[1]),t[2]),t[3]),t[4])
#define u6(t) ux(ux(ux(ux(ux(t[0]&0x1,t[1]),t[2]),t[3]),t[4]),t[5])

#define vx(c,u) c=0x80|((u)&0x3F)
#define v1(t,u) t[0]=u
#define v2(t,u) t[0]=0xC0|(u>>6);vx(t[1],u)
#define v3(t,u) t[0]=0xE0|(u>>12);vx(t[1],u>>6);vx(t[2],u)
#define v4(t,u) t[0]=0xF0|(u>>18);vx(t[1],u>>12);vx(t[2],u>>6);vx(t[3],u)
#define v5(t,u) t[0]=0xF8|(u>>24);vx(t[1],u>>18);vx(t[2],u>>12);vx(t[3],u>>6);vx(t[4],u)
#define v6(t,u) t[0]=0xFC|(u>>30);vx(t[1],u>>24);vx(t[2],u>>18);vx(t[3],u>>12);vx(t[4],u>>6);vx(t[5],u)

#define B1 0xFFFFFF80
#define B2 0xFFFFF800
#define B3 0xFFFF0000
#define B4 0xFFE00000
#define B5 0xFC000000
#define B6 0x80000000

#define BOM "\xEF\xBB\xBF"
#define BOMLEN 3

int u_bom(char *s,int n) {
  char *bom=(char*)(BOM+BOMLEN);
  if(n>=BOMLEN) {
    n=BOMLEN; s+=n;
    while(n--!=0) if(*(--s)!=*(--bom)) return 0;
    return BOMLEN;
  }
  return 0;
}

int u_get(int *up,char *s) {
  unsigned char *t=(unsigned char*)s;
  if(*t<0x80) {*up=u1(t); return 1;}
  if(*t<0xC0) return 0;
  if(*t<0xE0) {*up=u2(t); return (*up&B1)?2:0;}
  if(*t<0xF0) {*up=u3(t); return (*up&B2)?3:0;}
  if(*t<0xF8) {*up=u4(t); return (*up&B3)?4:0;}
  if(*t<0xFC) {*up=u5(t); return (*up&B4)?5:0;}
  if(*t<0xFE) {*up=u6(t); return (*up&B5)?6:0;}
  return 0;
}

int u_put(char *s,int u) {
  unsigned char *t=(unsigned char*)s;
  if(!(u&B1)) {v1(t,u); return 1;}
  if(!(u&B2)) {v2(t,u); return 2;}
  if(!(u&B3)) {v3(t,u); return 3;}
  if(!(u&B4)) {v4(t,u); return 4;}
  if(!(u&B5)) {v5(t,u); return 5;}
  if(!(u&B6)) {v6(t,u); return 6;}
  return 0;
}

int u_strlen(char *s) {int n=0; while(*(s+n)) ++n; return u_strnlen(s,n);}
int u_strnlen(char *s,int n) {
  int i,len=0,u;
  char *end=s+n;
  for(;;) {
    if(s==end) break;
    i=u_get(&u,s);
    if(i==0) {len=-1; break;}
    s+=i;
    if(s>end) {len=-1; break;}
    ++len;
  }
  return len;
}

int u_in_ranges(int u,int r[][2],int len) {
  int n=0,m=len-1,i;
  for(;;) {
    if(n>m) return 0;
    i=(n+m)/2;
    if(u<r[i][0]) m=i-1;
    else if(u>r[i][1]) n=i+1;
    else return 1;
  }
}
