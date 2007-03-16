/* $Id: xsd.c,v 1.47 2005/01/06 21:04:06 dvd Exp $ */

#include <limits.h> /*INT_MAX*/
#include <stdlib.h> /*atof,atol,strtol*/
#include <string.h> /*strlen*/
#include <math.h> /*HUGE_VAL*/
#include <assert.h>
#include "u.h"
#include "xmlc.h"
#include "s.h"
#include "erbit.h"
#include "rx.h"
#include "xsd_tm.h"
#include "er.h"
#include "xsd.h"

#define err(msg) (*er_vprintf)(msg"\n",ap)
void xsd_default_verror_handler(int erno,va_list ap) {
  (*er_printf)("XML Schema datatypes: ");
  if(erno&ERBIT_RX) {
    rx_default_verror_handler(erno&~ERBIT_RX,ap);
  } else {
    switch(erno) {
    case XSD_ER_TYP: err("unknown type %s"); break;
    case XSD_ER_PAR: err("unknown parameter %s"); break;
    case XSD_ER_PARVAL: err("invalid parameter value %s=\"%s\""); break;
    case XSD_ER_VAL: err("invalid typed value \"%s\" for type %s"); break;
    case XSD_ER_NPAT: err("no more than 16 patterns per type are supported"); break;
    case XSD_ER_WS: err("the builtin derived datatype that specifies the desired value for the whiteSpace facet should be used instead of 'whiteSpace'"); break;
    case XSD_ER_ENUM: err("'value' should be used instead of 'enumeration'"); break;
    default: assert(0);
    }
  }
}

void (*xsd_verror_handler)(int erno,va_list ap)=&xsd_default_verror_handler;

static void error_handler(int erno,...) {
  va_list ap; va_start(ap,erno); (*xsd_verror_handler)(erno,ap); va_end(ap);
}

static void verror_handler_rx(int erno,va_list ap) {(*xsd_verror_handler)(erno|ERBIT_RX,ap);}

static void windup(void);
static int initialized=0;
void xsd_init(void) {
  if(!initialized) { initialized=1;
    rx_init(); rx_verror_handler=&verror_handler_rx;
    windup();
  }
}

void xsd_clear(void) {
  windup();
}

static void windup(void) {
}

#define FCT_ENUMERATION 0
#define FCT_FRACTION_DIGITS 1
#define FCT_LENGTH 2
#define FCT_MAX_EXCLUSIVE 3
#define FCT_MAX_INCLUSIVE 4
#define FCT_MAX_LENGTH 5
#define FCT_MIN_EXCLUSIVE 6
#define FCT_MIN_INCLUSIVE 7
#define FCT_MIN_LENGTH 8
#define FCT_PATTERN 9
#define FCT_TOTAL_DIGITS 10
#define FCT_WHITE_SPACE 11
#define NFCT 12
static char *fcttab[NFCT]={
  "enumeration", "fractionDigits", "length", "maxExclusive", "maxInclusive", "maxLength",
  "minExclusive", "minInclusive", "minLength", "pattern", "totalDigits", "whiteSpace"};

#define FCT_IBOUNDS (1<<FCT_MIN_INCLUSIVE|1<<FCT_MAX_INCLUSIVE)
#define FCT_EBOUNDS (1<<FCT_MIN_EXCLUSIVE|1<<FCT_MAX_EXCLUSIVE)
#define FCT_BOUNDS (FCT_IBOUNDS|FCT_EBOUNDS)

#define WS_PRESERVE 0
#define WS_REPLACE 1
#define WS_COLLAPSE 2

static int (*match[])(char *r,char *s,int n)={&rx_match,&rx_rmatch,&rx_cmatch};

#define TYP_ENTITIES 0
#define TYP_ENTITY 1
#define TYP_ID 2
#define TYP_IDREF 3
#define TYP_IDREFS 4
#define TYP_NCNAME 5
#define TYP_NMTOKEN 6
#define TYP_NMTOKENS 7
#define TYP_NOTATION 8
#define TYP_NAME 9
#define TYP_QNAME 10
#define TYP_ANY_URI 11
#define TYP_BASE64_BINARY 12
#define TYP_BOOLEAN 13
#define TYP_BYTE 14
#define TYP_DATE 15
#define TYP_DATE_TIME 16
#define TYP_DECIMAL 17
#define TYP_DOUBLE 18
#define TYP_DURATION 19
#define TYP_FLOAT 20
#define TYP_G_DAY 21
#define TYP_G_MONTH 22
#define TYP_G_MONTH_DAY 23
#define TYP_G_YEAR 24
#define TYP_G_YEAR_MONTH 25
#define TYP_HEX_BINARY 26
#define TYP_INT 27
#define TYP_INTEGER 28
#define TYP_LANGUAGE 29
#define TYP_LONG 30
#define TYP_NEGATIVE_INTEGER 31
#define TYP_NON_NEGATIVE_INTEGER 32
#define TYP_NON_POSITIVE_INTEGER 33
#define TYP_NORMALIZED_STRING 34
#define TYP_POSITIVE_INTEGER 35
#define TYP_SHORT 36
#define TYP_STRING 37
#define TYP_TIME 38
#define TYP_TOKEN 39
#define TYP_UNSIGNED_BYTE 40
#define TYP_UNSIGNED_INT 41
#define TYP_UNSIGNED_LONG 42
#define TYP_UNSIGNED_SHORT 43
#define NTYP 44
static char *typtab[NTYP]={
"ENTITIES", "ENTITY", "ID", "IDREF", "IDREFS", "NCName", "NMTOKEN", "NMTOKENS",
"NOTATION", "Name", "QName", "anyURI", "base64Binary", "boolean", "byte", "date",
"dateTime", "decimal", "double", "duration", "float", "gDay", "gMonth",
"gMonthDay", "gYear", "gYearMonth", "hexBinary", "int", "integer", "language",
"long", "negativeInteger", "nonNegativeInteger", "nonPositiveInteger",
"normalizedString", "positiveInteger", "short", "string", "time", "token",
"unsignedByte", "unsignedInt", "unsignedLong", "unsignedShort"};

#define ERR_PARAMETER "invalid XML Schema datatype parameter '%s'"
#define ERR_DATATYPE "invalid XML Schema datatype name '%s'"
#define ERR_VALUE "invalid value '%s' for XML Schema datatype '%s'"

struct dura {int yr,mo,dy,hr,mi;double se;};
static void durainit(struct dura *d) {d->yr=d->mo=d->dy=d->hr=d->mi=0; d->se=0.0;}

static void s2dura(struct dura *dp,char *s,int n) {
  char *end=s+n,*np="0";
  int sign=1,time=0;
  durainit(dp);
  while(s!=end) {
    switch(*s) {
    case '-': sign=-1; break;
    case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
    case '8': case '9': case '.': np=s; break;
    case 'T': time=1; break;
    case 'Y': dp->yr=sign*atoi(np); break;
    case 'M': if(time) dp->mi=sign*atoi(np); else dp->mo=sign*atoi(np); break;
    case 'D': dp->dy=sign*atoi(np); break;
    case 'H': dp->hr=sign*atoi(np); break;
    case 'S': dp->se=sign*atof(np); break;
    }
    ++s;
  }
}

static int duracmp(char *s1,char *s2,int n) {
  struct dura d1,d2;
  s2dura(&d1,s1,strlen(s1)); s2dura(&d2,s2,n);
  if(d1.yr!=d2.yr) return d1.yr-d2.yr;
  if(d1.mo!=d2.mo) return d1.mo-d2.mo;
  if(d1.dy!=d2.dy) return d1.dy-d2.dy;
  if(d1.hr!=d2.hr) return d1.hr-d2.hr;
  if(d1.mi!=d2.mi) return d1.mi-d2.mi;
  if(d1.se!=d2.se) return d1.se<d2.se?-1:1;
  return 0;
}

static int dtcmpn(char *s1,char *s2,int n,char *fmt) {
  struct xsd_tm tm1,tm2;
  xsd_mktm(&tm1,fmt,s1); xsd_mktmn(&tm2,fmt,s2,n);
  return xsd_tmcmp(&tm1,&tm2);
}

static int toklenn(char *s,int n) {
  char *end=s+n;
  int u,len=0;
  SKIP_SPACE:
  for(;;) { if(s==end) return len?len-1:0;
    s+=u_get(&u,s);
    if(!xmlc_white_space(u)) break;
  }
  ++len;
  for(;;) { if(s==end) return len;
    ++len;
    s+=u_get(&u,s);
    if(xmlc_white_space(u)) goto SKIP_SPACE;
  }
}

static int tokcntn(char *s,int n) {
  char *end=s+n;
  int u,cnt=0;
  SKIP_SPACE:
  for(;;) { if(s==end) return cnt;
    s+=u_get(&u,s);
    if(!xmlc_white_space(u)) break;
  }
  ++cnt;
  for(;;) { if(s==end) return cnt;
    s+=u_get(&u,s);
    if(xmlc_white_space(u)) goto SKIP_SPACE;
  }
}

static int b64lenn(char *s,int n) {
  char *end=s+n;
  int l=0,len;
  for(;;) { if(end==s) break;
    --end;
    if(!xmlc_white_space(*end)&&*end!='=') {++end; break;}
  }
  while(s!=end) {if(!xmlc_white_space(*s)) ++l; ++s;}
  len=l/4*3;
  switch(l%4) {
  case 0: break;
  case 1: len=-1; break;
  case 2: len+=1; break;
  case 3: len+=2; break;
  }
  return len;
}

static int fdiglenn(char *s,int n) {
  char *end=s+n; int len=0;
  for(;;) { if(end==s) break;
    --end;
    if(*end!='0'&&!xmlc_white_space(*end)) {++end; break;}
  }
  for(;;) { if(s==end) break;
    if(*(s++)=='.') {
      while(s++!=end) ++len;
      break;
    }
  }
  return len;
}

static int diglenn(char *s,int n) {
  char *end=s+n; int len=0;
  for(;;) { if(s==end) break;
    if(!(xmlc_white_space(*s)||*s=='+'||*s=='-'||*s=='0')) break;
    ++s;
  }
  for(;;) { if(s==end||*s=='.'||xmlc_white_space(*s)) break;
    ++len; ++s;
  }
  if(len==0) len=1;
  if(*s=='.') len+=fdiglenn(s,end-s);
  return len;
}

#define NPAT 16

struct facets {
  int set;
  char *pattern[NPAT+1]; int npat;
  int length, minLength, maxLength, totalDigits, fractionDigits;
  char *maxExclusive, *maxInclusive, *minExclusive, *minInclusive;
  int whiteSpace;
};

/* PAT_DECIMAL is unsigned decimal, signed decimal matches PAT_FIXED */
#define PAT_ORDINAL "([0-9]+)"
#define PAT_FRACTIONAL "(\\.[0-9]+)"
#define PAT_DECIMAL "("PAT_ORDINAL"\\.?|"PAT_ORDINAL"?"PAT_FRACTIONAL")"

#define PAT_POSITIVE "\\+?"PAT_ORDINAL
#define PAT_NON_NEGATIVE "\\+?"PAT_ORDINAL
#define PAT_NON_POSITIVE "\\-"PAT_ORDINAL"|0+"
#define PAT_NEGATIVE "\\-"PAT_ORDINAL
#define PAT_INTEGER "([+\\-]?"PAT_ORDINAL")"

#define PAT_FIXED "([+\\-]?"PAT_DECIMAL")"
#define PAT_FLOATING PAT_FIXED"([Ee]"PAT_INTEGER")?|INF|-INF|NaN"

#define PAT_HEX_BINARY "[0-9a-fA-F]+"

#define PAT_BASE64 "([A-Za-z0-9+/] ?)"
#define PAT_BASE64_2 "([AQgw] ?)"
#define PAT_BASE64_1 "([AEIMQUYcgkosw048] ?)"
#define PAT_BASE64_BINARY \
   "("PAT_BASE64"{4})*" \
   "("PAT_BASE64 PAT_BASE64_2"= ?=" \
   "|"PAT_BASE64"{2}" PAT_BASE64_1"=)?"

#define PAT_ANY_URI "(([a-zA-Z][0-9a-zA-Z+\\-\\.]*:)?/{0,2}[0-9a-zA-Z;/?:@&=+$\\.\\-_!~*'()%]+)?(#[0-9a-zA-Z;/?:@&=+$\\.\\-_!~*'()%]+)?"

#define PAT_NCNAME "[\\i-[:]][\\c-[:]]*"
#define PAT_QNAME "("PAT_NCNAME":)?"PAT_NCNAME
#define PAT_NMTOKEN "\\c+"
#define PAT_NAME "\\i\\c*"
#define PAT_NCNAMES PAT_NCNAME"( "PAT_NCNAME")*"
#define PAT_NMTOKENS PAT_NMTOKEN"( "PAT_NMTOKEN")*"
#define PAT_NAMES PAT_NAME"( "PAT_NAME")*"

#define PAT_LANGUAGE "([a-zA-Z]{1,8}(-[a-zA-Z0-9]{1,8})*)"

#define PAT_DURAY "("PAT_ORDINAL"Y)"
#define PAT_DURAM "("PAT_ORDINAL"M)"
#define PAT_DURAD "("PAT_ORDINAL"D)"
#define PAT_DURADATE \
  "(" PAT_DURAY   PAT_DURAM"?"PAT_DURAD"?" \
  "|" PAT_DURAY"?"PAT_DURAM   PAT_DURAD"?" \
  "|" PAT_DURAY"?"PAT_DURAM"?"PAT_DURAD ")"
#define PAT_DURAH "("PAT_ORDINAL"H)"
#define PAT_DURAM "("PAT_ORDINAL"M)"
#define PAT_DURAS "("PAT_DECIMAL"S)"
#define PAT_DURATIME \
"(T(" PAT_DURAH   PAT_DURAM"?"PAT_DURAS"?" \
  "|" PAT_DURAM"?"PAT_DURAM   PAT_DURAS"?" \
  "|" PAT_DURAS"?"PAT_DURAM"?"PAT_DURAS "))"
#define PAT_DURATION "-?P("PAT_DURADATE PAT_DURATIME"|"PAT_DURADATE"|"PAT_DURATIME")"

#define PAT_ZONE "(Z|[+\\-](0[0-9]|1[0-4]):[0-5][0-9])"
#define PAT_YEAR0 "[0-9]{4,}"
#define PAT_MONTH0 "(0[1-9]|1[0-2])"
#define PAT_DAY0 "([0-2][0-9]|3[01])"
#define PAT_YEAR "-?"PAT_YEAR0 PAT_ZONE"?"
#define PAT_MONTH "--"PAT_MONTH0"--"PAT_ZONE"?"
#define PAT_DAY "---"PAT_DAY0 PAT_ZONE"?"
#define PAT_YEAR_MONTH "-?"PAT_YEAR0"-"PAT_MONTH0 PAT_ZONE"?"
#define PAT_MONTH_DAY "--"PAT_MONTH0"-"PAT_DAY0 PAT_ZONE"?"
#define PAT_DATE0 PAT_YEAR0"-"PAT_MONTH0"-"PAT_DAY0
#define PAT_TIME0 "([0-1][0-9]|2[0-3]):[0-5][0-9]:([0-5][0-9]|60)"PAT_FRACTIONAL"?"
#define PAT_DATE "-?"PAT_DATE0 PAT_ZONE"?"
#define PAT_TIME PAT_TIME0 PAT_ZONE"?"
#define PAT_DATE_TIME "-?"PAT_DATE0"T"PAT_TIME0 PAT_ZONE"?"

static void anchdec(int *plus,int *zero,char **beg,char **dp,char **end,char *s,int n) {
  char *end0=s+n;
  *beg=s; *zero=1; *plus=1;
  for(;;) { if(end0==*beg) break;
    --end0;
    if(!xmlc_white_space(*end0)) {++end0; break;}
  }
  *end=end0;
  for(;;) { if(*end==*beg) break;
    --*end;
    if(!(**end=='0'||**end=='+'||**end=='-')) {
      if(**end!='.') *zero=0; 
      ++*end; 
      break;
    }
  }
  *dp=*end;
  for(;;) { if(*beg==*end) break;
    if(**beg=='-') *plus=0;
    else if(!(**beg=='0'||**beg=='+'||xmlc_white_space(**beg))) {
      if(**beg!='.') *zero=0;
      for(;;) {
	if(*dp==*beg) {*dp=*end=end0; break;}
	--*dp;
	if(**dp=='.') break;
      }
      break;
    }
    ++*beg;
  }
}

static int deccmp(char *s1,int n1,char *s2,int n2) {
  int p1,p2,z1,z2,cmp;
  char *d1,*e1,*d2,*e2,*c1,*c2;
  anchdec(&p1,&z1,&s1,&d1,&e1,s1,n1); anchdec(&p2,&z2,&s2,&d2,&e2,s2,n2);
  if(z1&&z2) return 0;
  if(p1!=p2) return p1-p2;
  cmp=0;
  if(d1-s1!=d2-s2) cmp=d1-s1-(d2-s2);
  if(cmp!=0) return p1?cmp:-cmp;
  c1=s1; c2=s2;
  for(;;) {
    if(c1==d1) break;
    if(*c1!=*c2) {cmp=*c1-*c2; break;}
    ++c1; ++c2;
  }
  if(cmp!=0) return p1?cmp:-cmp;
  if(c1!=e1) ++c1; if(c2!=e2) ++c2;
  for(;;) {
    if(c1==e1) {cmp=-(c2!=e2); break;}
    if(c2==e2) {cmp=1; break;}
    if(*c1!=*c2) {cmp=*c1-*c2; break;}
    ++c1; ++c2;
  }
  return p1?cmp:-cmp;
}

static int chkdec(struct facets *fp,char *s,int n) {
  int ok=1;
  if(fp->set&(1<<FCT_MIN_EXCLUSIVE)) ok=ok&&deccmp(s,n,fp->minExclusive,strlen(fp->minExclusive))>0;
  if(fp->set&(1<<FCT_MIN_INCLUSIVE)) ok=ok&&deccmp(s,n,fp->minInclusive,strlen(fp->minInclusive))>=0;
  if(fp->set&(1<<FCT_MAX_INCLUSIVE)) ok=ok&&deccmp(s,n,fp->maxInclusive,strlen(fp->maxInclusive))<=0;
  if(fp->set&(1<<FCT_MAX_EXCLUSIVE)) ok=ok&&deccmp(s,n,fp->maxExclusive,strlen(fp->maxExclusive))<0;
  return ok;
}

static double atodn(char *s,int n) {
  return s_tokcmpn("-INF",s,n)==0?-HUGE_VAL
    : s_tokcmpn("INF",s,n)==0?HUGE_VAL
    : atof(s);
}
static double atod(char *s) {return atodn(s,strlen(s));}

static int chkdbl(struct facets *fp,char *s,int n) {
  int ok=1,nan=s_tokcmpn("NaN",s,n)==0;
  double d=atodn(s,n);
  if(fp->set&(1<<FCT_MIN_EXCLUSIVE)) ok=ok&&!nan&&d>atod(fp->minExclusive);
  if(fp->set&(1<<FCT_MIN_INCLUSIVE)) ok=ok&&!nan&&d>=atod(fp->minInclusive);
  if(fp->set&(1<<FCT_MAX_INCLUSIVE)) ok=ok&&!nan&&d<=atod(fp->maxInclusive);
  if(fp->set&(1<<FCT_MAX_EXCLUSIVE)) ok=ok&&!nan&&d<atod(fp->maxExclusive);
  return ok;
}

static int chktmlim(char *typ,char *fmt,char *val,int cmpmin,int cmpmax,struct xsd_tm *tmp) {
  struct xsd_tm tmf; int cmp;
  if(!xsd_allows(typ,"",val,strlen(val))) {(*error_handler)(XSD_ER_PARVAL); return 0;}
  xsd_mktm(&tmf,fmt,val);
  cmp=xsd_tmcmp(tmp,&tmf);
  return cmpmin<=cmp&&cmp<=cmpmax;
}

static int chktm(char *typ,char *fmt,struct facets *fp,char *s,int n) {
  int ok=1;
  struct xsd_tm tms;
  if(!xsd_allows(typ,"",s,n)) return 0;
  xsd_mktmn(&tms,fmt,s,n);
  if(fp->set&(1<<FCT_MIN_EXCLUSIVE)) ok=ok&&chktmlim(typ,fmt,fp->minExclusive,1,1,&tms);
  if(fp->set&(1<<FCT_MIN_INCLUSIVE)) ok=ok&&chktmlim(typ,fmt,fp->minInclusive,0,1,&tms);
  if(fp->set&(1<<FCT_MAX_INCLUSIVE)) ok=ok&&chktmlim(typ,fmt,fp->maxInclusive,-1,0,&tms);
  if(fp->set&(1<<FCT_MAX_EXCLUSIVE)) ok=ok&&chktmlim(typ,fmt,fp->maxExclusive,-1,-1,&tms);
  return ok;
}

int xsd_allows(char *typ,char *ps,char *s,int n) {
  int ok=1,length;
  int dt=s_tab(typ,typtab,NTYP);
  struct facets fct; fct.set=0; fct.npat=0;
  switch(dt) {
  case TYP_INTEGER:
    fct.pattern[fct.npat++]=PAT_INTEGER;
    dt=TYP_DECIMAL;
    break;
  case TYP_POSITIVE_INTEGER:
    fct.pattern[fct.npat++]=PAT_POSITIVE;
    dt=TYP_DECIMAL; fct.set|=1<<FCT_MIN_INCLUSIVE;
    fct.minInclusive="1";
    break;
  case TYP_NON_NEGATIVE_INTEGER:
    fct.pattern[fct.npat++]=PAT_NON_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=1<<FCT_MIN_INCLUSIVE;
    fct.minInclusive="0";
    break;
  case TYP_NON_POSITIVE_INTEGER:
    fct.pattern[fct.npat++]=PAT_NON_POSITIVE;
    dt=TYP_DECIMAL; fct.set|=1<<FCT_MAX_INCLUSIVE;
    fct.maxInclusive="0";
    break;
  case TYP_NEGATIVE_INTEGER:
    fct.pattern[fct.npat++]=PAT_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=1<<FCT_MAX_INCLUSIVE;
    fct.maxInclusive="-1";
    break;
  case TYP_BYTE:
    fct.pattern[fct.npat++]=PAT_INTEGER;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="-128"; fct.maxInclusive="127";
    break;
  case TYP_UNSIGNED_BYTE:
    fct.pattern[fct.npat++]=PAT_NON_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="0"; fct.maxInclusive="255";
    break;
  case TYP_SHORT:
    fct.pattern[fct.npat++]=PAT_INTEGER;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="-32768"; fct.maxInclusive="32767";
    break;
  case TYP_UNSIGNED_SHORT:
    fct.pattern[fct.npat++]=PAT_NON_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="0"; fct.maxInclusive="65535";
    break;
  case TYP_INT:
    fct.pattern[fct.npat++]=PAT_INTEGER;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="-2147483648"; fct.maxInclusive="2147483647";
    break;
  case TYP_UNSIGNED_INT:
    fct.pattern[fct.npat++]=PAT_NON_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="0"; fct.maxInclusive="4294967295";
    break;
  case TYP_LONG:
    fct.pattern[fct.npat++]=PAT_INTEGER;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="-9223372036854775808"; fct.maxInclusive="9223372036854775807";
    break;
  case TYP_UNSIGNED_LONG:
    fct.pattern[fct.npat++]=PAT_NON_NEGATIVE;
    dt=TYP_DECIMAL; fct.set|=FCT_IBOUNDS;
    fct.minInclusive="0"; fct.maxInclusive="18446744073709551615";
    break;
  }

  { int n;
    while((n=strlen(ps))) {
      char *key=ps,*val=key+n+1,*end,i;
      switch(i=s_tab(key,fcttab,NFCT)) {
      case FCT_LENGTH: fct.length=(int)strtol(val,&end,10); if(!*val||*end) (*error_handler)(XSD_ER_PARVAL,key,val); break;
      case FCT_MAX_LENGTH: fct.maxLength=(int)strtol(val,&end,10); if(!*val||*end) (*error_handler)(XSD_ER_PARVAL,key,val); break;
      case FCT_MIN_LENGTH: fct.minLength=(int)strtol(val,&end,10); if(!*val||*end) (*error_handler)(XSD_ER_PARVAL,key,val); break;
      case FCT_FRACTION_DIGITS: fct.fractionDigits=(int)strtol(val,&end,10); if(!*val||*end) (*error_handler)(XSD_ER_PARVAL,key,val); break;
      case FCT_TOTAL_DIGITS: fct.totalDigits=(int)strtol(val,&end,10); if(!*val||*end) (*error_handler)(XSD_ER_PARVAL,key,val); break;
      case FCT_PATTERN:
	if(fct.npat==NPAT) (*error_handler)(XSD_ER_NPAT); else {
	  fct.pattern[fct.npat++]=val;
	} break;
      case FCT_MAX_EXCLUSIVE: fct.maxExclusive=val; break;
      case FCT_MAX_INCLUSIVE: fct.maxInclusive=val; break;
      case FCT_MIN_EXCLUSIVE: fct.minExclusive=val; break;
      case FCT_MIN_INCLUSIVE: fct.minInclusive=val; break;
      case FCT_WHITE_SPACE: (*error_handler)(XSD_ER_WS); break;
      case FCT_ENUMERATION: (*error_handler)(XSD_ER_ENUM); break;
      case NFCT: (*error_handler)(XSD_ER_PAR,key); break;
      default: assert(0);
      }
      fct.set|=1<<i;
      ps=val+strlen(val)+1;
    }
  }

  fct.whiteSpace=WS_COLLAPSE;
  length=INT_MAX;
  switch(dt) {
 /*primitive*/
  case TYP_STRING: fct.whiteSpace=WS_PRESERVE;
    length=u_strnlen(s,n);
    break;
  case TYP_BOOLEAN:
    fct.pattern[fct.npat++]="true|false|1|0";
    break;
  case TYP_DECIMAL:
    fct.pattern[fct.npat++]=PAT_FIXED;
    if(fct.set&(1<<FCT_FRACTION_DIGITS)) ok=ok&&fdiglenn(s,n)<=fct.fractionDigits;
    if(fct.set&(1<<FCT_TOTAL_DIGITS)) ok=ok&&diglenn(s,n)<=fct.totalDigits;
    if(fct.set&FCT_BOUNDS) ok=ok&chkdec(&fct,s,n);
    break;
  case TYP_FLOAT: case TYP_DOUBLE: /* float and double is the same type */
    fct.pattern[fct.npat++]=PAT_FLOATING;
    if(fct.set&FCT_BOUNDS) ok=ok&chkdbl(&fct,s,n);
    break;
  case TYP_DURATION:
    fct.pattern[fct.npat++]=PAT_DURATION;
    break;
  case TYP_DATE_TIME:
    fct.pattern[fct.npat++]=PAT_DATE_TIME;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"ymdtz",&fct,s,n);
    break;
  case TYP_DATE:
    fct.pattern[fct.npat++]=PAT_DATE;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"ymdz",&fct,s,n);
    break;
  case TYP_TIME:
    fct.pattern[fct.npat++]=PAT_TIME;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"tz",&fct,s,n);
    break;
  case TYP_G_YEAR_MONTH:
    fct.pattern[fct.npat++]=PAT_YEAR_MONTH;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"ymz",&fct,s,n);
    break;
  case TYP_G_YEAR:
    fct.pattern[fct.npat++]=PAT_YEAR;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"yz",&fct,s,n);
    break;
  case TYP_G_MONTH_DAY:
    fct.pattern[fct.npat++]=PAT_MONTH_DAY;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"mdz",&fct,s,n);
    break;
  case TYP_G_DAY:
    fct.pattern[fct.npat++]=PAT_DAY;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"dz",&fct,s,n);
    break;
  case TYP_G_MONTH:
    fct.pattern[fct.npat++]=PAT_MONTH;
    if(fct.set&FCT_BOUNDS) ok=ok&chktm(typ,"mz",&fct,s,n);
    break;
  case TYP_HEX_BINARY:
    fct.pattern[fct.npat++]=PAT_HEX_BINARY;
    length=(toklenn(s,n)+1)/2;
    break;
  case TYP_BASE64_BINARY:
    fct.pattern[fct.npat++]=PAT_BASE64_BINARY;
    length=b64lenn(s,n);
    break;
  case TYP_ANY_URI:
    fct.pattern[fct.npat++]=PAT_ANY_URI;
    length=toklenn(s,n);
    break;
  case TYP_QNAME: case TYP_NOTATION:
    fct.pattern[fct.npat++]=PAT_QNAME;
    fct.set&=~(1<<FCT_LENGTH|1<<FCT_MIN_LENGTH|1<<FCT_MAX_LENGTH); /* the errata states that any value is valid */
    break;
 /*derived*/
  case TYP_NORMALIZED_STRING: fct.whiteSpace=WS_REPLACE;
    length=u_strnlen(s,n);
    break;
  case TYP_TOKEN:
    length=toklenn(s,n);
    break;
  case TYP_LANGUAGE:
    fct.pattern[fct.npat++]=PAT_LANGUAGE;
    length=toklenn(s,n);
    break;
  case TYP_NMTOKEN:
    fct.pattern[fct.npat++]=PAT_NMTOKEN;
    length=toklenn(s,n);
    break;
  case TYP_NMTOKENS:
    fct.pattern[fct.npat++]=PAT_NMTOKENS;
    length=tokcntn(s,n);
    break;
  case TYP_NAME:
    fct.pattern[fct.npat++]=PAT_NAME;
    length=toklenn(s,n);
    break;
  case TYP_NCNAME:
    fct.pattern[fct.npat++]=PAT_NCNAME;
    length=toklenn(s,n);
    break;
  case TYP_ID:
    fct.pattern[fct.npat++]=PAT_NCNAME;
    length=toklenn(s,n);
    break;
  case TYP_IDREF:
    fct.pattern[fct.npat++]=PAT_NCNAME;
    length=toklenn(s,n);
    break;
  case TYP_IDREFS:
    fct.pattern[fct.npat++]=PAT_NCNAMES;
    length=tokcntn(s,n);
    break;
  case TYP_ENTITY:
    fct.pattern[fct.npat++]=PAT_NCNAME;
    length=toklenn(s,n);
    break;
  case TYP_ENTITIES:
    fct.pattern[fct.npat++]=PAT_NCNAMES;
    length=tokcntn(s,n);
    break;
  case NTYP: (*error_handler)(XSD_ER_TYP,typ); break;
  default: assert(0);
  }

  while(fct.npat--) ok=ok&&match[fct.whiteSpace](fct.pattern[fct.npat],s,n);

  if(fct.set&(1<<FCT_LENGTH)) ok=ok&&length==fct.length;
  if(fct.set&(1<<FCT_MAX_LENGTH)) ok=ok&&length<=fct.maxLength;
  if(fct.set&(1<<FCT_MIN_LENGTH)) ok=ok&&length>=fct.minLength;

  return ok;
}

static int dblcmpn(char *val,char *s,char n) {
  double d1,d2;
  return s_tokcmpn(val,s,n)==0?0
    : s_tokcmpn(val,"NaN",3)==0||s_tokcmpn("NaN",s,n)==0?1
    : (d1=atod(val),d2=atodn(s,n),d1<d2?-1:d1>d2?1:0);
}

static int hexcmpn(char *s1,char *s2,int n) {
  char *end=s2+n;
  for(;;++s1,++s2) {
    while(*s1&&xmlc_white_space(*s1)) ++s1;
    while(s2!=end&&xmlc_white_space(*s2)) ++s2;
    if(s2==end) return *s1;
    if(!*s1) return -*s2;
    switch(*s1) {
    case 'a': case 'A': if(*s2=='a'||*s2=='A') continue;
    case 'b': case 'B': if(*s2=='b'||*s2=='B') continue;
    case 'c': case 'C': if(*s2=='c'||*s2=='C') continue;
    case 'd': case 'D': if(*s2=='d'||*s2=='D') continue;
    case 'e': case 'E': if(*s2=='e'||*s2=='E') continue;
    case 'f': case 'F': if(*s2=='f'||*s2=='F') continue;
    default: if(*s1!=*s2) return *s1-*s2;
    }
  }
}

static int b64cmpn(char *s1,char *s2,int n) {
  char *end=s2+n;
  for(;;++s1,++s2) {
    while(*s1&&xmlc_white_space(*s1)) ++s1;
    while(s2!=end&&xmlc_white_space(*s2)) ++s2;
    if(s2==end) return *s1;
    if(!*s1) return -*s2;
    if(*s1!=*s2) return *s1-*s2;
  }
}

static int nrmcmpn(char *s1,char *s2,int n) {
  char *end=s2+n;
  for(;;++s1,++s2) {
    if(s2==end) return *s1;
    if(!*s1) return -*s2;
    if(!(*s1==*s2||(xmlc_white_space(*s1)&&xmlc_white_space(*s2))))
      return *s1-*s2;
  }
}

static int qncmpn(char *s1,char *s2,int n2) { /* context is not passed over; compare local parts */
  char *ln1=s1,*ln2=s2;
  int n=n2;
  while(*ln1&&*ln1!=':') ++ln1;
  while(n!=0&&*ln2!=':') {++ln2; --n;}
  if(*ln1) {
    return n?s_tokcmpn(ln1+1,ln2+1,n-1):s_tokcmpn(ln1+1,s2,n2);
  } else {
    return n?s_tokcmpn(s1,ln2+1,n-1):s_tokcmpn(s1,s2,n2);
  }
}

int xsd_equal(char *typ,char *val,char *s,int n) {
  if(!xsd_allows(typ,"",val,strlen(val))) {
    (*error_handler)(XSD_ER_VAL,val);
    return 0;
  }
  if(!xsd_allows(typ,"",s,n)) return 0;
  switch(s_tab(typ,typtab,NTYP)) {
 /*primitive*/
  case TYP_STRING: return s_cmpn(val,s,n)==0;
  case TYP_BOOLEAN: return (s_tokcmpn("true",val,strlen(val))==0||s_tokcmpn("1",val,strlen(val))==0)==(s_tokcmpn("true",s,n)==0||s_tokcmpn("1",s,n)==0);
  case TYP_DECIMAL: return deccmp(val,strlen(val),s,n)==0;
  case TYP_FLOAT: case TYP_DOUBLE: return dblcmpn(val,s,n)==0;
  case TYP_DURATION: return duracmp(val,s,n)==0;
  case TYP_DATE_TIME: return dtcmpn(val,s,n,"ymdtz")==0;
  case TYP_DATE: return dtcmpn(val,s,n,"ymdz")==0;
  case TYP_TIME: return dtcmpn(val,s,n,"tz")==0;
  case TYP_G_YEAR_MONTH: return dtcmpn(val,s,n,"ymz")==0;
  case TYP_G_YEAR: return dtcmpn(val,s,n,"yz")==0;
  case TYP_G_MONTH_DAY: return dtcmpn(val,s,n,"mdz")==0;
  case TYP_G_DAY: return dtcmpn(val,s,n,"dz")==0;
  case TYP_G_MONTH: return dtcmpn(val,s,n,"mz")==0;
  case TYP_HEX_BINARY: return hexcmpn(val,s,n)==0;
  case TYP_BASE64_BINARY: return b64cmpn(val,s,n)==0;
  case TYP_ANY_URI: return s_tokcmpn(val,s,n)==0;
  case TYP_QNAME: case TYP_NOTATION:
    return qncmpn(val,s,n)==0;
 /*derived*/
  case TYP_NORMALIZED_STRING: return nrmcmpn(val,s,n)==0;
  case TYP_TOKEN:
  case TYP_LANGUAGE:
  case TYP_NMTOKEN:
  case TYP_NMTOKENS:
  case TYP_NAME:
  case TYP_NCNAME:
  case TYP_ID:
  case TYP_IDREF:
  case TYP_IDREFS:
  case TYP_ENTITY:
  case TYP_ENTITIES: return s_tokcmpn(val,s,n)==0;
  case TYP_INTEGER:
  case TYP_POSITIVE_INTEGER:
  case TYP_NON_NEGATIVE_INTEGER:
  case TYP_NON_POSITIVE_INTEGER:
  case TYP_NEGATIVE_INTEGER:
  case TYP_BYTE:
  case TYP_UNSIGNED_BYTE:
  case TYP_SHORT:
  case TYP_UNSIGNED_SHORT:
  case TYP_INT:
  case TYP_UNSIGNED_INT:
  case TYP_LONG:
  case TYP_UNSIGNED_LONG: return deccmp(val,strlen(val),s,n)==0;
  case NTYP: (*error_handler)(XSD_ER_TYP,typ); return 0;
  default: assert(0);
  }
  return 0;
}

void xsd_test() {
  rx_init();

  assert(toklenn("",0)==0);
  assert(toklenn("A",1)==1);
  assert(toklenn(" A  ",4)==1);
  assert(toklenn(" A  B  ",7)==3);

  assert(tokcntn("",0)==0);
  assert(tokcntn("A",1)==1);
  assert(tokcntn("AB CD",5)==2);
  assert(tokcntn("   AB  C ",9)==2);

  assert(diglenn(" +14.25",7)==4);
  assert(diglenn("1",1)==1);
  assert(diglenn("0",1)==1);
  assert(diglenn("+00.0",5)==1);

  assert(fdiglenn(".1",2)==1);
  assert(fdiglenn("+0.0140",7)==3);
  assert(fdiglenn("0",1)==0);

  assert(deccmp("0",1,"0.0",3)==0);
  assert(deccmp("1 ",2," 1",2)==0);
  assert(deccmp("0.",2,".0",2)==0);
  assert(deccmp("1",1,"1.0",3)==0);
  assert(deccmp("01.1",4,"1.10",4)==0);
  assert(deccmp("+1",2,"1.0",3)==0);
  assert(deccmp("+0.",3,"-0",2)==0);
  assert(deccmp("0",1,"0.1",3)<0);
  assert(deccmp("1.",2,".0",2)>0);
  assert(deccmp("+1",2,"-1",2)>0);

  assert(hexcmpn("","",0)==0);
  assert(hexcmpn("ABC123","ABC123",6)==0);
  assert(hexcmpn("aBCd","AbCd",4)==0);
  assert(hexcmpn("ABC 123"," ABC123",7)==0);
  assert(hexcmpn("ABC124","ABC123",6)>0);

  assert(rx_match(PAT_BASE64_BINARY,"",0));
  assert(rx_match(PAT_BASE64_BINARY,"YmFz",4));
  assert(rx_match(PAT_BASE64_BINARY,"YA==",4));
  assert(rx_match(PAT_BASE64_BINARY,"Y w = =",7));
  assert(rx_match(PAT_BASE64_BINARY,"YF8=",4));

  assert(!rx_match(PAT_BASE64_BINARY,"YmF@",4));
  assert(!rx_match(PAT_BASE64_BINARY,"YmFgH",5));
  assert(!rx_match(PAT_BASE64_BINARY,"Y===",4));
  assert(!rx_match(PAT_BASE64_BINARY,"YF=O",4));
  assert(!rx_match(PAT_BASE64_BINARY,"YFZ=",4));
 
  assert(b64cmpn("","",0)==0);
  assert(b64cmpn("ABC123","ABC123",6)==0);
  assert(b64cmpn("ABC 123"," ABC123",7)==0);
  assert(b64cmpn("ABC124","ABC123",6)>0);
  assert(b64cmpn("ABC123","abc123",6)<0);

  assert(nrmcmpn("A B","A B",3)==0);
  assert(nrmcmpn("A B","A C",3)<0);
  assert(nrmcmpn("A B","A\nB",3)==0);
  assert(nrmcmpn(" A","A ",2)<0);
}
