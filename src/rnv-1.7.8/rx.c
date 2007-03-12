/* $Id: rx.c,v 1.33 2004/02/25 00:00:32 dvd Exp $ */

#include <string.h> /*strlen,strcpy,strcmp*/
#include <assert.h>
#include "u.h" /*u_get,u_strlen*/
#include "xmlc.h"
#include "m.h"
#include "s.h"
#include "ht.h"
#include "ll.h"
#include "er.h"
#include "rx.h"

#define LEN_P RX_LEN_P
#define PRIME_P RX_PRIME_P
#define LIM_P RX_LIM_P
#define LEN_2 RX_LEN_2
#define PRIME_2 RX_PRIME_2
#define LEN_R RX_LEN_R
#define PRIME_R RX_PRIME_R

#define R_AVG_SIZE 16

/* it is good to have few patterns when deltas are memoized */
#define P_ERROR 0
#define P_NOT_ALLOWED  1
#define P_EMPTY 2
#define P_CHOICE 3
#define P_GROUP 4
#define P_ONE_OR_MORE 5 /*+*/
#define P_EXCEPT 6 /*single-single*/
#define P_RANGE 7 /*lower,upper inclusive*/
#define P_CLASS 8 /*complement is .-*/
#define P_ANY 9
#define P_CHAR 10

#define P_SIZE 3
#define P_AVG_SIZE 2

static int p_size[]={1,1,1,3,3,2,3,3,2,1,2};

#define P_TYP(i) (pattern[i]&0xF)
#define P_IS(i,x)  (x==P_TYP(i))
#define P_CHK(i,x)  assert(P_IS(i,x))

#define P_unop(TYP,p,p1) P_CHK(p,TYP); p1=pattern[p+1]
#define P_binop(TYP,p,p1,p2) P_unop(TYP,p,p1); p2=pattern[p+2]
#define NotAllowed(p) P_CHK(p,P_NotAllowed)
#define Empty(p) P_CHK(p,P_Empty)
#define Any(p) P_CHK(p,P_Any)
#define Choice(p,p1,p2) P_binop(P_CHOICE,p,p1,p2)
#define Group(p,p1,p2) P_binop(P_GROUP,p,p1,p2)
#define OneOrMore(p,p1) P_unop(P_ONE_OR_MORE,p,p1)
#define Except(p,p1,p2) P_binop(P_EXCEPT,p,p1,p2)
#define Range(p,cf,cl) P_binop(P_RANGE,p,cf,cl)
#define Class(p,cn) P_unop(P_CLASS,p,cn)
#define Char(p,c) P_unop(P_CHAR,p,c)

#define P_NUL 0x100

#define setNullable(x) if(x) pattern[i_p]|=P_NUL
#define nullable(p) (pattern[p]&P_NUL)

int rx_compact=0;
/* 'compact' in drv and rx do different things.
 In drv, it limits the size of the table of memoized deltas. In rx, it limits the size
 of the buffer for cached regular expressions; memoized deltas are always limited by LIM_M,
 since the whole repertoire of unicode characters can blow up the buffer.
 */

static char *regex;
static int *pattern;
static int (*r2p)[2];
static struct hashtable ht_r,ht_p,ht_2;
static int i_p,len_p,i_r,len_r,i_2,len_2;
static int empty,notAllowed,any;

static int accept_p(void) {
  int j;
  if((j=ht_get(&ht_p,i_p))==-1) {
    ht_put(&ht_p,j=i_p);
    i_p+=p_size[P_TYP(i_p)];
    if(i_p+P_SIZE>len_p) pattern=(int*)m_stretch(pattern,len_p=2*(i_p+P_SIZE),i_p,sizeof(int));
  }
  return j;
}

#define P_NEW(x) (pattern[i_p]=x)

#define P_newunop(TYP,p1) P_NEW(TYP); pattern[i_p+1]=p1
#define P_newbinop(TYP,p1,p2) P_newunop(TYP,p1); pattern[i_p+2]=p2
static int newNotAllowed(void) {P_NEW(P_NOT_ALLOWED); return accept_p();}
static int newEmpty(void) {P_NEW(P_EMPTY); setNullable(1); return accept_p();}
static int newAny(void) {P_NEW(P_ANY); return accept_p();}
static int newChoice(int p1,int p2) {P_newbinop(P_CHOICE,p1,p2); setNullable(nullable(p1)||nullable(p2)); return accept_p();}
static int newGroup(int p1,int p2) {P_newbinop(P_GROUP,p1,p2); setNullable(nullable(p1)&&nullable(p2)); return accept_p();}
static int newOneOrMore(int p1) {P_newunop(P_ONE_OR_MORE,p1); setNullable(nullable(p1)); return accept_p();}
static int newExcept(int p1,int p2) {P_newbinop(P_EXCEPT,p1,p2); return accept_p();}
static int newRange(int cf,int cl) {P_newbinop(P_RANGE,cf,cl); return accept_p();}
static int newClass(int cn) {P_newunop(P_CLASS,cn); return accept_p();}
static int newChar(int c) {P_newunop(P_CHAR,c); return accept_p();}

static int one_or_more(int p) {
  if(P_IS(p,P_EMPTY)) return p;
  if(P_IS(p,P_NOT_ALLOWED)) return p;
  return newOneOrMore(p);
}

static int group(int p1,int p2) {
  if(P_IS(p1,P_NOT_ALLOWED)) return p1;
  if(P_IS(p2,P_NOT_ALLOWED)) return p2;
  if(P_IS(p1,P_EMPTY)) return p2;
  if(P_IS(p2,P_EMPTY)) return p1;
  return newGroup(p1,p2);
}

static int samechoice(int p1,int p2) {
  if(P_IS(p1,P_CHOICE)) {
    int p11,p12; Choice(p1,p11,p12);
    return p12==p2||samechoice(p11,p2);
  } else return p1==p2;
}

static int choice(int p1,int p2) {
  if(P_IS(p1,P_NOT_ALLOWED)) return p2;
  if(P_IS(p2,P_NOT_ALLOWED)) return p1;
  if(P_IS(p2,P_CHOICE)) {
    int p21,p22; Choice(p2,p21,p22);
    p1=choice(p1,p21); return choice(p1,p22);
  }
  if(samechoice(p1,p2)) return p1;
  if(nullable(p1) && (P_IS(p2,P_EMPTY))) return p1;
  if(nullable(p2) && (P_IS(p1,P_EMPTY))) return p2;
  return newChoice(p1,p2);
}

static int cls(int cn) {
  if(cn<0) return newExcept(any,newClass(-cn));
  if(cn==0) return notAllowed;
  return newClass(cn);
}

static int equal_r(int r1,int r2) {return strcmp(regex+r1,regex+r2)==0;}
static int hash_r(int r) {return s_hval(regex+r);}

static int equal_p(int p1,int p2) {
  int *pp1=pattern+p1,*pp2=pattern+p2;
  if(P_TYP(p1)!=P_TYP(p2)) return 0;
  switch(p_size[P_TYP(p1)]) {
  case 3: if(pp1[2]!=pp2[2]) return 0;
  case 2: if(pp1[1]!=pp2[1]) return 0;
  case 1: return 1;
  default: assert(0);
  }
  return 0;
}
static int hash_p(int p) {
  int *pp=pattern+p; int h=0;
  switch(p_size[P_TYP(p)]) {
  case 1: h=pp[0]&0xF; break;
  case 2: h=(pp[0]&0xF)|(pp[1]<<4); break;
  case 3: h=(pp[0]&0xF)|((pp[1]^pp[2])<<4); break;
  default: assert(0);
  }
  return h*PRIME_P;
}

static int equal_2(int x1,int x2) {return r2p[x1][0]==r2p[x2][0];}
static int hash_2(int x) {return r2p[x][0]*PRIME_2;}

static int add_r(char *rx) {
  int len=strlen(rx)+1;
  if(i_r+len>len_r) regex=(char*)m_stretch(regex,len_r=2*(i_r+len),i_r,sizeof(char));
  strcpy(regex+i_r,rx);
  return len;
}

#define ERRPOS

#define err(msg) (*er_vprintf)(msg" in \"%s\" at offset %i\n",ap)
void rx_default_verror_handler(int erno,va_list ap) {
  (*er_printf)("regular expressions: ");
  switch(erno) {
  case RX_ER_BADCH: err("bad character"); break;
  case RX_ER_UNFIN: err("unfinished expression"); break;
  case RX_ER_NOLSQ: err("'[' expected"); break;
  case RX_ER_NORSQ: err("']' expected"); break;
  case RX_ER_NOLCU: err("'{' expected"); break;
  case RX_ER_NORCU: err("'}' expected"); break;
  case RX_ER_NOLPA: err("'(' expected"); break;
  case RX_ER_NORPA: err("')' expected"); break;
  case RX_ER_BADCL: err("unknown class"); break;
  case RX_ER_NODGT: err("digit expected"); break;
  case RX_ER_DNUOB: err("reversed bounds"); break;
  case RX_ER_NOTRC: err("range or class expected"); break;
  default: assert(0);
  }
}

void (*rx_verror_handler)(int erno,va_list ap)=&rx_default_verror_handler;

static void error_handler(int erno,...) {
  va_list ap; va_start(ap,erno); (*rx_verror_handler)(erno,ap); va_end(ap);
}

#define LEN_M RX_LEN_M
#define PRIME_M RX_PRIME_M
#define LIM_M RX_LIM_M

#define M_SIZE 3

#define M_SET(p) memo[i_m][M_SIZE-1]=p
#define M_RET(m) memo[m][M_SIZE-1]

static int (*memo)[M_SIZE];
static int i_m,len_m;
static struct hashtable ht_m;

static int new_memo(int p,int c) {
  int *me=memo[i_m];
  ht_deli(&ht_m,i_m);
  me[0]=p; me[1]=c;
  return ht_get(&ht_m,i_m);
}

static int equal_m(int m1,int m2) {
  int *me1=memo[m1],*me2=memo[m2];
  return (me1[0]==me2[0])&&(me1[1]==me2[1]);
}
static int hash_m(int m) {
  int *me=memo[m];
  return (me[0]^me[1])*PRIME_M;
}

static void accept_m(void) {
  if(ht_get(&ht_m,i_m)!=-1) ht_del(&ht_m,i_m);
  ht_put(&ht_m,i_m++);
  if(i_m>=LIM_M) i_m=0;
  if(i_m==len_m) memo=(int(*)[M_SIZE])m_stretch(memo,len_m=i_m*2,i_m,sizeof(int[M_SIZE]));
}

static void windup(void);
static int initialized=0;
void rx_init(void) {
  if(!initialized) { initialized=1;
    pattern=(int *)m_alloc(len_p=P_AVG_SIZE*LEN_P,sizeof(int));
    r2p=(int (*)[2])m_alloc(len_2=LEN_2,sizeof(int[2]));
    regex=(char*)m_alloc(len_r=R_AVG_SIZE*LEN_R,sizeof(char));
    memo=(int (*)[M_SIZE])m_alloc(len_m=LEN_M,sizeof(int[M_SIZE]));

    ht_init(&ht_p,LEN_P,&hash_p,&equal_p);
    ht_init(&ht_2,LEN_2,&hash_2,&equal_2);
    ht_init(&ht_r,LEN_R,&hash_r,&equal_r);
    ht_init(&ht_m,LEN_M,&hash_m,&equal_m);

    windup();
  }
}

void rx_clear(void) {
  ht_clear(&ht_p); ht_clear(&ht_2); ht_clear(&ht_r); ht_clear(&ht_m);
  windup();
}

static void windup(void) {
  i_p=i_r=i_2=i_m=0;
  pattern[0]=P_ERROR;  accept_p();
  empty=newEmpty(); notAllowed=newNotAllowed(); any=newAny();
}

#define SYM_END 0
#define SYM_CLS 1
#define SYM_ESC 2
#define SYM_CHR 3

static int r0,ri,sym,val,errors;

static void error(int erno) {
  if(!errors) error_handler(erno,regex+r0,u_strlen(regex+r0)-u_strlen(regex+ri));
  ++errors;
}

#include "rx_cls_u.c"

static int chclass(void) {
  int u,cl,rj;
  ri+=u_get(&u,regex+ri);
  if(u=='\0') {--ri; error(RX_ER_NOLCU); return 0;}
  if(u!='{') {error(RX_ER_NOLCU); return 0;}
  rj=ri;
  for(;;) {
    if(regex[rj]=='\0') {ri=rj; error(RX_ER_NORCU); return 0;}
    if(regex[rj]=='}') {
      if((cl=s_ntab(regex+ri,rj-ri,clstab,NUM_CLS_U))==NUM_CLS_U) {error(RX_ER_BADCL); cl=0;}
      ri=rj+1;
      return cl;
    }
    ++rj;
  }
}

#define CLS_NL (NUM_CLS_U+1)
#define CLS_S (NUM_CLS_U+2)
#define CLS_I (NUM_CLS_U+3)
#define CLS_C (NUM_CLS_U+4)
#define CLS_W (NUM_CLS_U+5)
#define NUM_CLS (NUM_CLS_U+6)

static void getsym(void) {
  int u;
  if(regex[ri]=='\0') sym=SYM_END; else {
    ri+=u_get(&u,regex+ri);
    if(u=='\\') {
      ri+=u_get(&u,regex+ri);
      switch(u) {
      case '\0': --ri; error(RX_ER_UNFIN); sym=SYM_END; break;
      case 'p': sym=SYM_CLS; val=chclass(); break;
      case 'P': sym=SYM_CLS; val=-chclass(); break;
      case 's': sym=SYM_CLS; val=CLS_S; break;
      case 'S': sym=SYM_CLS; val=-CLS_S; break;
      case 'i': sym=SYM_CLS; val=CLS_I; break;
      case 'I': sym=SYM_CLS; val=-CLS_I; break;
      case 'c': sym=SYM_CLS; val=CLS_C; break;
      case 'C': sym=SYM_CLS; val=-CLS_C; break;
      case 'd': sym=SYM_CLS; val=CLS_U_Nd; break;
      case 'D': sym=SYM_CLS; val=-CLS_U_Nd; break;
      case 'w': sym=SYM_CLS; val=CLS_W; break;
      case 'W': sym=SYM_CLS; val=-CLS_W; break;
      case 'n': sym=SYM_ESC; val=0xA; break;
      case 'r': sym=SYM_ESC; val=0xD; break;
      case 't': sym=SYM_ESC; val=0x9; break;
      case '\\': case '|': case '.': case '-': case '^': case '?': case '*': case '+':
      case '{': case '}': case '[': case ']': case '(': case ')':
	sym=SYM_ESC; val=u; break;
      default: error(RX_ER_BADCH); sym=SYM_ESC; val=u; break;
      }
    } else {
      switch(u) {
      case '.': sym=SYM_CLS; val=-CLS_NL; break;
      default: sym=SYM_CHR; val=u; break;
      }
    }
  }
}

static void chk_get(int v,int erno) {if(sym!=SYM_CHR||val!=v) error(erno); getsym();}


#define chkrch(val) if((val)=='['||(val)==']'||(val)=='-') error(RX_ER_NOTRC)

static int chgroup(void) {
  int p=notAllowed,c;
  for(;;) {
    switch(sym) {
    case SYM_CHR: chkrch(val);
    case SYM_ESC: c=val; getsym();
      if(sym==SYM_CHR&&val=='-') {
	if(regex[ri]=='[') {
	  p=choice(p,newChar(c));
	  goto END_OF_GROUP;
	} else {
	  getsym();
	  switch(sym) {
	  case SYM_CHR: chkrch(val);
	  case SYM_ESC: p=choice(p,newRange(c,val)); getsym(); break;
	  default: error(RX_ER_BADCH); getsym(); break;
	  }
	}
      } else {
	p=choice(p,newChar(c));
      }
      break;
    case SYM_CLS: p=choice(p,cls(val)); getsym(); break;
    case SYM_END: error(RX_ER_NORSQ); goto END_OF_GROUP;
    default: assert(0);
    }
    if(sym==SYM_CHR&&(val==']'||val=='-')) goto END_OF_GROUP;
  }
  END_OF_GROUP:;
  return p;
}

static int chexpr(void) {
  int p;
  if(sym==SYM_CHR&&val=='^') { getsym();
    p=newExcept(any,chgroup());
  } else {
    p=chgroup();
  }
  if(sym==SYM_CHR&&val=='-') { getsym();
    chk_get('[',RX_ER_NOLSQ); p=newExcept(p,chexpr()); chk_get(']',RX_ER_NORSQ);
  }
  return p;
}

static int expression(void);
static int atom(void) {
  int p=0;
  switch(sym) {
  case SYM_CHR:
    switch(val) {
    case '[': getsym(); p=chexpr(); chk_get(']',RX_ER_NORSQ); break;
    case '(': getsym(); p=expression(); chk_get(')',RX_ER_NORPA); break;
    case '{': case '?': case '*': case '+': case '|':
    case ')': case ']': case '}': error(RX_ER_BADCH); getsym(); break;
    default: p=newChar(val); getsym(); break;
    }
    break;
  case SYM_ESC: p=newChar(val); getsym(); break;
  case SYM_CLS: p=cls(val); getsym(); break;
  default: error(RX_ER_BADCH); getsym(); break;
  }
  return p;
}

static int number(void) {
  int n=0,m;
  for(;;) {
    if(sym!=SYM_CHR) goto END_OF_DIGITS;
    switch(val) {
    case '0': m=0; break;
    case '1': m=1; break;
    case '2': m=2; break;
    case '3': m=3; break;
    case '4': m=4; break;
    case '5': m=5; break;
    case '6': m=6; break;
    case '7': m=7; break;
    case '8': m=8; break;
    case '9': m=9; break;
    default: goto END_OF_DIGITS;
    }
    n=n*10+m;
    getsym();
  }
  END_OF_DIGITS:;
  return n;
}

static int quantifier(int p0) {
  int p=empty,n,n0;
  n=n0=number();
  while(n--) p=group(p,p0);
  if(sym==SYM_CHR) {
    if(val==',') {
      getsym();
      if(sym==SYM_CHR && val=='}') {
	p=group(p,choice(empty,one_or_more(p0)));
      } else {
	n=number()-n0; if(n<0) {error(RX_ER_DNUOB); n=0;}
	while(n--) p=group(p,choice(empty,p0));
      }
    }
  } else error(RX_ER_NODGT);
  return p;
}

static int piece(void) {
  int p;
  p=atom();
  if(sym==SYM_CHR) {
    switch(val) {
    case '{': getsym(); p=quantifier(p); chk_get('}',RX_ER_NOLCU); break;
    case '?': getsym(); p=choice(empty,p); break;
    case '*': getsym(); p=choice(empty,one_or_more(p)); break;
    case '+': getsym(); p=one_or_more(p); break;
    default: break;
    }
  }
  return p;
}

static int branch(void) {
  int p;
  p=empty;
  while(!(sym==SYM_END||(sym==SYM_CHR&&(val=='|'||val==')')))) p=group(p,piece());
  return p;
}

static int expression(void) {
  int p;
  p=branch();
  while(sym==SYM_CHR&&val=='|') {
    getsym();
    p=choice(p,branch());
  }
  return p;
}

static void bind(int r) {
  r0=ri=r; sym=-1; errors=0;
  getsym();
}

static int compile(char *rx) {
  int r=0,p=0,d_r;
  d_r=add_r(rx);
  if((r=ht_get(&ht_r,i_r))==-1) {
    if(rx_compact&&i_p>=P_AVG_SIZE*LIM_P) {rx_clear(); d_r=add_r(rx);}
    ht_put(&ht_r,r=i_r);
    i_r+=d_r;
    bind(r); p=expression(); if(sym!=SYM_END) error(RX_ER_BADCH);
    r2p[i_2][0]=r; r2p[i_2][1]=p;
    ht_put(&ht_2,i_2++);
    if(i_2==len_2) r2p=(int(*)[2])m_stretch(r2p,len_2=2*i_2,i_2,sizeof(int[2]));
  } else {
    r2p[i_2][0]=r;
    p=r2p[ht_get(&ht_2,i_2)][1];
  }
  return p;
}

#include "rx_cls_ranges.c"

static int in_class(int c,int cn) {
  switch(cn) {
  case 0: return 0;
  case CLS_U_C: return in_class(c,CLS_U_Cc)||in_class(c,CLS_U_Cf)||in_class(c,CLS_U_Co);
  case CLS_U_Cc: return u_in_ranges(c,CcRanges,sizeof(CcRanges)/sizeof(int[2]));
  case CLS_U_Cf: return u_in_ranges(c,CfRanges,sizeof(CfRanges)/sizeof(int[2]));
  case CLS_U_Co: return u_in_ranges(c,CoRanges,sizeof(CoRanges)/sizeof(int[2]));
  case CLS_U_IsAlphabeticPresentationForms: return u_in_ranges(c,IsAlphabeticPresentationFormsRanges,sizeof(IsAlphabeticPresentationFormsRanges)/sizeof(int[2]));
  case CLS_U_IsArabic: return u_in_ranges(c,IsArabicRanges,sizeof(IsArabicRanges)/sizeof(int[2]));
  case CLS_U_IsArabicPresentationForms_A: return u_in_ranges(c,IsArabicPresentationForms_ARanges,sizeof(IsArabicPresentationForms_ARanges)/sizeof(int[2]));
  case CLS_U_IsArabicPresentationForms_B: return u_in_ranges(c,IsArabicPresentationForms_BRanges,sizeof(IsArabicPresentationForms_BRanges)/sizeof(int[2]));
  case CLS_U_IsArmenian: return u_in_ranges(c,IsArmenianRanges,sizeof(IsArmenianRanges)/sizeof(int[2]));
  case CLS_U_IsArrows: return u_in_ranges(c,IsArrowsRanges,sizeof(IsArrowsRanges)/sizeof(int[2]));
  case CLS_U_IsBasicLatin: return u_in_ranges(c,IsBasicLatinRanges,sizeof(IsBasicLatinRanges)/sizeof(int[2]));
  case CLS_U_IsBengali: return u_in_ranges(c,IsBengaliRanges,sizeof(IsBengaliRanges)/sizeof(int[2]));
  case CLS_U_IsBlockElements: return u_in_ranges(c,IsBlockElementsRanges,sizeof(IsBlockElementsRanges)/sizeof(int[2]));
  case CLS_U_IsBopomofo: return u_in_ranges(c,IsBopomofoRanges,sizeof(IsBopomofoRanges)/sizeof(int[2]));
  case CLS_U_IsBopomofoExtended: return u_in_ranges(c,IsBopomofoExtendedRanges,sizeof(IsBopomofoExtendedRanges)/sizeof(int[2]));
  case CLS_U_IsBoxDrawing: return u_in_ranges(c,IsBoxDrawingRanges,sizeof(IsBoxDrawingRanges)/sizeof(int[2]));
  case CLS_U_IsBraillePatterns: return u_in_ranges(c,IsBraillePatternsRanges,sizeof(IsBraillePatternsRanges)/sizeof(int[2]));
  case CLS_U_IsByzantineMusicalSymbols: return u_in_ranges(c,IsByzantineMusicalSymbolsRanges,sizeof(IsByzantineMusicalSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsCJKCompatibility: return u_in_ranges(c,IsCJKCompatibilityRanges,sizeof(IsCJKCompatibilityRanges)/sizeof(int[2]));
  case CLS_U_IsCJKCompatibilityForms: return u_in_ranges(c,IsCJKCompatibilityFormsRanges,sizeof(IsCJKCompatibilityFormsRanges)/sizeof(int[2]));
  case CLS_U_IsCJKCompatibilityIdeographs: return u_in_ranges(c,IsCJKCompatibilityIdeographsRanges,sizeof(IsCJKCompatibilityIdeographsRanges)/sizeof(int[2]));
  case CLS_U_IsCJKCompatibilityIdeographsSupplement: return u_in_ranges(c,IsCJKCompatibilityIdeographsSupplementRanges,sizeof(IsCJKCompatibilityIdeographsSupplementRanges)/sizeof(int[2]));
  case CLS_U_IsCJKRadicalsSupplement: return u_in_ranges(c,IsCJKRadicalsSupplementRanges,sizeof(IsCJKRadicalsSupplementRanges)/sizeof(int[2]));
  case CLS_U_IsCJKSymbolsandPunctuation: return u_in_ranges(c,IsCJKSymbolsandPunctuationRanges,sizeof(IsCJKSymbolsandPunctuationRanges)/sizeof(int[2]));
  case CLS_U_IsCJKUnifiedIdeographs: return u_in_ranges(c,IsCJKUnifiedIdeographsRanges,sizeof(IsCJKUnifiedIdeographsRanges)/sizeof(int[2]));
  case CLS_U_IsCJKUnifiedIdeographsExtensionA: return u_in_ranges(c,IsCJKUnifiedIdeographsExtensionARanges,sizeof(IsCJKUnifiedIdeographsExtensionARanges)/sizeof(int[2]));
  case CLS_U_IsCJKUnifiedIdeographsExtensionB: return u_in_ranges(c,IsCJKUnifiedIdeographsExtensionBRanges,sizeof(IsCJKUnifiedIdeographsExtensionBRanges)/sizeof(int[2]));
  case CLS_U_IsCherokee: return u_in_ranges(c,IsCherokeeRanges,sizeof(IsCherokeeRanges)/sizeof(int[2]));
  case CLS_U_IsCombiningDiacriticalMarks: return u_in_ranges(c,IsCombiningDiacriticalMarksRanges,sizeof(IsCombiningDiacriticalMarksRanges)/sizeof(int[2]));
  case CLS_U_IsCombiningHalfMarks: return u_in_ranges(c,IsCombiningHalfMarksRanges,sizeof(IsCombiningHalfMarksRanges)/sizeof(int[2]));
  case CLS_U_IsCombiningMarksforSymbols: return u_in_ranges(c,IsCombiningMarksforSymbolsRanges,sizeof(IsCombiningMarksforSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsControlPictures: return u_in_ranges(c,IsControlPicturesRanges,sizeof(IsControlPicturesRanges)/sizeof(int[2]));
  case CLS_U_IsCurrencySymbols: return u_in_ranges(c,IsCurrencySymbolsRanges,sizeof(IsCurrencySymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsCyrillic: return u_in_ranges(c,IsCyrillicRanges,sizeof(IsCyrillicRanges)/sizeof(int[2]));
  case CLS_U_IsDeseret: return u_in_ranges(c,IsDeseretRanges,sizeof(IsDeseretRanges)/sizeof(int[2]));
  case CLS_U_IsDevanagari: return u_in_ranges(c,IsDevanagariRanges,sizeof(IsDevanagariRanges)/sizeof(int[2]));
  case CLS_U_IsDingbats: return u_in_ranges(c,IsDingbatsRanges,sizeof(IsDingbatsRanges)/sizeof(int[2]));
  case CLS_U_IsEnclosedAlphanumerics: return u_in_ranges(c,IsEnclosedAlphanumericsRanges,sizeof(IsEnclosedAlphanumericsRanges)/sizeof(int[2]));
  case CLS_U_IsEnclosedCJKLettersandMonths: return u_in_ranges(c,IsEnclosedCJKLettersandMonthsRanges,sizeof(IsEnclosedCJKLettersandMonthsRanges)/sizeof(int[2]));
  case CLS_U_IsEthiopic: return u_in_ranges(c,IsEthiopicRanges,sizeof(IsEthiopicRanges)/sizeof(int[2]));
  case CLS_U_IsGeneralPunctuation: return u_in_ranges(c,IsGeneralPunctuationRanges,sizeof(IsGeneralPunctuationRanges)/sizeof(int[2]));
  case CLS_U_IsGeometricShapes: return u_in_ranges(c,IsGeometricShapesRanges,sizeof(IsGeometricShapesRanges)/sizeof(int[2]));
  case CLS_U_IsGeorgian: return u_in_ranges(c,IsGeorgianRanges,sizeof(IsGeorgianRanges)/sizeof(int[2]));
  case CLS_U_IsGothic: return u_in_ranges(c,IsGothicRanges,sizeof(IsGothicRanges)/sizeof(int[2]));
  case CLS_U_IsGreek: return u_in_ranges(c,IsGreekRanges,sizeof(IsGreekRanges)/sizeof(int[2]));
  case CLS_U_IsGreekExtended: return u_in_ranges(c,IsGreekExtendedRanges,sizeof(IsGreekExtendedRanges)/sizeof(int[2]));
  case CLS_U_IsGujarati: return u_in_ranges(c,IsGujaratiRanges,sizeof(IsGujaratiRanges)/sizeof(int[2]));
  case CLS_U_IsGurmukhi: return u_in_ranges(c,IsGurmukhiRanges,sizeof(IsGurmukhiRanges)/sizeof(int[2]));
  case CLS_U_IsHalfwidthandFullwidthForms: return u_in_ranges(c,IsHalfwidthandFullwidthFormsRanges,sizeof(IsHalfwidthandFullwidthFormsRanges)/sizeof(int[2]));
  case CLS_U_IsHangulCompatibilityJamo: return u_in_ranges(c,IsHangulCompatibilityJamoRanges,sizeof(IsHangulCompatibilityJamoRanges)/sizeof(int[2]));
  case CLS_U_IsHangulJamo: return u_in_ranges(c,IsHangulJamoRanges,sizeof(IsHangulJamoRanges)/sizeof(int[2]));
  case CLS_U_IsHangulSyllables: return u_in_ranges(c,IsHangulSyllablesRanges,sizeof(IsHangulSyllablesRanges)/sizeof(int[2]));
  case CLS_U_IsHebrew: return u_in_ranges(c,IsHebrewRanges,sizeof(IsHebrewRanges)/sizeof(int[2]));
  case CLS_U_IsHiragana: return u_in_ranges(c,IsHiraganaRanges,sizeof(IsHiraganaRanges)/sizeof(int[2]));
  case CLS_U_IsIPAExtensions: return u_in_ranges(c,IsIPAExtensionsRanges,sizeof(IsIPAExtensionsRanges)/sizeof(int[2]));
  case CLS_U_IsIdeographicDescriptionCharacters: return u_in_ranges(c,IsIdeographicDescriptionCharactersRanges,sizeof(IsIdeographicDescriptionCharactersRanges)/sizeof(int[2]));
  case CLS_U_IsKanbun: return u_in_ranges(c,IsKanbunRanges,sizeof(IsKanbunRanges)/sizeof(int[2]));
  case CLS_U_IsKangxiRadicals: return u_in_ranges(c,IsKangxiRadicalsRanges,sizeof(IsKangxiRadicalsRanges)/sizeof(int[2]));
  case CLS_U_IsKannada: return u_in_ranges(c,IsKannadaRanges,sizeof(IsKannadaRanges)/sizeof(int[2]));
  case CLS_U_IsKatakana: return u_in_ranges(c,IsKatakanaRanges,sizeof(IsKatakanaRanges)/sizeof(int[2]));
  case CLS_U_IsKhmer: return u_in_ranges(c,IsKhmerRanges,sizeof(IsKhmerRanges)/sizeof(int[2]));
  case CLS_U_IsLao: return u_in_ranges(c,IsLaoRanges,sizeof(IsLaoRanges)/sizeof(int[2]));
  case CLS_U_IsLatin_1Supplement: return u_in_ranges(c,IsLatin_1SupplementRanges,sizeof(IsLatin_1SupplementRanges)/sizeof(int[2]));
  case CLS_U_IsLatinExtended_A: return u_in_ranges(c,IsLatinExtended_ARanges,sizeof(IsLatinExtended_ARanges)/sizeof(int[2]));
  case CLS_U_IsLatinExtended_B: return u_in_ranges(c,IsLatinExtended_BRanges,sizeof(IsLatinExtended_BRanges)/sizeof(int[2]));
  case CLS_U_IsLatinExtendedAdditional: return u_in_ranges(c,IsLatinExtendedAdditionalRanges,sizeof(IsLatinExtendedAdditionalRanges)/sizeof(int[2]));
  case CLS_U_IsLetterlikeSymbols: return u_in_ranges(c,IsLetterlikeSymbolsRanges,sizeof(IsLetterlikeSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsMalayalam: return u_in_ranges(c,IsMalayalamRanges,sizeof(IsMalayalamRanges)/sizeof(int[2]));
  case CLS_U_IsMathematicalAlphanumericSymbols: return u_in_ranges(c,IsMathematicalAlphanumericSymbolsRanges,sizeof(IsMathematicalAlphanumericSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsMathematicalOperators: return u_in_ranges(c,IsMathematicalOperatorsRanges,sizeof(IsMathematicalOperatorsRanges)/sizeof(int[2]));
  case CLS_U_IsMiscellaneousSymbols: return u_in_ranges(c,IsMiscellaneousSymbolsRanges,sizeof(IsMiscellaneousSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsMiscellaneousTechnical: return u_in_ranges(c,IsMiscellaneousTechnicalRanges,sizeof(IsMiscellaneousTechnicalRanges)/sizeof(int[2]));
  case CLS_U_IsMongolian: return u_in_ranges(c,IsMongolianRanges,sizeof(IsMongolianRanges)/sizeof(int[2]));
  case CLS_U_IsMusicalSymbols: return u_in_ranges(c,IsMusicalSymbolsRanges,sizeof(IsMusicalSymbolsRanges)/sizeof(int[2]));
  case CLS_U_IsMyanmar: return u_in_ranges(c,IsMyanmarRanges,sizeof(IsMyanmarRanges)/sizeof(int[2]));
  case CLS_U_IsNumberForms: return u_in_ranges(c,IsNumberFormsRanges,sizeof(IsNumberFormsRanges)/sizeof(int[2]));
  case CLS_U_IsOgham: return u_in_ranges(c,IsOghamRanges,sizeof(IsOghamRanges)/sizeof(int[2]));
  case CLS_U_IsOldItalic: return u_in_ranges(c,IsOldItalicRanges,sizeof(IsOldItalicRanges)/sizeof(int[2]));
  case CLS_U_IsOpticalCharacterRecognition: return u_in_ranges(c,IsOpticalCharacterRecognitionRanges,sizeof(IsOpticalCharacterRecognitionRanges)/sizeof(int[2]));
  case CLS_U_IsOriya: return u_in_ranges(c,IsOriyaRanges,sizeof(IsOriyaRanges)/sizeof(int[2]));
  case CLS_U_IsPrivateUse: return u_in_ranges(c,IsPrivateUseRanges,sizeof(IsPrivateUseRanges)/sizeof(int[2]));
  case CLS_U_IsRunic: return u_in_ranges(c,IsRunicRanges,sizeof(IsRunicRanges)/sizeof(int[2]));
  case CLS_U_IsSinhala: return u_in_ranges(c,IsSinhalaRanges,sizeof(IsSinhalaRanges)/sizeof(int[2]));
  case CLS_U_IsSmallFormVariants: return u_in_ranges(c,IsSmallFormVariantsRanges,sizeof(IsSmallFormVariantsRanges)/sizeof(int[2]));
  case CLS_U_IsSpacingModifierLetters: return u_in_ranges(c,IsSpacingModifierLettersRanges,sizeof(IsSpacingModifierLettersRanges)/sizeof(int[2]));
  case CLS_U_IsSpecials: return u_in_ranges(c,IsSpecialsRanges,sizeof(IsSpecialsRanges)/sizeof(int[2]));
  case CLS_U_IsSuperscriptsandSubscripts: return u_in_ranges(c,IsSuperscriptsandSubscriptsRanges,sizeof(IsSuperscriptsandSubscriptsRanges)/sizeof(int[2]));
  case CLS_U_IsSyriac: return u_in_ranges(c,IsSyriacRanges,sizeof(IsSyriacRanges)/sizeof(int[2]));
  case CLS_U_IsTags: return u_in_ranges(c,IsTagsRanges,sizeof(IsTagsRanges)/sizeof(int[2]));
  case CLS_U_IsTamil: return u_in_ranges(c,IsTamilRanges,sizeof(IsTamilRanges)/sizeof(int[2]));
  case CLS_U_IsTelugu: return u_in_ranges(c,IsTeluguRanges,sizeof(IsTeluguRanges)/sizeof(int[2]));
  case CLS_U_IsThaana: return u_in_ranges(c,IsThaanaRanges,sizeof(IsThaanaRanges)/sizeof(int[2]));
  case CLS_U_IsThai: return u_in_ranges(c,IsThaiRanges,sizeof(IsThaiRanges)/sizeof(int[2]));
  case CLS_U_IsTibetan: return u_in_ranges(c,IsTibetanRanges,sizeof(IsTibetanRanges)/sizeof(int[2]));
  case CLS_U_IsUnifiedCanadianAboriginalSyllabics: return u_in_ranges(c,IsUnifiedCanadianAboriginalSyllabicsRanges,sizeof(IsUnifiedCanadianAboriginalSyllabicsRanges)/sizeof(int[2]));
  case CLS_U_IsYiRadicals: return u_in_ranges(c,IsYiRadicalsRanges,sizeof(IsYiRadicalsRanges)/sizeof(int[2]));
  case CLS_U_IsYiSyllables: return u_in_ranges(c,IsYiSyllablesRanges,sizeof(IsYiSyllablesRanges)/sizeof(int[2]));
  case CLS_U_L: return in_class(c,CLS_U_Ll)||in_class(c,CLS_U_Lm)||in_class(c,CLS_U_Lo)||in_class(c,CLS_U_Lt)||in_class(c,CLS_U_Lu);
  case CLS_U_Ll: return u_in_ranges(c,LlRanges,sizeof(LlRanges)/sizeof(int[2]));
  case CLS_U_Lm: return u_in_ranges(c,LmRanges,sizeof(LmRanges)/sizeof(int[2]));
  case CLS_U_Lo: return u_in_ranges(c,LoRanges,sizeof(LoRanges)/sizeof(int[2]));
  case CLS_U_Lt: return u_in_ranges(c,LtRanges,sizeof(LtRanges)/sizeof(int[2]));
  case CLS_U_Lu: return u_in_ranges(c,LuRanges,sizeof(LuRanges)/sizeof(int[2]));
  case CLS_U_M: return in_class(c,CLS_U_Mc)||in_class(c,CLS_U_Me)||in_class(c,CLS_U_Mn);
  case CLS_U_Mc: return u_in_ranges(c,McRanges,sizeof(McRanges)/sizeof(int[2]));
  case CLS_U_Me: return u_in_ranges(c,MeRanges,sizeof(MeRanges)/sizeof(int[2]));
  case CLS_U_Mn: return u_in_ranges(c,MnRanges,sizeof(MnRanges)/sizeof(int[2]));
  case CLS_U_N: return in_class(c,CLS_U_Nd)||in_class(c,CLS_U_Nl)||in_class(c,CLS_U_No);
  case CLS_U_Nd: return u_in_ranges(c,NdRanges,sizeof(NdRanges)/sizeof(int[2]));
  case CLS_U_Nl: return u_in_ranges(c,NlRanges,sizeof(NlRanges)/sizeof(int[2]));
  case CLS_U_No: return u_in_ranges(c,NoRanges,sizeof(NoRanges)/sizeof(int[2]));
  case CLS_U_P: return in_class(c,CLS_U_Pc)||in_class(c,CLS_U_Pd)||in_class(c,CLS_U_Pe)||in_class(c,CLS_U_Pf)||in_class(c,CLS_U_Pi)||in_class(c,CLS_U_Po)||in_class(c,CLS_U_Ps);
  case CLS_U_Pc: return u_in_ranges(c,PcRanges,sizeof(PcRanges)/sizeof(int[2]));
  case CLS_U_Pd: return u_in_ranges(c,PdRanges,sizeof(PdRanges)/sizeof(int[2]));
  case CLS_U_Pe: return u_in_ranges(c,PeRanges,sizeof(PeRanges)/sizeof(int[2]));
  case CLS_U_Pf: return u_in_ranges(c,PfRanges,sizeof(PfRanges)/sizeof(int[2]));
  case CLS_U_Pi: return u_in_ranges(c,PiRanges,sizeof(PiRanges)/sizeof(int[2]));
  case CLS_U_Po: return u_in_ranges(c,PoRanges,sizeof(PoRanges)/sizeof(int[2]));
  case CLS_U_Ps: return u_in_ranges(c,PsRanges,sizeof(PsRanges)/sizeof(int[2]));
  case CLS_U_S: return in_class(c,CLS_U_Sc)||in_class(c,CLS_U_Sk)||in_class(c,CLS_U_Sm)||in_class(c,CLS_U_So);
  case CLS_U_Sc: return u_in_ranges(c,ScRanges,sizeof(ScRanges)/sizeof(int[2]));
  case CLS_U_Sk: return u_in_ranges(c,SkRanges,sizeof(SkRanges)/sizeof(int[2]));
  case CLS_U_Sm: return u_in_ranges(c,SmRanges,sizeof(SmRanges)/sizeof(int[2]));
  case CLS_U_So: return u_in_ranges(c,SoRanges,sizeof(SoRanges)/sizeof(int[2]));
  case CLS_U_Z: return in_class(c,CLS_U_Zl)||in_class(c,CLS_U_Zp)||in_class(c,CLS_U_Zs);
  case CLS_U_Zl: return u_in_ranges(c,ZlRanges,sizeof(ZlRanges)/sizeof(int[2]));
  case CLS_U_Zp: return u_in_ranges(c,ZpRanges,sizeof(ZpRanges)/sizeof(int[2]));
  case CLS_U_Zs: return u_in_ranges(c,ZsRanges,sizeof(ZsRanges)/sizeof(int[2]));
  case CLS_NL: return c=='\n'||c=='\r';
  case CLS_S: return xmlc_white_space(c);
  case CLS_I: return xmlc_base_char(c)||xmlc_ideographic(c)||c=='_'||c==':';
  case CLS_C: return in_class(c,CLS_I)||xmlc_digit(c)||xmlc_combining_char(c)||xmlc_extender(c)||c=='.'||c=='-';
  case CLS_W: return !(in_class(c,CLS_U_P)||in_class(c,CLS_U_Z)||in_class(c,CLS_U_C));
  default: assert(0);
  }
  return 0;
}


static int drv(int p,int c) {
  int p1,p2,cf,cl,cn,ret,m;
  assert(!P_IS(p,P_ERROR));
  m=new_memo(p,c);
  if(m!=-1) return M_RET(m);
  switch(P_TYP(p)) {
  case P_NOT_ALLOWED: case P_EMPTY: ret=notAllowed; break;
  case P_CHOICE: Choice(p,p1,p2); ret=choice(drv(p1,c),drv(p2,c)); break;
  case P_GROUP: Group(p,p1,p2); {int p11=group(drv(p1,c),p2); ret=nullable(p1)?choice(p11,drv(p2,c)):p11;} break;
  case P_ONE_OR_MORE: OneOrMore(p,p1); ret=group(drv(p1,c),choice(empty,p)); break;
  case P_EXCEPT: Except(p,p1,p2); ret=nullable(drv(p1,c))&&!nullable(drv(p2,c))?empty:notAllowed; break;
  case P_RANGE: Range(p,cf,cl); ret=cf<=c&&c<=cl?empty:notAllowed; break;
  case P_CLASS: Class(p,cn); ret=in_class(c,cn)?empty:notAllowed; break;
  case P_ANY: ret=empty; break;
  case P_CHAR: Char(p,cf); ret=c==cf?empty:notAllowed; break;
  default: ret=0; assert(0);
  }
  new_memo(p,c); M_SET(ret);
  accept_m();
  return ret;
}

int rx_check(char *rx) {(void)compile(rx); return !errors;}

int rx_match(char *rx,char *s,int n) {
  int p=compile(rx);
  if(!errors) {
    char *end=s+n;
    int u;
    for(;;) {
      if(p==notAllowed) return 0;
      if(s==end) return nullable(p);
      s+=u_get(&u,s);
      p=drv(p,u);
    }
  } else return 0;
}

int rx_rmatch(char *rx,char *s,int n) {
  int p=compile(rx);
  if(!errors) {
    char *end=s+n;
    int u;
    for(;;) {
      if(p==notAllowed) return 0;
      if(s==end) return nullable(p);
      s+=u_get(&u,s);
      if(xmlc_white_space(u)) u=' ';
      p=drv(p,u);
    }
  } else return 0;
}

int rx_cmatch(char *rx,char *s,int n) {
  int p=compile(rx);
  if(!errors) {
    char *end=s+n;
    int u;
    SKIP_SPACE: for(;;) {
      if(s==end) return nullable(p);
      s+=u_get(&u,s);
      if(!xmlc_white_space(u)) break;
    }
    for(;;) {
      if(p==notAllowed) return 0;
      if(xmlc_white_space(u)) { u=' ';
	p=drv(p,u);
	if(p==notAllowed) {
	  for(;;) {
	    if(s==end) return 1;
	    s+=u_get(&u,s);
	    if(!xmlc_white_space(u)) return 0;
	  }
	} else goto SKIP_SPACE;
      }
      p=drv(p,u);
      if(s==end) goto SKIP_SPACE;
      s+=u_get(&u,s);
    }
  } else return 0;
}

