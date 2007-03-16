/* $Id: arx.c,v 1.37 2004/03/22 21:06:37 dvd Exp $ */
/* Regular Associations for XML

arx grammar:

arx = grammars route*
grammars = "grammars"  "{" type2string+ "}"
type2string =  type "=" literal
type = nmtoken
route = match|nomatch|valid|invalid
match = "=~" regexp "=>" type
nomatch = "!~" regexp "=>" type
valid = "valid" "{" rng "}" "=>" type
invalid = "!valid" "{" rng "}" "=>" type

literal=string in '"', '"' inside quoted by '\'
regexp=string in '/', '/' inside quoted by '\'
rng=relax ng compact syntax

comments start with # and continue till end of line
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include UNISTD_H
#include <fcntl.h>
#include <stdarg.h>
#include <errno.h>
#include <assert.h>
#include EXPAT_H
#include "u.h"
#include "m.h"
#include "s.h"
#include "xmlc.h"
#include "ht.h"
#include "erbit.h"
#include "rnl.h"
#include "rnv.h"
#include "rx.h"
#include "er.h"
#include "ary.h"

extern int rn_notAllowed;

/* rules */
#define VALID 1
#define INVAL 2
#define MATCH 3
#define NOMAT 4

#define LEN_2 16
#define LEN_R 64
#define LEN_S 64
#define S_AVG_SIZE 64

#define LEN_V 64

#define LEN_T 1024
#define LIM_T 65536

#define BUFSIZE 1024

static char *xml;
static int len_2,len_r,len_s,i_2,i_r,i_s;
static int (*t2s)[2],(*rules)[3];
static char *string; static struct hashtable ht_s;
static int path2abs;

/* arx parser */
static char *arxfn;
static int arxfd, i_b,len_b, cc, line,col,prevline,rnc, sym,len_v, errors;
static char buf[BUFSIZE];
static char *value;

/* xml validator */
static XML_Parser expat=NULL;
static int current,previous;
static int mixed=0;
static int ok,wf,any;
static char *text; static int len_txt;
static int n_txt;


static int add_s(char *s) {
  int len=strlen(s)+1,j;
  if(i_s+len>len_s) string=(char*)m_stretch(
    string,len_s=2*(i_s+len),i_s,sizeof(char));
  strcpy(string+i_s,s);
  if((j=ht_get(&ht_s,i_s))==-1) {
    ht_put(&ht_s,j=i_s);
    i_s+=len;
  }
  return j;
}

static int hash_s(int i) {return s_hval(string+i);}
static int equal_s(int s1,int s2) {return strcmp(string+s1,string+s2)==0;}

static void silent_verror_handler(int erno,va_list ap) {
  if(erno&ERBIT_DRV) rnv_default_verror_handler(erno,ap); /* low-level diagnostics */
}

static void windup(void);
static int initialized=0;
static void init(void) {
  if(!initialized) {initialized=1;
    rnl_init(); rnv_init();
    rnv_verror_handler=&silent_verror_handler;
    string=(char*)m_alloc(len_v=LEN_S*S_AVG_SIZE,sizeof(char));
    t2s=(int(*)[2])m_alloc(len_2=LEN_2,sizeof(int[2]));
    rules=(int(*)[3])m_alloc(len_r=LEN_R,sizeof(int[3]));
    ht_init(&ht_s,LEN_S,&hash_s,&equal_s);
    value=(char*)m_alloc(len_v=LEN_V,sizeof(char));
    text=(char*)m_alloc(len_txt=LEN_T,sizeof(char));
    windup();
  }
}

static void clear(void) {
  if(len_txt>LIM_T) {m_free(text); text=(char*)m_alloc(len_txt=LEN_T,sizeof(char));}
  ht_clear(&ht_s);
  windup();
}

static void windup(void) {
  text[n_txt=0]='\0';
  i_2=1; i_r=i_s=0;
}

/* parser */
#define SYM_EOF 0
#define SYM_GRMS 1
#define SYM_IDNT 2
#define SYM_LTRL 3
#define SYM_RGXP 4
#define SYM_RENG 5
#define SYM_MTCH 6
#define SYM_NMTC 7
#define SYM_VALD 8
#define SYM_NVAL 9
#define SYM_LCUR 10
#define SYM_RCUR 11
#define SYM_ASGN 12
#define SYM_INVL 13

static char *sym2str(int sym) {
  switch(sym) {
  case SYM_EOF: return "end of file";
  case SYM_GRMS: return "'grammars'";
  case SYM_IDNT: return "identifier";
  case SYM_LTRL: return "literal";
  case SYM_RGXP: return "regular expression";
  case SYM_RENG: return "Relax NG";
  case SYM_MTCH: return "'=~'";
  case SYM_NMTC: return "'!~'";
  case SYM_VALD: return "'valid'";
  case SYM_NVAL: return "'!valid'";
  case SYM_LCUR: return "'{'";
  case SYM_RCUR: return "'}'";
  case SYM_ASGN: return "'='";
  case SYM_INVL: return "invalid character";
  default: assert(0);
  }
  return NULL;
}

#define ARX_ER_IO 0
#define ARX_ER_SYN 1
#define ARX_ER_EXP 2
#define ARX_ER_REX 3
#define ARX_ER_RNG 4
#define ARX_ER_NOQ 5
#define ARX_ER_TYP 6

/* there is nothing in the grammar I need utf-8 processing for */
#define err(msg) (*er_vprintf)(msg"\n",ap)
static void verror_handler(int erno,va_list ap) {
  (*er_printf)("%s:%i:%i: error: ",arxfn,line,col);
  switch(erno) {
  case ARX_ER_IO: err("I/O error: %s"); break;
  case ARX_ER_SYN: err("syntax error"); break;
  case ARX_ER_EXP: err("%s expected, %s found"); break;
  case ARX_ER_REX: err("invalid regular expression"); break;
  case ARX_ER_RNG: err("invalid Relax NG grammar"); break;
  case ARX_ER_NOQ: err("unterminated literal or regular expression"); break;
  case ARX_ER_TYP: err("undeclared type '%s'"); break;
  }
}

static void error(int erno,...) {
  if(line!=prevline) {
    va_list ap; va_start(ap,erno); verror_handler(erno,ap); va_end(ap);
    prevline=line;
  }
  ++errors;
}

static void getcc(void) {
  for(;;) { int cc0=cc;
    if(i_b==len_b) {i_b=0; if((len_b=read(arxfd,buf,BUFSIZE))==-1) error(ARX_ER_IO,strerror(errno));}
    cc=i_b>=len_b?-1:((unsigned char*)buf)[i_b++];
    if(cc==-1) {if(cc0=='\n') break; else cc='\n';}
    if(cc=='\n' && cc0=='\r') continue;
    if(cc0=='\n' || cc0=='\r') {++line; col=0;} else ++col;
    break;
  }
}

static int nmtoken(int cc) {return cc>0x7F||xmlc_base_char(cc)||xmlc_digit(cc)||cc=='_'||cc=='.'||cc=='-'||cc==':';}
static int getid(void) {
  if(nmtoken(cc)) {
    int i=0;
    do {
      value[i++]=cc;
      if(i==len_v) value=(char*)m_stretch(value,len_v=2*i,i,sizeof(char));
      getcc();
    } while(nmtoken(cc));
    value[i]='\0';
    return 1;
  } else return 0;
}

static void getq(void) {
  int cq=cc;
  int i=0;
  for(;;) {
    getcc();
    if(cc==cq) {
      if(i!=0&&value[i-1]=='\\') --i; else {getcc(); break;}
    } else if(cc<' ') {error(ARX_ER_NOQ); break;}
    value[i++]=cc;
    if(i==len_v) value=(char*)m_stretch(value,len_v=2*i,i,sizeof(char));
  }
  value[i]='\0';
}

static void getrng(void) {
  int ircur=-1,i=0;
  int cc0;
  for(;;) {
    cc0=cc; getcc();
    if(cc=='}') ircur=i;
    else if(cc=='>') {if(cc0=='=') {getcc(); break;}} /* use => as terminator */
    else if(cc==-1) {error(ARX_ER_EXP,"=>",sym2str(SYM_EOF)); break;}
    value[i++]=cc;
    if(i==len_v) value=(char*)m_stretch(value,len_v=2*i,i,sizeof(char));
  }
  if(ircur==-1) {error(ARX_ER_EXP,sym2str(SYM_RCUR),sym2str(SYM_EOF)); ircur=0;}
  value[ircur]='\0';
}

static void getsym(void) {
  for(;;) {
    if(0<=cc&&cc<=' ') {getcc(); continue;}
    switch(cc) {
    case -1: sym=SYM_EOF; return;
    case '#': do getcc(); while(cc!='\n'&&cc!='\r'); getcc(); continue;
    case '{':
      if(sym==SYM_VALD||sym==SYM_NVAL) {
	getrng(); sym=SYM_RENG;
      } else {
	getcc(); sym=SYM_LCUR;
      }
      return;
    case '}': getcc(); sym=SYM_RCUR; return;
    case '!': getcc();
      if(cc=='~') {
	getcc(); sym=SYM_NMTC;
      } else {
	if(getid()) {
	  if(strcmp("valid",value)!=0) {error(ARX_ER_EXP,sym2str(SYM_NVAL),value);} sym=SYM_NVAL;
	} else {error(ARX_ER_SYN); sym=SYM_INVL;}
      }
      return;
    case '=': getcc();
      switch(cc) {
      case '~': getcc(); sym=SYM_MTCH; return;
      case '>': getcc(); if(sym!=SYM_RGXP) error(ARX_ER_SYN); continue;
      default: sym=SYM_ASGN; return;
      }
    case '"': getq(); sym=SYM_LTRL; return;
    case '/': getq(); sym=SYM_RGXP; return;
    default:
      if(getid()) {
	sym=strcmp("grammars",value)==0?SYM_GRMS
	 : strcmp("valid",value)==0?SYM_VALD:SYM_IDNT;
      } else {getcc(); error(ARX_ER_SYN); sym=SYM_INVL;}
      return;
    }
  }
}

static int chksym(int x) {
  if(sym!=x) {error(ARX_ER_EXP,sym2str(x),sym2str(sym)); return 0;}
  return 1;
}

static void chk_get(int x) {
  (void)chksym(x); getsym();
}

static int typ2str(void) {
  int i=i_2,typ=add_s(value);
  t2s[0][0]=typ; for(;;) if(t2s[--i][0]==typ) break;
  if(i==0) error(ARX_ER_TYP,value);
  return t2s[i][1];
}

static int arx(char *fn) {
  if((arxfd=open(arxfn=fn,O_RDONLY))==-1) {
    (*er_printf)("error (%s): %s\n",arxfn,strerror(errno));
    return 0;
  } else {
    errors=0;
    len_b=read(arxfd,buf,BUFSIZE); i_b=u_bom(buf,len_b);
    prevline=-1; line=1; col=0; rnc=0;
    cc=' '; getsym();
    chk_get(SYM_GRMS); chk_get(SYM_LCUR);
    do {
      if(i_2==len_2) t2s=(int(*)[2])m_stretch(t2s,len_2=i_2*2,i_2,sizeof(int[2]));
      if(chksym(SYM_IDNT)) t2s[i_2][0]=add_s(value);
      getsym();
      chk_get(SYM_ASGN);
      if(chksym(SYM_LTRL)) {
	if(path2abs) {
	  int len=strlen(arxfn)+strlen(value)+1;
	  if(len>len_v) {value=(char*)m_stretch(value,len,len_v,sizeof(char)); len_v=len;}
	  s_abspath(value,arxfn);
	}
	t2s[i_2][1]=add_s(value);
      }
      getsym();
      ++i_2;
    } while(sym==SYM_IDNT);
    chk_get(SYM_RCUR);
    for(;;) {
      if(i_r==len_r) rules=(int(*)[3])m_stretch(rules,len_r=i_r*2,i_r,sizeof(int[3]));
      switch(sym) {
      case SYM_MTCH: rules[i_r][0]=MATCH; goto REGEXP;
      case SYM_NMTC: rules[i_r][0]=NOMAT; goto REGEXP;
      REGEXP: getsym();
	if(chksym(SYM_RGXP)) {
	  if(!rx_check(value)) error(ARX_ER_REX);
	  rules[i_r][1]=add_s(value);
	}
	getsym();
	if(chksym(SYM_IDNT)) rules[i_r][2]=typ2str();
	goto NEXT;
      case SYM_VALD: rules[i_r][0]=VALID; goto RNG;
      case SYM_NVAL: rules[i_r][0]=INVAL; goto RNG;
      RNG: getsym();
	if(chksym(SYM_RENG)) {
	  char *rncfn=(char*)m_alloc(strlen(arxfn)+strlen("#rnc[]")+12,sizeof(char));
	  sprintf(rncfn,"%s#rnc[%i]",arxfn,rnc++);
	  if(!(rules[i_r][1]=rnl_s(rncfn,value,strlen(value)))) error(ARX_ER_RNG);
	  m_free(rncfn);
	}
	getsym();
	if(chksym(SYM_IDNT)) rules[i_r][2]=typ2str();
	goto NEXT;
      default: goto LAST;
      }
      NEXT: ++i_r; getsym();
    }
    LAST: chk_get(SYM_EOF);
    close(arxfd);
    return !errors;
  }
}

static void flush_text(void) {
  ok=rnv_text(&current,&previous,text,n_txt,mixed)&&ok;
  text[n_txt=0]='\0';
}

static void start_element(void *userData,const char *name,const char **attrs) {
  if(current!=rn_notAllowed) {
    mixed=1;
    flush_text();
    ok=rnv_start_tag(&current,&previous,(char*)name,(char**)attrs)&&ok;
    mixed=0; any=any||ary_isany(current);
  }
}

static void end_element(void *userData,const char *name) {
  if(current!=rn_notAllowed) {
    flush_text();
    ok=rnv_end_tag(&current,&previous,(char*)name)&&ok;
    mixed=1;
  }
}

static void characters(void *userData,const char *s,int len) {
  if(current!=rn_notAllowed) {
    int newlen_txt=n_txt+len+1;
    if(newlen_txt<=LIM_T&&LIM_T<len_txt) newlen_txt=LIM_T;
    else if(newlen_txt<len_txt) newlen_txt=len_txt;
    if(len_txt!=newlen_txt) text=(char*)m_stretch(text,len_txt=newlen_txt,n_txt,sizeof(char));
    memcpy(text+n_txt,s,len); n_txt+=len; text[n_txt]='\0'; /* '\0' guarantees that the text is bounded, and strto[ld] work for data */
  }
}

static void validate(int start,int fd) {
  void *buf; int len;
  previous=current=start;

  expat=XML_ParserCreateNS(NULL,':');
  XML_SetElementHandler(expat,&start_element,&end_element);
  XML_SetCharacterDataHandler(expat,&characters);
  ok=1; any=0;
  for(;;) {
    buf=XML_GetBuffer(expat,BUFSIZE);
    len=read(fd,buf,BUFSIZE);
    if(len<0) {
      (*er_printf)("error (%s): %s\n",xml,strerror(errno));
      wf=ok=0; break;
    }
    if(!XML_ParseBuffer(expat,len,len==0)) wf=ok=0;
    if(!ok||any||len==0) break;
  }
  XML_ParserFree(expat);
  return;
}

static void version(void) {(*er_printf)("arx version %s\n",ARX_VERSION);}
static void usage(void) {(*er_printf)("usage: arx {-[nvh?]} document.xml arx.conf {arx.conf}\n");}

int main(int argc,char **argv) {
  int fd;
  init();

  path2abs=1;
  while(*(++argv)&&**argv=='-') {
    int i=1;
    for(;;) {
      switch(*(*argv+i)) {
      case '\0': goto END_OF_OPTIONS;
      case 'h': case '?': usage(); return 1;
      case 'n': path2abs=0; break;
      case 'v': version(); break;
      default: (*er_printf)("unknown option '-%c'\n",*(*argv+i)); break;
      }
      ++i;
    }
    END_OF_OPTIONS:;
  }

  if(!(*(argv)&&*(argv+1))) {usage(); return 1;}

  xml=*(argv++); if((wf=(fd=open(xml,O_RDONLY))!=-1)) close(fd);
  do {
    if(arx(*(argv++))) {
      int i;
      for(i=0;i!=i_r;++i) {
	switch(rules[i][0]) {
	case VALID: if((ok=wf)) {validate(rules[i][1],fd=open(xml,O_RDONLY)); close(fd);} break;
	case INVAL: if((ok=wf)) {validate(rules[i][1],fd=open(xml,O_RDONLY)); close(fd); ok=wf&&!ok;} break;
	case MATCH: ok=rx_match(string+rules[i][1],xml,strlen(xml)); break;
	case NOMAT: ok=!rx_match(string+rules[i][1],xml,strlen(xml)); break;
	default: assert(0);
	}
	if(ok) {
	  printf("%s\n",string+rules[i][2]);
	  return EXIT_SUCCESS;
	}
      }
    }
    clear();
  } while(*argv);
  return EXIT_FAILURE;
}
