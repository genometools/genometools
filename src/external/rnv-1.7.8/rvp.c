/* $Id: rvp.c,v 1.23 2004/11/09 11:29:31 dvd Exp $ */

/* validation pipe:
 synopsis

   rvp -qsdevh grammar.rnc

 reads from 0, writes to 1, 2 for grammar parse errors only, then redirected.
   -q switches to numerical error codes
   -s takes less space but more time
   -d plugs in an external type checker
   -e the argument is a Scheme program providing a datatype library
   -v displays version
   -h help message
 exit code: 0 on valid, non-zero on invalid

 protocol
  query ::= (start | quit | start-tag-open | attribute | start-tag-close | text | end-tag) z.
   quit ::= "quit".
   start ::= "start" [gramno].
   start-tag-open ::= "start-tag-open" patno name.
   attribute ::= "attribute" patno name value.
   start-tag-close :: = "start-tag-close" patno name.
   text ::= ("text"|"mixed") patno text.
   end-tag ::= "end-tag" patno name.
  response ::= (ok | er | error) z.
   ok ::= "ok" patno.
   er ::= "er" patno erno.
   error ::= "error" patno erno error.
  z ::= "\0" .

  conventions:
    last colon in name separates namespace uri and local part
    -q?er:error
    error==0 yields message 'protocol error' and happens when a query is not understood
    start assumes gramno=0 if the argument is omitted
*/

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <fcntl.h>  /*open,close*/
#include <sys/types.h>
#include UNISTD_H   /*open,read,close*/
#include <string.h> /*strerror*/
#include <setjmp.h>
#include <errno.h>
#include <assert.h>
#include "m.h"
#include "s.h"
#include "erbit.h"
#include "drv.h"
#include "rnl.h"
#include "rnv.h"
#include "dxl.h"
#include "dsl.h"
#include "er.h"

extern int rn_notAllowed, drv_compact, rx_compact;

#define ATT 0
#define ENT 1
#define MIX 2
#define QUIT 3
#define START 4
#define STC 5
#define STO 6
#define TXT 7
#define NKWD 8
char *kwdtab[NKWD]={
  "attribute",
  "end-tag",
  "mixed",
  "quit",
  "start",
  "start-tag-close",
  "start-tag-open",
  "text"
};

#define OK "ok %u"
#define ER "er %u"
#define ERROR "error %u"

#define LEN_B 1024

static FILE *nstderr;
static int explain=1, lasterr, *starts, n_st;
static int len_q,n_q; char *quebuf;
static int erp[2]; /* *erp to read error messages */
static jmp_buf IOER;

static void verror_handler(int erno,va_list ap) {
  lasterr=erno;
  rnv_default_verror_handler(erno&~ERBIT_RNV,ap);
}
static void verror_handler_rnv(int erno,va_list ap) {verror_handler(erno|ERBIT_RNV,ap);}

static int initialized=0;
static void init(void) {
  if(!initialized) {initialized=1;
    rnl_init();
    rnv_init(); rnv_verror_handler=&verror_handler_rnv;
    drv_add_dtl(DXL_URL,&dxl_equal,&dxl_allows);
    drv_add_dtl(DSL_URL,&dsl_equal,&dsl_allows);
    quebuf=(char*)m_alloc(len_q=LEN_B,sizeof(char));
  }
}

static int tok(int i) {
  for(;;) {
    switch(quebuf[i]) {
    case '\t': case '\n': case '\r': case ' ': break;
    default: return i;
    }
    ++i;
  }
}

static int endtok(int i) {
  for(;;) {
    switch(quebuf[i]) {
    case '\0': case '\t': case '\n': case '\r': case ' ': return i;
    default: break;
    }
    ++i;
  }
}

static void writeall(int fd,char *buf,int len) {
  int ofs=0;
  do {
    int n=write(fd,buf+ofs,len);
    if(n==-1) longjmp(IOER,1);
    ofs+=n; len-=n;
  } while(len);
}

static void resp(int ok,int patno,int prevno) {
  int len;
  static char buf[LEN_B];
  char *f=(char*)(ok?OK:explain?ERROR:ER);
  len=sprintf(buf,f,patno); assert(len<LEN_B);
  writeall(1,buf,len);
  if(!ok) {
    len=sprintf(buf," %u",lasterr); assert(len<LEN_B);
    writeall(1,buf,len);
    if(explain) {buf[0]=' '; writeall(1,buf,1);}
  }
  for(;;) { /* read always, write if verbose */
    len=read(erp[0],buf,LEN_B);
    if(len<0) {if(errno==EAGAIN) break; else longjmp(IOER,1);}
    if(len==0) break;
    if(!ok&&explain&&prevno!=rn_notAllowed) writeall(1,buf,len);
  }
  buf[0]='\0'; writeall(1,buf,1);
}

static int query(void) {
  int i,j,n,dn, kwd, patno,prevno, ok=0;
  char *name;
  n=0;
  for(;;) {
    if(n==n_q) {
      if(len_q-n_q<LEN_B) quebuf=(char*)m_stretch(quebuf,len_q=n_q+LEN_B,n_q,sizeof(char));
      dn=read(0,quebuf+n_q,LEN_B);
      if(dn<0) longjmp(IOER,1);
      if(dn==0) {errno=EIO; longjmp(IOER,1);}
      n_q+=dn;
    }
    if(quebuf[n++]=='\0') break;
  }

  j=endtok(i=tok(0));
  if((kwd=s_ntab(quebuf+i,j-i,kwdtab,NKWD))==QUIT) {resp(1,0,0); return 0;}
  switch(kwd) {
  case START:
    j=endtok((i=tok(j)));
    patno=0; while(i!=j) patno=patno*10+quebuf[i++]-'0';
    if(patno>=n_st) goto PROTER;
    ok=1; patno=starts[patno];
    break;
  case STO: case ATT: case STC: case TXT: case MIX: case ENT:
    j=endtok((i=tok(j))); if(i==j) goto PROTER;
    patno=0; do patno=patno*10+quebuf[i++]-'0'; while(i!=j);
    if(patno==0) goto PROTER; /* 0 is ERROR, not allowed */
    switch(kwd) {
    case STO: case ATT: case STC: case ENT:
      j=endtok((i=tok(j))); if(i==j||(kwd==ATT&&quebuf[j]=='\0')) goto PROTER;
      name=quebuf+i; quebuf[j]='\0';
      switch(kwd) {
      case STO: ok=rnv_start_tag_open(&patno,&prevno,name); break;
      case ATT: ok=rnv_attribute(&patno,&prevno,name,quebuf+j+1); break;
      case STC: ok=rnv_start_tag_close(&patno,&prevno,name); break;
      case ENT: ok=rnv_end_tag(&patno,&prevno,name); break;
      }
      break;
    case TXT: case MIX:
      if(quebuf[j]) ++j; i=j; while(quebuf[j]) ++j;
      ok=rnv_text(&patno,&prevno,quebuf+i,j-i,kwd==MIX);
      break;
    }
    break;

  case NKWD: PROTER: (*er_printf)("protocol error\n"); lasterr=0; patno=0; ok=0; break;
  default: assert(0);
  }
  resp(ok,patno,prevno);

  i=0; while(n!=n_q) quebuf[i++]=quebuf[n++]; n_q=i;
  return 1;
}

static void version(void) {(*er_printf)("rvp version %s\n",RVP_VERSION);}
static void usage(void) {(*er_printf)("usage: rvp {-[qs"
#if DXL_EXC
"d"
#endif
#if DSL_SCM
"e"
#endif
"vh?]} {schema.rnc}\n");}

int main(int argc,char **argv) {
  int i, ok;
  init();

  --argc;
  while(*(++argv)&&**argv=='-') {
    --argc; i=1;
    for(;;) {
      switch(*(*argv+i)) {
      case '\0': goto END_OF_OPTIONS;
      case 'h': case '?': usage(); return 0;
      case 'v': version(); break;
      case 's': drv_compact=1; rx_compact=1; break;
#if DXL_EXC
      case 'd': dxl_cmd=*(argv+1); if(*(argv+1)) ++argv; goto END_OF_OPTIONS;
#endif
#if DSL_SCM
      case 'e': dsl_ld(*(argv+1)); if(*(argv+1)) ++argv; goto END_OF_OPTIONS;
#endif
      case 'q': explain=0; break;
      default: (*er_printf)("unknown option '-%c'\n",*(*argv+i)); break;
      }
      ++i;
    }
    END_OF_OPTIONS:;
  }

  if(*argv==NULL) {usage(); return 1;}

  starts=(int*)m_alloc(argc,sizeof(int));
  ok=1; n_st=0;
  do {
    ok=(starts[n_st++]=rnl_fn(*(argv++)))&&ok;
  } while(*argv);
  if(ok) {
    int fd2;

    nstderr=stderr;
    if(setjmp(IOER)) {
      fprintf(nstderr,"%s\n",strerror(errno));
      return EXIT_FAILURE;
    }

    if((fd2=dup(2))==-1) longjmp(IOER,1);
    nstderr=fdopen(fd2,"w");
    if(pipe(erp)==-1||dup2(erp[1],2)==-1) longjmp(IOER,1);
    fcntl(erp[0],F_SETFL,O_NONBLOCK);
    setbuf(stderr,NULL);

    while(query());
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
