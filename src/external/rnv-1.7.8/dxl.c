/* $Id: dxl.c,v 1.13 2004/01/28 12:34:12 dvd Exp $ */

#include <stdlib.h>
#include "dxl.h"

char *dxl_cmd=NULL;

#if DXL_EXC

#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include "m.h"
#include "er.h"

int dxl_allows(char *typ,char *ps,char *s,int n) {
  int pid,status;

  if(!dxl_cmd) return 0;
  if((pid=fork())==0) {
    char **argv; int argc;
    char *p; int arg, i;

    argc=5; p=ps; arg=0;
    for(;;) {
      if(*p=='\0') {
	if(arg) {arg=0; ++argc;} else break;
      } else arg=1;
      ++p;
    }
    argv=(char**)m_alloc(argc,sizeof(char*));
    argv[--argc]=NULL;
    argv[--argc]=(char*)m_alloc(n+1,sizeof(char)); argv[argc][n]='\0'; strncpy(argv[argc],s,n);
    argv[0]=dxl_cmd; argv[1]="allows"; argv[2]=typ;
    i=3; if(i<argc) {
      for(;;) {
	argv[i++]=ps;
	if(i==argc) break;
	while(*(ps++));
      }
    }
    execv(dxl_cmd,argv);
    (*er_printf)("dxl: cannot execute %s: %s\n",dxl_cmd,strerror(errno));
 } else if(pid>0) {
    wait(&status);
    return !WEXITSTATUS(status);
  }
  (*er_printf)("dxl: %s\n",strerror(errno));
  return 0;
}

int dxl_equal(char *typ,char *val,char *s,int n) {
  int pid,status;

  if(!dxl_cmd) return 0;
  if((pid=fork())==0) {
    char *argv[]={NULL,"equal",NULL,NULL,NULL,NULL};
    argv[0]=dxl_cmd; argv[2]=typ; argv[3]=val;

   argv[4]=(char*)m_alloc(n+1,sizeof(char)); argv[4][n]='\0'; strncpy(argv[4],s,n);
     execvp(dxl_cmd,argv);
    (*er_printf)("dxl: cannot execute %s\n",dxl_cmd,strerror(errno));
  } else if(pid>0) {
    wait(&status);
    return !WEXITSTATUS(status);
  }
  (*er_printf)("dxl: %s\n",strerror(errno));
  return 0;
}

#else

int dxl_allows(__attribute__ ((unused)) char *typ,
               __attribute__ ((unused)) char *ps,
               __attribute__ ((unused)) char *s,
               __attribute__ ((unused)) int n)
{
  return 0;
}
int dxl_equal(__attribute__ ((unused)) char *typ,
              __attribute__ ((unused)) char *val,
              __attribute__ ((unused)) char *s,
              __attribute__ ((unused)) int n)
{
  return 0;
}

#endif
