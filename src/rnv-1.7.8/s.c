/* $Id: s.c,v 1.4 2005/01/05 09:46:25 dvd Exp $ */

#include <string.h> /*strcpy,strlen*/
#include <assert.h>
#include "xmlc.h"
#include "m.h"
#include "s.h"

int s_cmpn(char *s1,char *s2,int n2) {
  char *end=s2+n2;
  for(;;++s1,++s2) {
    if(s2==end) return *s1;
    if(*s1=='\0') return -*s2;
    if(*s1!=*s2) return *s1-*s2;
  }
}

int s_tokcmpn(char *s1,char *s2,int n2) {
  char *end2=s2+n2;
 /* all white space characters are one byte long */
  while(xmlc_white_space(*s1)) ++s1;
  while(s2!=end2&&xmlc_white_space(*s2)) ++s2;
  for(;;) {
    if(s2==end2) {
      while(xmlc_white_space(*s1)) ++s1;
      return *s1;
    }
    if(*s1=='\0') {
      while(s2!=end2&&xmlc_white_space(*s2)) ++s2;
      return s2==end2?0:-*s2;
    }
    if(xmlc_white_space(*s1)&&xmlc_white_space(*s2)) {
      do ++s1; while(xmlc_white_space(*s1));
      do ++s2; while(s2!=end2&&xmlc_white_space(*s2));
    } else {
      if(*s1!=*s2) return *s1-*s2;
      ++s1; ++s2;
    }
  }
}

int s_hval(char *s) {
  int h=0;
  while(*s) h=h*31+*(s++);
  return h;
}

char *s_clone(char *s) {
  return strcpy((char*)m_alloc(strlen(s)+1,sizeof(char)),s);
}

char *s_abspath(char *r,char *b) {
  if(*r!='/') {
    char *c=b,*sep=(char*)0;
    for(;;) {if(!(*c)) break; if(*c++=='/') sep=c;}
    if(sep) {
      char *p=r,*q;
      while(*p++); q=p+(sep-b);
      do *(--q)=*(--p); while(p!=r);
      while(b!=sep) *r++=*b++;
    }
  }
  return r;
}

int s_tab(char *s,char *tab[],int size) {return s_ntab(s,strlen(s),tab,size);}
int s_ntab(char *s,int len,char *tab[],int size) {
  int n=0,m=size-1,i,cmp;
  for(;;) {
    if(n>m) return size;
    i=(n+m)/2;
    if((cmp=s_cmpn(tab[i],s,len))==0) return i; else {if(cmp>0) m=i-1; else n=i+1;}
  }
}

void s_test() {
  assert(s_cmpn("","",0)==0);
  assert(s_cmpn("/xyz","/xyz",4)==0);
  assert(s_cmpn("xyz","yz",2)<0);
  assert(s_cmpn("xyz","xxyz",4)>0);

  { char r[256];
    s_abspath(strcpy(r,"/x"),"/y");
    assert(strcmp(r,"/x")==0);
    s_abspath(strcpy(r,"x"),"/y");
    assert(strcmp(r,"/x")==0);
    s_abspath(strcpy(r,"x"),"/y/");
    assert(strcmp(r,"/y/x")==0);
    s_abspath(strcpy(r,"x"),"y/");
    assert(strcmp(r,"y/x")==0);
    s_abspath(strcpy(r,"x"),"y");
    assert(strcmp(r,"x")==0);
    s_abspath(strcpy(r,""),"y");
    assert(strcmp(r,"")==0);
  }

  assert(s_tokcmpn("","",0)==0);
  assert(s_tokcmpn(""," ",1)==0);
  assert(s_tokcmpn("A","A",1)==0);
  assert(s_tokcmpn(" A   B","A B  ",5)==0);
  assert(s_tokcmpn("AB","A B",3)>0);
  assert(s_tokcmpn("","A",1)<0);
}
