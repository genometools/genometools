/* $Id: rnx.c,v 1.24 2004/02/25 00:00:32 dvd Exp $ */

#include <stdlib.h> /*NULL*/
#include <string.h> /*strcat*/
#include "m.h"
#include "s.h"
#include "rn.h"
#include "ll.h"
#include "rnx.h"

#define LEN_EXP RNX_LEN_EXP
#define LIM_EXP RNX_LIM_EXP

int rnx_n_exp,*rnx_exp=NULL;
static int len_exp;

static int initialized=0;
void rnx_init(void) {
  if(!initialized) { initialized=1;
    rnx_exp=(int*)m_alloc(len_exp=LEN_EXP,sizeof(int));
  }
}

void rnx_clear(void) {}

static void expected(int p,int first,int req) {
  int p1,p2,px=0,i;
  if(req && rn_nullable(p)) return;
  switch(RN_P_TYP(p)) {
  case RN_P_ERROR: break;
  case RN_P_NOT_ALLOWED: break;
  case RN_P_EMPTY: break;
  case RN_P_TEXT: px=p; break;
  case RN_P_CHOICE: rn_Choice(p,p1,p2);
    expected(p1,first,req); expected(p2,first,req); break;
  case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2);
    expected(p1,first,req); expected(p2,first,req); break;
  case RN_P_GROUP: rn_Group(p,p1,p2);
    expected(p1,first,req); expected(p2,first&&rn_nullable(p1),req); break;
  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); expected(p1,first,req); break;
  case RN_P_LIST: rn_List(p,p1); expected(p1,first,req); break;
  case RN_P_DATA: px=p; break;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2);
    expected(p1,first,req); break;
  case RN_P_VALUE: px=p; break;
  case RN_P_ATTRIBUTE: px=p; break;
  case RN_P_ELEMENT: px=p; break;
  case RN_P_AFTER: rn_After(p,p1,p2);
    expected(p1,first,req); if(rn_nullable(p1)) px=p; break;
  case RN_P_REF: break;
  default: assert(0);
  }
  if(px&&(first||RN_P_IS(px,RN_P_ATTRIBUTE))) {
    for(i=0;i!=rnx_n_exp;++i) {
      if(rnx_exp[i]==px) {px=0; break;}
    }
    if(px) {
      if(rnx_n_exp==len_exp) rnx_exp=(int*)m_stretch(rnx_exp,len_exp=2*rnx_n_exp,rnx_n_exp,sizeof(int));
      rnx_exp[rnx_n_exp++]=px;
    }
  }
}
void rnx_expected(int p,int req) {
  if(req) {
    if(len_exp>LIM_EXP) {
      m_free(rnx_exp);
      rnx_exp=(int*)m_alloc(len_exp=LIM_EXP,sizeof(int));
    }
    rnx_n_exp=0;
  }
  expected(p,1,req);
}

char *rnx_p2str(int p) {
  char *s=NULL,*s1;
  int dt,ps,val,nc,p1;
  switch(RN_P_TYP(p)) {
  case RN_P_ERROR: s=s_clone("error"); break;
  case RN_P_NOT_ALLOWED: s=s_clone("notAllowed"); break;
  case RN_P_EMPTY: s=s_clone("empty"); break;
  case RN_P_TEXT: s=s_clone("text"); break;
  case RN_P_CHOICE: s=s_clone("choice (|)"); break;
  case RN_P_INTERLEAVE: s=s_clone("interleave (&)"); break;
  case RN_P_GROUP: s=s_clone("group (,)"); break;
  case RN_P_ONE_OR_MORE: s=s_clone("one or more (+)"); break;
  case RN_P_LIST: s=s_clone("list"); break;
  case RN_P_DATA: rn_Data(p,dt,ps);
    s1=rnx_nc2str(dt);
    s=(char*)m_alloc(strlen("data ")+1+strlen(s1),sizeof(char));
    strcpy(s,"data "); strcat(s,s1);
    m_free(s1);
    break;
  case RN_P_DATA_EXCEPT: s=s_clone("dataExcept (-)");  break;
  case RN_P_VALUE: rn_Value(p,dt,val);
    s1=rnx_nc2str(dt);
    s=(char*)m_alloc(strlen("value \"\" ")+1+strlen(s1)+strlen(rn_string+val),sizeof(char));
    strcpy(s,"value "); strcat(s,s1); strcat(s," \""); strcat(s,rn_string+val); strcat(s,"\"");
    m_free(s1);
    break;
  case RN_P_ATTRIBUTE: rn_Attribute(p,nc,p1);
    s1=rnx_nc2str(nc);
    s=(char*)m_alloc(strlen("attribute ")+1+strlen(s1),sizeof(char));
    strcpy(s,"attribute "); strcat(s,s1);
    m_free(s1);
    break;
  case RN_P_ELEMENT: rn_Element(p,nc,p1);
    s1=rnx_nc2str(nc);
    s=(char*)m_alloc(strlen("element ")+1+strlen(s1),sizeof(char));
    strcpy(s,"element "); strcat(s,s1);
    m_free(s1);
    break;
  case RN_P_REF: s=s_clone("ref"); break;
  case RN_P_AFTER: s=s_clone("after"); break;
  default: assert(0);
  }
  return s;
}

char *rnx_nc2str(int nc) {
  char *s=NULL,*s1,*s2;
  int nc1,nc2,uri,name;
  switch(RN_NC_TYP(nc)) {
  case RN_NC_ERROR: s=s_clone("?"); break;
  case RN_NC_NSNAME:
    rn_NsName(nc,uri);
    s=(char*)m_alloc(strlen(rn_string+uri)+3,sizeof(char));
    strcpy(s,rn_string+uri); strcat(s,":*");
    break;
  case RN_NC_QNAME:
    rn_QName(nc,uri,name);
    s=(char*)m_alloc(strlen(rn_string+uri)+strlen(rn_string+name)+2,sizeof(char));
    strcpy(s,rn_string+uri); strcat(s,"^"); strcat(s,rn_string+name);
    break;
  case RN_NC_ANY_NAME: s=s_clone("*"); break;
  case RN_NC_EXCEPT:
    rn_NameClassExcept(nc,nc1,nc2);
    s1=rnx_nc2str(nc1); s2=rnx_nc2str(nc2);
    s=(char*)m_alloc(strlen(s1)+strlen(s2)+2,sizeof(char));
    strcpy(s,s1); strcat(s,"-"); strcat(s,s2);
    m_free(s1); m_free(s2);
    break;
  case RN_NC_CHOICE:
    rn_NameClassChoice(nc,nc1,nc2);
    s1=rnx_nc2str(nc1); s2=rnx_nc2str(nc2);
    s=(char*)m_alloc(strlen(s1)+strlen(s2)+2,sizeof(char));
    strcpy(s,s1); strcat(s,"|"); strcat(s,s2);
    m_free(s1); m_free(s2);
    break;
  case RN_NC_DATATYPE:
    rn_Datatype(nc,uri,name);
    s=(char*)m_alloc(strlen(rn_string+uri)+strlen(rn_string+name)+2,sizeof(char));
    strcpy(s,rn_string+uri); strcat(s,"^"); strcat(s,rn_string+name);
    break;
  default: assert(0);
  }
  return s;
}
