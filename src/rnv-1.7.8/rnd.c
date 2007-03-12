/* $Id: rnd.c,v 1.33 2004/02/25 00:00:32 dvd Exp $ */

#include <stdlib.h>
#include <assert.h>
#include "m.h"
#include "rn.h"
#include "rnx.h"
#include "ll.h"
#include "er.h"
#include "rnd.h"

#define LEN_F RND_LEN_F

static int len_f,n_f;
static int *flat;
static int errors;

#define err(msg) (*er_vprintf)("error: "msg"\n",ap)
void rnd_default_verror_handler(int er_no,va_list ap) {
  switch(er_no) {
  case RND_ER_LOOPST: err("loop in start pattern"); break;
  case RND_ER_LOOPEL: err("loop in pattern for element '%s'"); break;
  case RND_ER_CTYPE: err("content of element '%s' does not have a content-type"); break;
  case RND_ER_BADSTART: err("bad path in start pattern"); break;
  case RND_ER_BADMORE: err("bad path before '*' or '+' in element '%s'"); break;
  case RND_ER_BADEXPT: err("bad path after '-' in element '%s'"); break;
  case RND_ER_BADLIST: err("bad path after 'list' in element '%s'"); break;
  case RND_ER_BADATTR: err("bad path in attribute '%s' of element '%s'"); break;
  default: assert(0);
  }
}

void (*rnd_verror_handler)(int er_no,va_list ap)=&rnd_default_verror_handler;

static int initialized=0;
void rnd_init(void) {
  if(!initialized) {
    rn_init();
    initialized=1;
  }
}

void rnd_clear(void) {}

static void error(int er_no,...) {
  va_list ap; va_start(ap,er_no); (*rnd_verror_handler)(er_no,ap); va_end(ap);
  ++errors;
}

static int de(int p) {
  int p0=p,p1;
  RN_P_CHK(p,RN_P_REF);
  for(;;) {
    rn_Ref(p,p1);
    if(!RN_P_IS(p1,RN_P_REF)||p1==p0) break;
    p=p1;
  }
  return p1;
}

static void flatten(int p) { if(!rn_marked(p)) {flat[n_f++]=p; rn_mark(p);}}

static void deref(int start) {
  int p,p1,p2,nc,i,changed;

  flat=(int*)m_alloc(len_f=LEN_F,sizeof(int)); n_f=0;
  if(RN_P_IS(start,RN_P_REF)) start=de(start);
  flatten(start);

  i=0;
  do {
    p=flat[i++];
    switch(RN_P_TYP(p)) {
    case RN_P_NOT_ALLOWED: case RN_P_EMPTY: case RN_P_TEXT: case RN_P_DATA: case RN_P_VALUE:
      break;

    case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
    case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); goto BINARY;
    case RN_P_GROUP: rn_Group(p,p1,p2); goto BINARY;
    case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); goto BINARY;
    BINARY:
      changed=0;
      if(RN_P_IS(p1,RN_P_REF)) {p1=de(p1); changed=1;}
      if(RN_P_IS(p2,RN_P_REF)) {p2=de(p2); changed=1;}
      if(changed) {rn_del_p(p); rn_pattern[p+1]=p1; rn_pattern[p+2]=p2; rn_add_p(p);}
      if(n_f+2>len_f) flat=(int*)m_stretch(flat,len_f=2*(n_f+2),n_f,sizeof(int));
      flatten(p1); flatten(p2);
      break;

    case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); goto UNARY;
    case RN_P_LIST: rn_List(p,p1); goto UNARY;
    case RN_P_ATTRIBUTE: rn_Attribute(p,nc,p1); goto UNARY;
    case RN_P_ELEMENT: rn_Element(p,nc,p1); goto UNARY;
    UNARY:
      changed=0;
      if(RN_P_IS(p1,RN_P_REF)) {p1=de(p1); changed=1;}
      if(changed) {rn_del_p(p); rn_pattern[p+1]=p1; rn_add_p(p);}
      if(n_f+1>len_f) flat=(int*)m_stretch(flat,len_f=2*(n_f+1),n_f,sizeof(int));
      flatten(p1);
      break;

    case RN_P_REF: /* because of a loop, but will be handled in rnd_loops */
      break;

    default:
      assert(0);
    }
  } while(i!=n_f);
  for(i=0;i!=n_f;++i) rn_unmark(flat[i]);
}

static int loop(int p) {
  int nc,p1,p2,ret=1;
  if(rn_marked(p)) return 1;
  rn_mark(p);
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED: case RN_P_EMPTY: case RN_P_TEXT: case RN_P_DATA: case RN_P_VALUE:
  case RN_P_ELEMENT:
    ret=0; break;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); goto BINARY;
  case RN_P_GROUP: rn_Group(p,p1,p2); goto BINARY;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); goto BINARY;
  BINARY:
    ret=loop(p1)||loop(p2); break;

  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); goto UNARY;
  case RN_P_LIST: rn_List(p,p1); goto UNARY;
  case RN_P_ATTRIBUTE:  rn_Attribute(p,nc,p1); goto UNARY;
  UNARY:
    ret=loop(p1); break;

  case RN_P_REF: ret=1; break;

  default: assert(0);
  }
  rn_unmark(p);
  return ret;
}

static void loops(void) {
  int i=0,p=flat[i],nc=-1,p1;
  for(;;) {
    if(loop(p)) {
      if(i==0) error(RND_ER_LOOPST); else {
	char *s=rnx_nc2str(nc);
	error(RND_ER_LOOPEL,s);
	m_free(s);
      }
    }
    for(;;) {++i;
      if(i==n_f) return;
      p=flat[i];
      if(RN_P_IS(p,RN_P_ELEMENT)) {
	rn_Element(p,nc,p1); p=p1;
	break;
      }
    }
  }
}

static void ctype(int p) {
  int p1,p2,nc;
  if(!rn_contentType(p)) {
    switch(RN_P_TYP(p)) {
    case RN_P_NOT_ALLOWED: rn_setContentType(p,RN_P_FLG_CTE,0); break;
    case RN_P_EMPTY: rn_setContentType(p,RN_P_FLG_CTE,0); break;
    case RN_P_TEXT: rn_setContentType(p,RN_P_FLG_CTC,0); break;
    case RN_P_CHOICE: rn_Choice(p,p1,p2); ctype(p1); ctype(p2);
      rn_setContentType(p,rn_contentType(p1),rn_contentType(p2)); break;
    case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); ctype(p1); ctype(p2);
      if(rn_groupable(p1,p2)) rn_setContentType(p,rn_contentType(p1),rn_contentType(p2)); break;
    case RN_P_GROUP: rn_Group(p,p1,p2); ctype(p1); ctype(p2);
      if(rn_groupable(p1,p2)) rn_setContentType(p,rn_contentType(p1),rn_contentType(p2)); break;
    case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); ctype(p1);
      if(rn_groupable(p1,p1)) rn_setContentType(p,rn_contentType(p1),0); break;
    case RN_P_LIST: rn_setContentType(p,RN_P_FLG_CTS,0); break;
    case RN_P_DATA: rn_setContentType(p,RN_P_FLG_CTS,0); break;
    case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); ctype(p1); ctype(p2);
      if(rn_contentType(p2)) rn_setContentType(p,RN_P_FLG_CTS,0); break;
    case RN_P_VALUE: rn_setContentType(p,RN_P_FLG_CTS,0); break;
    case RN_P_ATTRIBUTE: rn_Attribute(p,nc,p1); ctype(p1);
      if(rn_contentType(p1)) rn_setContentType(p,RN_P_FLG_CTE,0); break;
    case RN_P_ELEMENT: rn_setContentType(p,RN_P_FLG_CTC,0); break;
    default: assert(0);
    }
  }
}

static void ctypes(void) {
  int i,p,p1,nc;
  for(i=0;i!=n_f;++i) {
    p=flat[i];
    if(RN_P_IS(p,RN_P_ELEMENT)) {
      rn_Element(p,nc,p1);
      ctype(p1);
      if(!rn_contentType(p1)) {
	char *s=rnx_nc2str(nc);
	error(RND_ER_CTYPE,s);
	m_free(s);
      }
    }
  }
}

static int bad_start(int p) {
  int p1,p2;
  switch(RN_P_TYP(p)) {
  case RN_P_EMPTY: case RN_P_TEXT:
  case RN_P_INTERLEAVE: case RN_P_GROUP: case RN_P_ONE_OR_MORE:
  case RN_P_LIST: case RN_P_DATA: case RN_P_DATA_EXCEPT: case RN_P_VALUE:
  case RN_P_ATTRIBUTE:
    return 1;
  case RN_P_NOT_ALLOWED:
  case RN_P_ELEMENT:
    return 0;
  case RN_P_CHOICE: rn_Choice(p,p1,p2);
    return bad_start(p1)||bad_start(p2);
  default: assert(0);
  }
  return 1;
}

static int bad_data_except(int p) {
  int p1,p2;
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED:
  case RN_P_VALUE: case RN_P_DATA:
    return 0;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_DATA_EXCEPT: rn_Choice(p,p1,p2); goto BINARY;
  BINARY: return bad_data_except(p1)||bad_data_except(p2);

  case RN_P_EMPTY: case RN_P_TEXT:
  case RN_P_INTERLEAVE: case RN_P_GROUP: case RN_P_ONE_OR_MORE:
  case RN_P_LIST:
  case RN_P_ATTRIBUTE: case RN_P_ELEMENT:
    return 1;
  default: assert(0);
  }
  return 1;
}

static int bad_one_or_more(int p,int in_group) {
  int nc,p1,p2;
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED: case RN_P_EMPTY: case RN_P_TEXT:
  case RN_P_DATA: case RN_P_VALUE:
  case RN_P_ELEMENT:
    return 0;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); in_group=1; goto BINARY;
  case RN_P_GROUP: rn_Group(p,p1,p2); in_group=1; goto BINARY;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); goto BINARY;
  BINARY: return  bad_one_or_more(p1,in_group)||bad_one_or_more(p2,in_group);

  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); goto UNARY;
  case RN_P_LIST: rn_List(p,p1); goto UNARY;
  case RN_P_ATTRIBUTE: if(in_group) return 1;
    rn_Attribute(p,nc,p1); goto UNARY;
  UNARY: return  bad_one_or_more(p1,in_group);
  default: assert(0);
  }
  return 1;
}

static int bad_list(int p) {
  int p1,p2;
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED: case RN_P_EMPTY:
  case RN_P_DATA: case RN_P_VALUE:
    return 0;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_GROUP: rn_Group(p,p1,p2); goto BINARY;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); goto BINARY;
  BINARY: return bad_list(p1)||bad_list(p2);

  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); goto UNARY;
  case RN_P_LIST: rn_List(p,p1); goto UNARY;
  UNARY: return bad_list(p1);

  case RN_P_TEXT:
  case RN_P_INTERLEAVE:
  case RN_P_ATTRIBUTE:
  case RN_P_ELEMENT:
    return 1;
  default: assert(0);
  }
  return 1;
}

static int bad_attribute(int p) {
  int p1,p2;
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED: case RN_P_EMPTY: case RN_P_TEXT:
  case RN_P_DATA: case RN_P_VALUE:
    return 0;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); goto BINARY;
  case RN_P_GROUP: rn_Group(p,p1,p2); goto BINARY;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2); goto BINARY;
  BINARY: return bad_attribute(p1)||bad_attribute(p2);


  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); goto UNARY;
  case RN_P_LIST: rn_List(p,p1); goto UNARY;
  UNARY: return bad_attribute(p1);

  case RN_P_ATTRIBUTE: case RN_P_ELEMENT:
    return 1;
  default: assert(0);
  }
  return 1;
}

static void path(int p,int nc) {
  int p1,p2,nc1;
  switch(RN_P_TYP(p)) {
  case RN_P_NOT_ALLOWED: case RN_P_EMPTY: case RN_P_TEXT:
  case RN_P_DATA: case RN_P_VALUE:
  case RN_P_ELEMENT:
    break;

  case RN_P_CHOICE: rn_Choice(p,p1,p2); goto BINARY;
  case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); goto BINARY;
  case RN_P_GROUP: rn_Group(p,p1,p2); goto BINARY;
  case RN_P_DATA_EXCEPT: rn_DataExcept(p,p1,p2);
    if(bad_data_except(p2)) {char *s=rnx_nc2str(nc); error(RND_ER_BADEXPT,s); m_free(s);}
    goto BINARY;
  BINARY: path(p1,nc); path(p2,nc); break;

  case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1);
    if(bad_one_or_more(p1,0)) {char *s=rnx_nc2str(nc); error(RND_ER_BADMORE,s); m_free(s);}
    goto UNARY;
  case RN_P_LIST: rn_List(p,p1);
    if(bad_list(p1)) {char *s=rnx_nc2str(nc); error(RND_ER_BADLIST,s); m_free(s);}
    goto UNARY;
  case RN_P_ATTRIBUTE: rn_Attribute(p,nc1,p1);
    if(bad_attribute(p1)) {char *s=rnx_nc2str(nc),*s1=rnx_nc2str(nc1); error(RND_ER_BADATTR,s1,s); m_free(s1); m_free(s);}
    goto UNARY;
  UNARY: path(p1,nc); break;

  default: assert(0);
  }
}

static void paths(void) {
  int i,p,p1,nc;
  if(bad_start(flat[0])) error(RND_ER_BADSTART);
  for(i=0;i!=n_f;++i) {
    p=flat[i];
    if(RN_P_IS(p,RN_P_ELEMENT)) {
      rn_Element(p,nc,p1);
      path(p1,nc);
    }
  }
}

static void restrictions(void) {
  loops(); if(errors) return; /* loops can cause endless loops in subsequent calls */
  ctypes();
  paths();
}

static void nullables(void) {
  int i,p,p1,p2,changed;
  do {
    changed=0;
    for(i=0;i!=n_f;++i) {
      p=flat[i];
      if(!rn_nullable(p)) {
	switch(RN_P_TYP(p)) {
	case RN_P_NOT_ALLOWED:
	case RN_P_DATA: case RN_P_DATA_EXCEPT: case RN_P_VALUE: case RN_P_LIST:
	case RN_P_ATTRIBUTE: case RN_P_ELEMENT:
	  break;

	case RN_P_CHOICE: rn_Choice(p,p1,p2); rn_setNullable(p,rn_nullable(p1)||rn_nullable(p2)); break;
	case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); rn_setNullable(p,rn_nullable(p1)&&rn_nullable(p2)); break;
	case RN_P_GROUP: rn_Group(p,p1,p2);  rn_setNullable(p,rn_nullable(p1)&&rn_nullable(p2)); break;

	case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); rn_setNullable(p,rn_nullable(p1)); break;

	default: assert(0);
	}
	changed=changed||rn_nullable(p);
      }
    }
  } while(changed);
}

static void cdatas(void) {
  int i,p,p1,p2,changed;
  do {
    changed=0;
    for(i=0;i!=n_f;++i) {
      p=flat[i];
      if(!rn_cdata(p)) {
	switch(RN_P_TYP(p)) {
	case RN_P_NOT_ALLOWED: case RN_P_EMPTY:
	case RN_P_ATTRIBUTE: case RN_P_ELEMENT:
	  break;

	case RN_P_CHOICE: rn_Choice(p,p1,p2); rn_setCdata(p,rn_cdata(p1)||rn_cdata(p2)); break;
	case RN_P_INTERLEAVE: rn_Interleave(p,p1,p2); rn_setCdata(p,rn_cdata(p1)||rn_cdata(p2)); break;
	case RN_P_GROUP: rn_Group(p,p1,p2);  rn_setCdata(p,rn_cdata(p1)||rn_cdata(p2)); break;

	case RN_P_ONE_OR_MORE: rn_OneOrMore(p,p1); rn_setCdata(p,rn_cdata(p1)); break;

	default: assert(0);
	}
	changed=changed||rn_cdata(p);
      }
    }
  } while(changed);
}

static void traits(void) {
  nullables();
  cdatas();
}

static int release(void) {
  int start=flat[0];
  m_free(flat); flat=NULL;
  return start;
}

int rnd_fixup(int start) {
  errors=0; deref(start);
  if(!errors) {restrictions(); if(!errors) traits();}
  start=release(); return errors?0:start;
}
