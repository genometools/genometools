/* $Id: ary.c,v 1.4 2004/02/25 00:00:32 dvd Exp $ */

#include "rn.h"
#include "ary.h"

/*
ary_isany::Pattern->Bool
ary_isany p =
  let
    isanycontent
      p@(OneOrMore
	  (Choice
	    (Choice
	      (Element AnyName p1)
	      (Attribute AnyName Text))
	    Text)) = p == p1
    isanycontent _ = False
    isanymixed (OneOrMore (Choice (Element AnyName p1) Text)) = isanycontent p1
    isanymixed _ = False
  in
     case p of
       (After p1 Empty) -> isanymixed p1
       (After p1 p2) -> isanymixed p1 && ary_isany p2
       _ -> False
*/

static int isanycont(int p) {
  int p0,nc,p1,p2,i,res,flat[3];
  p0=p; if(!RN_P_IS(p0,RN_P_ONE_OR_MORE)) return 0;
  rn_OneOrMore(p0,p1);
  p0=p1; if(!RN_P_IS(p0,RN_P_CHOICE)) return 0;
  rn_Choice(p0,p1,p2); flat[0]=p2;
  p0=p1; if(!RN_P_IS(p0,RN_P_CHOICE)) return 0;
  rn_Choice(p0,p1,p2); flat[1]=p1; flat[2]=p2;
  res=0;
  for(i=0;i!=3;++i) {
    p0=flat[i];
    switch(RN_P_TYP(p0)) {
    case RN_P_ELEMENT: rn_Element(p0,nc,p1);
      if(!(RN_NC_IS(nc,RN_NC_ANY_NAME)&&p==p1)) return 0;
      res|=1; break;
    case RN_P_ATTRIBUTE: rn_Attribute(p0,nc,p1);
      if(!(RN_NC_IS(nc,RN_NC_ANY_NAME)&&p1==rn_text)) return 0;
      res|=2; break;
    case RN_P_TEXT: break;
    default: return 0;
    }
  }
  return res==3;
}

static int isanymix(int p) {
  int p0,nc,p1,p2,i,res,flat[2];
  p0=p; if(!RN_P_IS(p0,RN_P_ONE_OR_MORE)) return 0;
  rn_OneOrMore(p0,p1);
  p0=p1; if(!RN_P_IS(p0,RN_P_CHOICE)) return 0;
  rn_Choice(p0,p1,p2); flat[0]=p1; flat[1]=p2;
  res=0;
  for(i=0;i!=2;++i) {
    p0=flat[i];
    switch(RN_P_TYP(p0)) {
    case RN_P_ELEMENT: rn_Element(p0,nc,p1);
      if(!(RN_NC_IS(nc,RN_NC_ANY_NAME)&& isanycont(p1))) return 0;
      res|=1; break;
    case RN_P_TEXT: break;
    default: return 0;
    }
  }
  return res==1;
}

int ary_isany(int p) {
  int p1,p2;
  if(!RN_P_IS(p,RN_P_AFTER)) return 0;
  rn_After(p,p1,p2); return isanymix(p1)&&(p2==rn_empty||ary_isany(p2));
}
