/* $Id: rn.h,v 1.35 2004/02/25 00:00:32 dvd Exp $ */

#ifndef RN_H
#define RN_H 1

#include <assert.h>

/* Patterns */
#define RN_P_ERROR 0
#define RN_P_NOT_ALLOWED 1
#define RN_P_EMPTY 2
#define RN_P_TEXT 3
#define RN_P_CHOICE 4
#define RN_P_INTERLEAVE 5
#define RN_P_GROUP 6
#define RN_P_ONE_OR_MORE 7
#define RN_P_LIST 8
#define RN_P_DATA 9
#define RN_P_DATA_EXCEPT 10
#define RN_P_VALUE 11
#define RN_P_ATTRIBUTE 12
#define RN_P_ELEMENT 13
#define RN_P_REF 14
#define RN_P_AFTER 15

/*
Patterns and nameclasses are stored in arrays of integers.
an integer is either an index in the same or another array,
or a value that denotes record type etc.

Each record has a macro that accesses its fields by assigning
them to variables in the local scope, and a creator.
*/

/* Pattern Bindings */
#define RN_P_TYP(i) (rn_pattern[i]&0xFF)
#define RN_P_IS(i,x)  (x==RN_P_TYP(i))
#define RN_P_CHK(i,x)  assert(RN_P_IS(i,x))

#define RN_P_FLG_NUL 0x00000100
#define RN_P_FLG_TXT 0x00000200
#define RN_P_FLG_CTE 0x00000400
#define RN_P_FLG_CTC 0x00000800
#define RN_P_FLG_CTS 0x00001000
#define RN_P_FLG_ERS 0x40000000
#define RN_P_FLG_MRK 0x80000000

#define rn_marked(i) (rn_pattern[i]&RN_P_FLG_MRK)
#define rn_mark(i) (rn_pattern[i]|=RN_P_FLG_MRK)
#define rn_unmark(i) (rn_pattern[i]&=~RN_P_FLG_MRK)

#define rn_nullable(i) (rn_pattern[i]&RN_P_FLG_NUL)
#define rn_setNullable(i,x) if(x) rn_pattern[i]|=RN_P_FLG_NUL

#define rn_cdata(i) rn_pattern[i]&RN_P_FLG_TXT
#define rn_setCdata(i,x) if(x) rn_pattern[i]|=RN_P_FLG_TXT

/* assert: p1 at 1, p2 at 2 */

#define rn_NotAllowed(i) RN_P_CHK(i,RN_P_NOT_ALLOWED)
#define rn_Empty(i) RN_P_CHK(i,RN_P_EMPTY)
#define rn_Text(i) RN_P_CHK(i,RN_P_TEXT)
#define rn_Choice(i,p1,p2) RN_P_CHK(i,RN_P_CHOICE); p1=rn_pattern[i+1]; p2=rn_pattern[i+2]
#define rn_Interleave(i,p1,p2) RN_P_CHK(i,RN_P_INTERLEAVE); p1=rn_pattern[i+1]; p2=rn_pattern[i+2]
#define rn_Group(i,p1,p2) RN_P_CHK(i,RN_P_GROUP); p1=rn_pattern[i+1]; p2=rn_pattern[i+2]
#define rn_OneOrMore(i,p1) RN_P_CHK(i,RN_P_ONE_OR_MORE); p1=rn_pattern[i+1]
#define rn_List(i,p1) RN_P_CHK(i,RN_P_LIST); p1=rn_pattern[i+1]
#define rn_Data(i,dt,ps) RN_P_CHK(i,RN_P_DATA); dt=rn_pattern[i+1]; ps=rn_pattern[i+2]
#define rn_DataExcept(i,p1,p2) RN_P_CHK(i,RN_P_DATA_EXCEPT); p1=rn_pattern[i+1]; p2=rn_pattern[i+2]
#define rn_Value(i,dt,s) RN_P_CHK(i,RN_P_VALUE); dt=rn_pattern[i+1]; s=rn_pattern[i+2]
#define rn_Attribute(i,nc,p1) RN_P_CHK(i,RN_P_ATTRIBUTE);  p1=rn_pattern[i+1]; nc=rn_pattern[i+2]
#define rn_Element(i,nc,p1) RN_P_CHK(i,RN_P_ELEMENT); p1=rn_pattern[i+1]; nc=rn_pattern[i+2]
#define rn_After(i,p1,p2) RN_P_CHK(i,RN_P_AFTER); p1=rn_pattern[i+1]; p2=rn_pattern[i+2]
#define rn_Ref(i,p) RN_P_CHK(i,RN_P_REF); p=rn_pattern[i+1]

/* Name Classes */
#define RN_NC_ERROR 0
#define RN_NC_QNAME 1
#define RN_NC_NSNAME 2
#define RN_NC_ANY_NAME 3
#define RN_NC_EXCEPT 4
#define RN_NC_CHOICE 5
#define RN_NC_DATATYPE 6

/* Name Class Bindings  */
#define RN_NC_TYP(i) (rn_nameclass[i]&0xFF)
#define RN_NC_IS(i,x) (x==RN_NC_TYP(i))
#define RN_NC_CHK(i,x) assert(RN_NC_IS(i,x))

#define rn_QName(i,uri,name) RN_NC_CHK(i,RN_NC_QNAME); uri=rn_nameclass[i+1]; name=rn_nameclass[i+2]
#define rn_NsName(i,uri) RN_NC_CHK(i,RN_NC_NSNAME); uri=rn_nameclass[i+1]
#define rn_AnyName(i) RN_NC_CHK(i,RN_NC_ANY_NAME)
#define rn_NameClassExcept(i,nc1,nc2) RN_NC_CHK(i,RN_NC_EXCEPT); nc1=rn_nameclass[i+1]; nc2=rn_nameclass[i+2]
#define rn_NameClassChoice(i,nc1,nc2) RN_NC_CHK(i,RN_NC_CHOICE); nc1=rn_nameclass[i+1]; nc2=rn_nameclass[i+2]
#define rn_Datatype(i,lib,typ) RN_NC_CHK(i,RN_NC_DATATYPE); lib=rn_nameclass[i+1]; typ=rn_nameclass[i+2]

extern int rn_empty,rn_text,rn_notAllowed,rn_dt_string,rn_dt_token,rn_xsd_uri;

extern char *rn_string;

extern int *rn_pattern;
extern int *rn_nameclass;

extern void rn_new_schema(void);

extern int rn_contentType(int i);
extern void rn_setContentType(int i,int t1,int t2);
extern int rn_groupable(int p1,int p2);

extern void rn_del_p(int i);
extern void rn_add_p(int i);

extern int rn_newString(char *s);

extern int rn_newNotAllowed(void);
extern int rn_newEmpty(void);
extern int rn_newText(void);
extern int rn_newChoice(int p1,int p2);
extern int rn_newInterleave(int p1,int p2);
extern int rn_newGroup(int p1,int p2);
extern int rn_newOneOrMore(int p1);
extern int rn_newList(int p1);
extern int rn_newData(int dt,int ps);
extern int rn_newDataExcept(int p1,int p2);
extern int rn_newValue(int dt,int s);
extern int rn_newAttribute(int nc,int p1);
extern int rn_newElement(int nc,int p1);
extern int rn_newAfter(int p1,int p2);
extern int rn_newRef(void);

extern int rn_one_or_more(int p);
extern int rn_group(int p1,int p2);
extern int rn_choice(int p1,int p2);
extern int rn_ileave(int p1,int p2);
extern int rn_after(int p1,int p2);

extern int rn_newAnyName(void);
extern int rn_newAnyNameExcept(int nc);
extern int rn_newQName(int uri,int name);
extern int rn_newNsName(int uri);
extern int rn_newNameClassExcept(int nc1,int nc2);
extern int rn_newNameClassChoice(int nc1,int nc2);
extern int rn_newDatatype(int lib,int typ);

extern int rn_i_ps(void);
extern void rn_add_pskey(char *s);
extern void rn_add_psval(char *s);
extern void rn_end_ps(void);

extern void rn_init(void);
extern void rn_clear(void);

extern void rn_compress(int *starts,int n);
extern int rn_compress_last(int start);

#endif
