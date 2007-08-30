/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <limits.h>
#include "arraydef.h"
#include "seqpos-def.h"
#include "sarr-def.h"
#include "measure-time-if.h"

#include "sfx-map.pr"
#include "sfx-apfxlen.pr"
#include "sfx-suffixer.pr"

#define ISLEFTDIVERSE   (Uchar) (state->alphabetsize)
#define INITIALCHAR     (Uchar) (state->alphabetsize+1)

#define CHECKCHAR(CC)\
        if (father->commonchar != (CC) || (CC) >= ISLEFTDIVERSE)\
        {\
          father->commonchar = ISLEFTDIVERSE;\
        }

#define NODEPOSLISTENTRY(NN,SYM)\
        (NN)->nodeposlist[SYM]

#define NODEPOSLISTLENGTH(NN,SYM)\
        NODEPOSLISTENTRY(NN,SYM).length

#define NODEPOSLISTSTART(NN,SYM)\
        NODEPOSLISTENTRY(NN,SYM).start

typedef struct
{
  unsigned long start,
                length;
} Listtype;

typedef struct
{
  Uchar commonchar;
  unsigned long uniquecharposstart,
                uniquecharposlength; /* uniquecharpos[start..start+len-1] */
  Listtype *nodeposlist;
} Dfsinfo;

DECLAREARRAYSTRUCT(Seqpos);

typedef struct
{
  bool initialized;
  unsigned int searchlength,
               alphabetsize;
  Seqpos depth;            /* value changes with each new match */
  int(*output)(void *,Seqpos,Seqpos,Seqpos);
  void *outinfo;
  ArraySeqpos uniquechar,
              poslist[UCHAR_MAX+1];
} Dfsstate;

#include "esa-dfs.pr"

static Dfsinfo *allocateDfsinfo(Dfsstate *state,Env *env)
{
  Dfsinfo *dfsinfo;

  ALLOCASSIGNSPACE(dfsinfo,NULL,Dfsinfo,1);
  ALLOCASSIGNSPACE(dfsinfo->nodeposlist,NULL,Listtype,state->alphabetsize);
  return dfsinfo;
}

static void freeDfsinfo(Dfsinfo *dfsinfo,/*@unused@*/ Dfsstate *state,Env *env)
{
  FREESPACE(dfsinfo->nodeposlist);
  FREESPACE(dfsinfo);
}

static void add2poslist(Dfsstate *state,Dfsinfo *ninfo,unsigned int base,
                        Seqpos leafnumber,Env *env)
{
  ArraySeqpos *ptr;

  if (base >= state->alphabetsize)
  {
    ninfo->uniquecharposlength++;
    STOREINARRAY(&state->uniquechar,Seqpos,4,leafnumber);
  } else
  {
    ptr = &state->poslist[base];
    STOREINARRAY(ptr,Seqpos,4,leafnumber);
    NODEPOSLISTLENGTH(ninfo,base)++;
  }
}

static void concatlists(Dfsstate *state,Dfsinfo *father,Dfsinfo *son)
{
  unsigned int base;

  for (base = 0; base < state->alphabetsize; base++)
  {
    NODEPOSLISTLENGTH(father,base) += NODEPOSLISTLENGTH(son,base);
  }
  father->uniquecharposlength += son->uniquecharposlength;
}

static int cartproduct1(Dfsstate *state,const Dfsinfo *ninfo,unsigned int base,
                        Seqpos leafnumber)
{
  Listtype *pl;
  Seqpos *spptr, *start;

  pl = &NODEPOSLISTENTRY(ninfo,base);
  start = state->poslist[base].spaceSeqpos + pl->start;
  for (spptr = start; spptr < start + pl->length; spptr++)
  {
    if (state->output(state->outinfo,state->depth,leafnumber,*spptr) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static int cartproduct2(Dfsstate *state,
                        const Dfsinfo *ninfo1, unsigned int base1,
                        const Dfsinfo *ninfo2, unsigned int base2)
{
  Listtype *pl1, *pl2;
  Seqpos *start1, *start2, *spptr1, *spptr2;

  pl1 = &NODEPOSLISTENTRY(ninfo1,base1);
  start1 = state->poslist[base1].spaceSeqpos + pl1->start;
  pl2 = &NODEPOSLISTENTRY(ninfo2,base2);
  start2 = state->poslist[base2].spaceSeqpos + pl2->start;
  for (spptr1 = start1; spptr1 < start1 + pl1->length; spptr1++)
  {
    for (spptr2 = start2; spptr2 < start2 + pl2->length; spptr2++)
    {
      if (state->output(state->outinfo,state->depth,*spptr1,*spptr2) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

static void setpostabto0(Dfsstate *state)
{
  unsigned int base;

  if (!state->initialized)
  {
    for (base = 0; base < state->alphabetsize; base++)
    {
      state->poslist[base].nextfreeSeqpos = 0;
    }
    state->uniquechar.nextfreeSeqpos = 0;
    state->initialized = true;
  }
}

static int processleafedge(bool firstsucc,
                           Seqpos fatherdepth,
                           Dfsinfo *father,
                           Uchar leftchar,
                           Seqpos leafnumber,
                           Dfsstate *state,
                           Env *env)
{
  unsigned int base;
  Seqpos *start, *spptr;

#ifdef DEBUG
  printf("processleafedge " FormatSeqpos " firstsucc=%s, "
         " depth(father)= " FormatSeqpos "\n",
         PRINTSeqposcast(leafnumber),
         firstsucc ? "true" : "false",
         PRINTSeqposcast(fatherdepth));
#endif
  if (fatherdepth < (Seqpos) state->searchlength)
  {
    setpostabto0(state);
    return 0;
  }
  state->initialized = false;
  state->depth = fatherdepth;
#ifdef DEBUG
  printf("processleafedge: leftchar %u\n",(unsigned int) leftchar);
#endif
  if (firstsucc)
  {
    father->commonchar = leftchar;
    father->uniquecharposlength = 0;
    father->uniquecharposstart = state->uniquechar.nextfreeSeqpos;
    for (base = 0; base < state->alphabetsize; base++)
    {
      NODEPOSLISTSTART(father,base) = state->poslist[base].nextfreeSeqpos;
      NODEPOSLISTLENGTH(father,base) = 0;
    }
    add2poslist(state,father,(unsigned int) leftchar,leafnumber,env);
    return 0;
  }
  if (father->commonchar != ISLEFTDIVERSE)
  {
    CHECKCHAR(leftchar);
  }
  if (father->commonchar == ISLEFTDIVERSE)
  {
    for (base = 0; base < state->alphabetsize; base++)
    {
      if (leftchar != (Uchar) base)
      {
        if (cartproduct1(state,father,base,leafnumber) != 0)
        {
          return -1;
        }
      }
    }
    start = state->uniquechar.spaceSeqpos +
            father->uniquecharposstart;
    for (spptr = start; spptr < start + father->uniquecharposlength; spptr++)
    {
      if (state->output(state->outinfo,state->depth,leafnumber,*spptr) != 0)
      {
        return -2;
      }
    }
  }
  add2poslist(state,father,(unsigned int) leftchar,leafnumber,env);
  return 0;
}

static int processbranchedge(bool firstsucc,
                             Seqpos fatherdepth,
                             Dfsinfo *father,
                             Dfsinfo *son,
                             Dfsstate *state,
                             /*@unused@*/ Env *env)
{
  unsigned int chfather, chson;
  Seqpos *start, *spptr, *fptr, *fstart;

#ifdef DEBUG
  printf("processbranchedge firstsucc=%s, "
         " depth(father)= " FormatSeqpos "\n",
         firstsucc ? "true" : "false",
         PRINTSeqposcast(fatherdepth));
#endif
  if (fatherdepth < (Seqpos) state->searchlength)
  {
    setpostabto0(state);
    return 0;
  }
  state->initialized = false;
  state->depth = fatherdepth;
  if (firstsucc)
  {
    return 0;
  }
  if (father->commonchar != ISLEFTDIVERSE)
  {
    assert(son != NULL);
#ifdef DEBUG
    printf("commonchar=%u\n",(unsigned int) son->commonchar);
#endif
    if (son->commonchar != ISLEFTDIVERSE)
    {
      CHECKCHAR(son->commonchar);
    } else
    {
      father->commonchar = ISLEFTDIVERSE;
    }
  }
  if (father->commonchar == ISLEFTDIVERSE)
  {
    start = state->uniquechar.spaceSeqpos + son->uniquecharposstart;
    for (chfather = 0; chfather < state->alphabetsize; chfather++)
    {
      for (chson = 0; chson < state->alphabetsize; chson++)
      {
        if (chson != chfather)
        {
          if (cartproduct2(state,father,chfather,son,chson) != 0)
          {
            return -1;
          }
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (cartproduct1(state,father,chfather,*spptr) != 0)
        {
          return -2;
        }
      }
    }
    fstart = state->uniquechar.spaceSeqpos +
             father->uniquecharposstart;
    for (fptr = fstart; fptr < fstart + father->uniquecharposlength; fptr++)
    {
      for (chson = 0; chson < state->alphabetsize; chson++)
      {
        if (cartproduct1(state,son,chson,*fptr) != 0)
        {
          return -3;
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (state->output(state->outinfo,state->depth,*fptr,*spptr) != 0)
        {
          return -4;
        }
      }
    }
  }
  concatlists(state,father,son);
  return 0;
}

int enumeratemaxpairs(Sequentialsuffixarrayreader *ssar,
                      unsigned int searchlength,
                      int(*output)(void *,Seqpos,Seqpos,Seqpos),
                      void *outinfo,
                      Env *env)
{
  unsigned int base;
  ArraySeqpos *ptr;
  Dfsstate state;
  bool haserr = false;

  state.alphabetsize = alphabetsizeSequentialsuffixarrayreader(ssar);
  state.searchlength = searchlength;
  state.output = output;
  state.outinfo = outinfo;
  state.initialized = false;

  INITARRAY(&state.uniquechar,Seqpos);
  for (base = 0; base < state.alphabetsize; base++)
  {
    ptr = &state.poslist[base];
    INITARRAY(ptr,Seqpos);
  }
  if (depthfirstesa(ssar,
                    (Uchar) (state.alphabetsize+1),
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    processbranchedge,
                    NULL,
                    NULL,
                    NULL,
                    &state,
                    env) != 0)
  {
    haserr = true;
  }
  FREEARRAY(&state.uniquechar,Seqpos);
  for (base = 0; base < state.alphabetsize; base++)
  {
    ptr = &state.poslist[base];
    FREEARRAY(ptr,Seqpos);
  }
  return haserr ? -1 : 0;
}

typedef struct
{
  unsigned int minlength;
  Encodedsequence *encseq;
  int (*processmaxmatch)(void *,Seqpos,Seqpos,Seqpos);
  void *processmaxmatchinfo;
} Substringmatchinfo;

static int processsuftab(/*@unused@*/ void *info,
                         /*@unused@*/ const Seqpos *suftabpart,
                         /*@unused@*/ Readmode readmode,
                         /*@unused@*/ Seqpos widthofpart,
                         /*@unused@*/ Env *env)
{
  return 0;
}

int sarrselfsubstringmatch(const Uchar *dbseq,
                           Seqpos dblen,
                           const Uchar *query,
                           unsigned long querylen,
                           unsigned int minlength,
                           const Alphabet *alpha,
                           int (*processmaxmatch)(void *,Seqpos,
                                                  Seqpos,Seqpos),
                           void *processmaxmatchinfo,
                           Env *env)
{
  Specialcharinfo samplespecialcharinfo;
  Substringmatchinfo ssi;
  unsigned int numofchars;
  bool haserr = false;

  ssi.encseq = plain2encodedsequence(true,
                                     &samplespecialcharinfo,
                                     dbseq,
                                     dblen,
                                     query,
                                     querylen,
                                     alpha,
                                     env);
  ssi.minlength = minlength;
  ssi.processmaxmatch = processmaxmatch;
  ssi.processmaxmatchinfo = processmaxmatchinfo;
  numofchars = getnumofcharsAlphabet(alpha);
  if (suffixerator(processsuftab,
                   &ssi,
                   samplespecialcharinfo.specialcharacters,
                   samplespecialcharinfo.specialranges,
                   ssi.encseq,
                   Forwardmode,
                   numofchars,
                   recommendedprefixlength(numofchars,dblen),
                   (unsigned int) 1, /* parts */
                   NULL,
		   env) != 0)
  {
    haserr = true;
  }
  freeEncodedsequence(&ssi.encseq,env);
  return haserr ? -1 : 0;
}
