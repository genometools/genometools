/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "libgtcore/arraydef.h"
#include "libgtcore/unused.h"
#include "spacedef.h"
#include "seqpos-def.h"
#include "esa-seqread.h"
#include "verbose-def.h"

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
  ArraySeqpos uniquechar,
              *poslist;
  const Encodedsequence *encseq;
  Readmode readmode;
  int(*processmaxpairs)(void *,Seqpos,Seqpos,Seqpos,Error *);
  void *processmaxpairsinfo;
} Dfsstate;

#include "esa-dfs.pr"

static Dfsinfo *allocateDfsinfo(Dfsstate *state)
{
  Dfsinfo *dfsinfo;

  ALLOCASSIGNSPACE(dfsinfo,NULL,Dfsinfo,1);
  ALLOCASSIGNSPACE(dfsinfo->nodeposlist,NULL,Listtype,state->alphabetsize);
  return dfsinfo;
}

static void freeDfsinfo(Dfsinfo *dfsinfo, UNUSED Dfsstate *state)
{
  FREESPACE(dfsinfo->nodeposlist);
  FREESPACE(dfsinfo);
}

static void add2poslist(Dfsstate *state,Dfsinfo *ninfo,unsigned int base,
                        Seqpos leafnumber)
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
                        Seqpos leafnumber,Error *err)
{
  Listtype *pl;
  Seqpos *spptr, *start;

  pl = &NODEPOSLISTENTRY(ninfo,base);
  start = state->poslist[base].spaceSeqpos + pl->start;
  for (spptr = start; spptr < start + pl->length; spptr++)
  {
    if (state->processmaxpairs(state->processmaxpairsinfo,
                               state->depth,leafnumber,*spptr,err) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static int cartproduct2(Dfsstate *state,
                        const Dfsinfo *ninfo1, unsigned int base1,
                        const Dfsinfo *ninfo2, unsigned int base2,
                        Error *err)
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
      if (state->processmaxpairs(state->processmaxpairsinfo,
                                 state->depth,*spptr1,*spptr2,err) != 0)
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
                           Seqpos leafnumber,
                           Dfsstate *state,
                           Error *err)
{
  unsigned int base;
  Seqpos *start, *spptr;
  Uchar leftchar;

#ifdef SKDEBUG
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
  if (leafnumber == 0)
  {
    leftchar = INITIALCHAR;
  } else
  {
    leftchar = getencodedchar(state->encseq, /* Random access */
                              leafnumber-1,
                              state->readmode);
  }
  state->initialized = false;
  state->depth = fatherdepth;
#ifdef SKDEBUG
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
    add2poslist(state,father,(unsigned int) leftchar,leafnumber);
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
        if (cartproduct1(state,father,base,leafnumber,err) != 0)
        {
          return -1;
        }
      }
    }
    start = state->uniquechar.spaceSeqpos +
            father->uniquecharposstart;
    for (spptr = start; spptr < start + father->uniquecharposlength; spptr++)
    {
      if (state->processmaxpairs(state->processmaxpairsinfo,
                                 state->depth,leafnumber,*spptr,err) != 0)
      {
        return -2;
      }
    }
  }
  add2poslist(state,father,(unsigned int) leftchar,leafnumber);
  return 0;
}

static int processbranchedge(bool firstsucc,
                             Seqpos fatherdepth,
                             Dfsinfo *father,
                             Dfsinfo *son,
                             Dfsstate *state,
                             /*@unused@*/ Error *err)
{
  unsigned int chfather, chson;
  Seqpos *start, *spptr, *fptr, *fstart;

#ifdef SKDEBUG
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
#ifdef SKDEBUG
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
          if (cartproduct2(state,father,chfather,son,chson,err) != 0)
          {
            return -1;
          }
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (cartproduct1(state,father,chfather,*spptr,err) != 0)
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
        if (cartproduct1(state,son,chson,*fptr,err) != 0)
        {
          return -3;
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (state->processmaxpairs(state->processmaxpairsinfo,
                                   state->depth,*fptr,*spptr,err) != 0)
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
                      unsigned int alphabetsize,
                      const Encodedsequence *encseq,
                      Readmode readmode,
                      unsigned int searchlength,
                      int(*processmaxpairs)(void *,Seqpos,Seqpos,
                                            Seqpos,Error *),
                      void *processmaxpairsinfo,
                      Verboseinfo *verboseinfo,
                      Error *err)
{
  unsigned int base;
  ArraySeqpos *ptr;
  Dfsstate state;
  bool haserr = false;

  state.alphabetsize = alphabetsize;
  state.searchlength = searchlength;
  state.processmaxpairs = processmaxpairs;
  state.processmaxpairsinfo = processmaxpairsinfo;
  state.initialized = false;
  state.encseq = encseq;
  state.readmode = readmode;

  INITARRAY(&state.uniquechar,Seqpos);
  ALLOCASSIGNSPACE(state.poslist,NULL,ArraySeqpos,alphabetsize);
  for (base = 0; base < state.alphabetsize; base++)
  {
    ptr = &state.poslist[base];
    INITARRAY(ptr,Seqpos);
  }
  if (depthfirstesa(ssar,
                    allocateDfsinfo,
                    freeDfsinfo,
                    processleafedge,
                    processbranchedge,
                    /*
                    NULL,
                    NULL,
                    NULL,
                    */
                    &state,
                    verboseinfo,
                    err) != 0)
  {
    haserr = true;
  }
  FREEARRAY(&state.uniquechar,Seqpos);
  for (base = 0; base < state.alphabetsize; base++)
  {
    ptr = &state.poslist[base];
    FREEARRAY(ptr,Seqpos);
  }
  FREESPACE(state.poslist);
  return haserr ? -1 : 0;
}

int callenummaxpairs(const Str *indexname,
                     unsigned int userdefinedleastlength,
                     bool scanfile,
                     int(*processmaxpairs)(void *,Seqpos,Seqpos,
                                           Seqpos,Error *),
                     void *processmaxpairsinfo,
                     Verboseinfo *verboseinfo,
                     Error *err)
{
  bool haserr = false;
  Sequentialsuffixarrayreader *ssar;

  error_check(err);
  ssar = newSequentialsuffixarrayreaderfromfile(indexname,
                                                SARR_LCPTAB |
                                                SARR_SUFTAB |
                                                SARR_ESQTAB,
                                                scanfile
                                                  ? SEQ_scan : SEQ_mappedboth,
                                                err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr &&
      enumeratemaxpairs(ssar,
                        getnumofcharsAlphabet(
                                 alphabetSequentialsuffixarrayreader(ssar)),
                        encseqSequentialsuffixarrayreader(ssar),
                        readmodeSequentialsuffixarrayreader(ssar),
                        userdefinedleastlength,
                        processmaxpairs,
                        processmaxpairsinfo,
                        verboseinfo,
                        err) != 0)
  {
    haserr = true;
  }
  if (ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar);
  }
  return haserr ? -1 : 0;
}
