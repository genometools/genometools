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

#include "core/arraydef.h"
#include "core/unused_api.h"
#include "esa-seqread.h"
#include "esa-maxpairs.h"

#define ISLEFTDIVERSE   (GtUchar) (state->alphabetsize)
#define INITIALCHAR     (GtUchar) (state->alphabetsize+1)

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

typedef struct /* information stored for each node of the lcp interval tree */
{
  GtUchar commonchar;
  unsigned long uniquecharposstart,
                uniquecharposlength; /* uniquecharpos[start..start+len-1] */
  Listtype *nodeposlist;
} GtBUinfo_maxpairs;

typedef struct  /* global information */
{
  bool initialized;
  unsigned int searchlength,
               alphabetsize;
  GtArrayGtUlong uniquechar,
              *poslist;
  const GtEncseq *encseq;
  GtReadmode readmode;
  Processmaxpairs processmaxpairs;
  void *processmaxpairsinfo;
} GtBUstate_maxpairs;

static void initBUinfo_maxpairs(GtBUinfo_maxpairs *buinfo,
                                GtBUstate_maxpairs *state)
{
  buinfo->nodeposlist = gt_malloc(sizeof (*buinfo->nodeposlist) *
                                  state->alphabetsize);
}

static void freeBUinfo_maxpairs(GtBUinfo_maxpairs *buinfo,
                                GT_UNUSED GtBUstate_maxpairs *state)
{
  gt_free(buinfo->nodeposlist);
}

static void add2poslist_maxpairs(GtBUstate_maxpairs *state,
                                 GtBUinfo_maxpairs *ninfo,
                                 unsigned int base,unsigned long leafnumber)
{
  GtArrayGtUlong *ptr;

  if (base >= state->alphabetsize)
  {
    ninfo->uniquecharposlength++;
    GT_STOREINARRAY(&state->uniquechar,GtUlong,4,leafnumber);
  } else
  {
    ptr = &state->poslist[base];
    GT_STOREINARRAY(ptr,GtUlong,4,leafnumber);
    NODEPOSLISTLENGTH(ninfo,base)++;
  }
}

static void concatlists_maxpairs(GtBUstate_maxpairs *state,
                                 GtBUinfo_maxpairs *father,
                                 GtBUinfo_maxpairs *son)
{
  unsigned int base;

  for (base = 0; base < state->alphabetsize; base++)
  {
    NODEPOSLISTLENGTH(father,base) += NODEPOSLISTLENGTH(son,base);
  }
  father->uniquecharposlength += son->uniquecharposlength;
}

static int cartproduct1_maxpairs(GtBUstate_maxpairs *state,
                                 unsigned long fatherdepth,
                                 const GtBUinfo_maxpairs *ninfo,
                                 unsigned int base,
                                 unsigned long leafnumber,GtError *err)
{
  Listtype *pl;
  unsigned long *spptr, *start;

  pl = &NODEPOSLISTENTRY(ninfo,base);
  start = state->poslist[base].spaceGtUlong + pl->start;
  for (spptr = start; spptr < start + pl->length; spptr++)
  {
    if (state->processmaxpairs(state->processmaxpairsinfo,state->encseq,
                               fatherdepth,leafnumber,*spptr,err) != 0)
    {
      return -1;
    }
  }
  return 0;
}

static int cartproduct2_maxpairs(GtBUstate_maxpairs *state,
                                 unsigned long fatherdepth,
                                 const GtBUinfo_maxpairs *ninfo1,
                                 unsigned int base1,
                                 const GtBUinfo_maxpairs *ninfo2,
                                 unsigned int base2,
                                 GtError *err)
{
  Listtype *pl1, *pl2;
  unsigned long *start1, *start2, *spptr1, *spptr2;

  pl1 = &NODEPOSLISTENTRY(ninfo1,base1);
  start1 = state->poslist[base1].spaceGtUlong + pl1->start;
  pl2 = &NODEPOSLISTENTRY(ninfo2,base2);
  start2 = state->poslist[base2].spaceGtUlong + pl2->start;
  for (spptr1 = start1; spptr1 < start1 + pl1->length; spptr1++)
  {
    for (spptr2 = start2; spptr2 < start2 + pl2->length; spptr2++)
    {
      if (state->processmaxpairs(state->processmaxpairsinfo,state->encseq,
                                 fatherdepth,*spptr1,*spptr2,err) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

static void setpostabto0_maxpairs(GtBUstate_maxpairs *state)
{
  unsigned int base;

  if (!state->initialized)
  {
    for (base = 0; base < state->alphabetsize; base++)
    {
      state->poslist[base].nextfreeGtUlong = 0;
    }
    state->uniquechar.nextfreeGtUlong = 0;
    state->initialized = true;
  }
}

static int processleafedge_maxpairs(bool firstsucc,
                                    unsigned long fatherdepth,
                                    GtBUinfo_maxpairs *father,
                                    unsigned long leafnumber,
                                    GtBUstate_maxpairs *state,
                                    GtError *err)
{
  unsigned int base;
  unsigned long *start, *spptr;
  GtUchar leftchar;

#ifdef SKDEBUG
  printf("%s %lu firstsucc=%s, __func__," " depth(father)= %lu\n",
         leafnumber,
         firstsucc ? "true" : "false",
         fatherdepth);
#endif
  if (fatherdepth < (unsigned long) state->searchlength)
  {
    setpostabto0_maxpairs(state);
    return 0;
  }
  if (leafnumber == 0)
  {
    leftchar = INITIALCHAR;
  } else
  {
    /* Random access */
    leftchar = gt_encseq_get_encoded_char(state->encseq,
                                          leafnumber-1,
                                          state->readmode);
  }
  state->initialized = false;
#ifdef SKDEBUG
  printf("%s: leftchar %u\n",__func__,(unsigned int) leftchar);
#endif
  if (firstsucc)
  {
    father->commonchar = leftchar;
    father->uniquecharposlength = 0;
    father->uniquecharposstart = state->uniquechar.nextfreeGtUlong;
    for (base = 0; base < state->alphabetsize; base++)
    {
      NODEPOSLISTSTART(father,base) = state->poslist[base].nextfreeGtUlong;
      NODEPOSLISTLENGTH(father,base) = 0;
    }
    add2poslist_maxpairs(state,father,(unsigned int) leftchar,leafnumber);
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
      if (leftchar != (GtUchar) base)
      {
        if (cartproduct1_maxpairs(state,fatherdepth,father,base,leafnumber,
                                  err) != 0)
        {
          return -1;
        }
      }
    }
    start = state->uniquechar.spaceGtUlong +
            father->uniquecharposstart;
    for (spptr = start; spptr < start + father->uniquecharposlength; spptr++)
    {
      if (state->processmaxpairs(state->processmaxpairsinfo,state->encseq,
                                 fatherdepth,leafnumber,*spptr,err) != 0)
      {
        return -2;
      }
    }
  }
  add2poslist_maxpairs(state,father,(unsigned int) leftchar,leafnumber);
  return 0;
}

static int processbranchingedge_maxpairs(bool firstsucc,
                                         unsigned long fatherdepth,
                                         GtBUinfo_maxpairs *father,
                                         GT_UNUSED unsigned long sondepth,
                                         GT_UNUSED unsigned long sonwidth,
                                         GtBUinfo_maxpairs *son,
                                         GtBUstate_maxpairs *state,
                                         GtError *err)
{
  unsigned int chfather, chson;
  unsigned long *start, *spptr, *fptr, *fstart;

#ifdef SKDEBUG
  printf("%s firstsucc=%s, depth(father)= %lu\n",
          __func__,firstsucc ? "true" : "false",fatherdepth);
#endif
  if (fatherdepth < (unsigned long) state->searchlength)
  {
    setpostabto0_maxpairs(state);
    return 0;
  }
  state->initialized = false;
  if (firstsucc)
  {
    return 0;
  }
  if (father->commonchar != ISLEFTDIVERSE)
  {
    gt_assert(son != NULL);
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
    start = state->uniquechar.spaceGtUlong + son->uniquecharposstart;
    for (chfather = 0; chfather < state->alphabetsize; chfather++)
    {
      for (chson = 0; chson < state->alphabetsize; chson++)
      {
        if (chson != chfather)
        {
          if (cartproduct2_maxpairs(state,fatherdepth,father,chfather,
                                    son,chson,err) != 0)
          {
            return -1;
          }
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (cartproduct1_maxpairs(state,fatherdepth,father,chfather,*spptr,
                                  err) != 0)
        {
          return -2;
        }
      }
    }
    fstart = state->uniquechar.spaceGtUlong +
             father->uniquecharposstart;
    for (fptr = fstart; fptr < fstart + father->uniquecharposlength; fptr++)
    {
      for (chson = 0; chson < state->alphabetsize; chson++)
      {
        if (cartproduct1_maxpairs(state,fatherdepth,son,chson,*fptr,err) != 0)
        {
          return -3;
        }
      }
      for (spptr = start; spptr < start + son->uniquecharposlength; spptr++)
      {
        if (state->processmaxpairs(state->processmaxpairsinfo,state->encseq,
                                   fatherdepth,*fptr,*spptr,err) != 0)
        {
          return -4;
        }
      }
    }
  }
  concatlists_maxpairs(state,father,son);
  return 0;
}

#include "esa-bottomup-maxpairs.inc"

int gt_enumeratemaxpairs(Sequentialsuffixarrayreader *ssar,
                         const GtEncseq *encseq,
                         GtReadmode readmode,
                         unsigned int searchlength,
                         Processmaxpairs processmaxpairs,
                         void *processmaxpairsinfo,
                         GtError *err)
{
  unsigned int base;
  GtArrayGtUlong *ptr;
  GtBUstate_maxpairs *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  state->alphabetsize = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  state->searchlength = searchlength;
  state->processmaxpairs = processmaxpairs;
  state->processmaxpairsinfo = processmaxpairsinfo;
  state->initialized = false;
  state->encseq = encseq;
  state->readmode = readmode;

  GT_INITARRAY(&state->uniquechar,GtUlong);
  state->poslist = gt_malloc(sizeof (*state->poslist) * state->alphabetsize);
  for (base = 0; base < state->alphabetsize; base++)
  {
    ptr = &state->poslist[base];
    GT_INITARRAY(ptr,GtUlong);
  }
  if (gt_esa_bottomup_maxpairs(ssar, state, err) != 0)
  {
    haserr = true;
  }
  GT_FREEARRAY(&state->uniquechar,GtUlong);
  for (base = 0; base < state->alphabetsize; base++)
  {
    ptr = &state->poslist[base];
    GT_FREEARRAY(ptr,GtUlong);
  }
  gt_free(state->poslist);
  gt_free(state);
  return haserr ? -1 : 0;
}
