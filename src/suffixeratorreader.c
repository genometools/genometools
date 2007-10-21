/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/

#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <libgtcore/filelengthvalues.h>
#include <libgtcore/minmax.h>
#include "libgtcore/seqiterator.h"
#include <libgtcore/str.h>
#include <libgtcore/strarray.h>
#include <libgtmatch/seqpos-def.h>
#include <libgtcore/symboldef.h>
#include <libgtmatch/fillsci.pr>
#include <libgtmatch/opensfxfile.pr>
#include <libgtmatch/sfx-lcpval.h>
#include <libgtmatch/esafileend.h>
#include <libgtmatch/sfx-apfxlen.pr>
#include <libgtmatch/intcode-def.h>
#include <libgtmatch/sfx-outprj.pr>
#include <libgtmatch/sfx-cmpsuf.pr>

#include "suffixerator-interface.h"

typedef struct
{
  FILE *outfpsuftab,
       *outfplcptab,
       *outfpllvtab,
       *outfpbwttab;
  Seqpos pageoffset,
         numoflargelcpvalues,
         maxbranchdepth;
  Encodedsequence *encseq;
  DefinedSeqpos longest;
  Lcpvalueiterator *lvi;
} Outfileinfo;

struct sfxIReaderState
{
  Seqpos nextReadPos;
  int readFlag;
};

struct sfxInterface
{
  Outfileinfo outfileinfo;
  Filelengthvalues *filelengthtab;
  Sfxiterator *sfi;
  Alphabet *alpha;
  Seqpos length;
  unsigned long numofsequences;
  Specialcharinfo specialcharinfo;
  Suffixeratoroptions so;

  Seqpos lastGeneratedLen, lastGeneratedStart;
  const Seqpos *lastGeneratedSufTabSegment;
  /**< pointer to part of suffix array
   * returned by last iterator call*/
  Seqpos *prevGeneratedSufTabSegments; /**< holds cache of suffix array
                                  * parts returned by previous to last
                                  * iterator calls as far back as
                                  * required by the readers */
  Seqpos prevGeneratedStart;
  size_t prevGeneratedSize,     /**< maximum number of Seqpos values
                                 * the cache may hold  */
    prevGeneratedLen;           /**< number of Seqpos values the
                                 * cache holds currently */
  size_t numReaders;
  int allRequests;
  struct sfxIReaderState *readers;
  bool specialsuffixes;
};

static Seqpos
SfxIReadAdvance(sfxInterface *iface,
                Seqpos requestMaxPos,
                Env *env);


extern sfxInterface *
newSfxInterface(Suffixeratoroptions *so,
                Env *env)
{
  return newSfxInterfaceWithReaders(so, 0, NULL, NULL, env);
}

static void initOutfileInfo(Outfileinfo *outfileinfo)
{
  outfileinfo->outfpsuftab = NULL;
  outfileinfo->outfplcptab = NULL;
  outfileinfo->outfpllvtab = NULL;
  outfileinfo->outfpbwttab = NULL;
  outfileinfo->pageoffset = 0;
  outfileinfo->numoflargelcpvalues = 0;
  outfileinfo->maxbranchdepth = 0;
  outfileinfo->longest.defined = false;
  outfileinfo->longest.valueseqpos = 0;
}

static unsigned long *initcharacterdistribution(const Alphabet *alpha,Env *env)
{
  unsigned long *characterdistribution;
  unsigned int mapsize, idx;

  mapsize = getmapsizeAlphabet(alpha);
  characterdistribution = env_ma_malloc(env,
                                        sizeof(unsigned long) * (mapsize-1));
  for (idx=0; idx<mapsize-1; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
}

static int outputsequencedescription(const Str *indexname,
                                     const StrArray *filenametab,
                                     Env *env)
{
  FILE *desfp;
  bool haserr = false;

  desfp = opensfxfile(indexname,DESTABSUFFIX,"wb",env);
  if (desfp == NULL)
  {
    haserr = true;
  } else
  {
    SeqIterator *seqit;
    char *desc = NULL;
    unsigned long seqlen;
    int retval;

    seqit = seqiterator_new(filenametab,NULL,false,env);
    while (!haserr)
    {
      retval = seqiterator_next(seqit,
                                NULL,
                                &seqlen,
                                &desc,
                                env);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      if (fputs(desc,desfp) == EOF)
      {
        env_error_set(env,"cannot write description to file %s.%s",
                          str_get(indexname),DESTABSUFFIX);
        haserr = true;
      }
      (void) putc((int) '\n',desfp);
      env_ma_free(desc, env);
    }
    env_fa_xfclose(desfp,env);
    seqiterator_delete(seqit,env);
  }
  return haserr ? -1 : 0;
}

static void showcharacterdistribution(
                   const  Alphabet *alpha,
                   const unsigned long *characterdistribution)
{
  unsigned int mapsize, idx;

  mapsize = getmapsizeAlphabet(alpha);
  assert(characterdistribution != NULL);
  for (idx=0; idx<mapsize-1; idx++)
  {
    printf("# ocurrences(%c)=%lu\n",(int) getprettysymbol(alpha,idx),
                                    characterdistribution[idx]);
  }
}

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)                                 \
  ((FLAG)?(PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",env)):((void*)1))

#define newSfxInterfaceWithReadersErrRet()        \
  do {                                            \
    if(iface->alpha)                              \
      freeAlphabet(&iface->alpha, env);           \
    if(iface->filelengthtab)                      \
      env_ma_free(iface->filelengthtab, env);     \
    if(iface) env_ma_free(iface, env);            \
                                                  \
  } while(0)

extern sfxInterface *
newSfxInterfaceWithReaders(Suffixeratoroptions *so,
                           size_t numReaders,
                           enum sfxDataRequest *requestors,
                           listenerID *ids, Env *env)
{
  Encodedsequence *encseq = NULL;
  unsigned long *characterdistribution = NULL;
  unsigned int numofchars = 0;
  sfxInterface *iface = NULL;

  env_error_check(env);

  iface = env_ma_calloc(env, 1, sizeof(iface));

  memcpy(&iface->so, so, sizeof(*so));
  
  if(!(iface->alpha = assigninputalphabet(so->isdna, so->isprotein,
                                   so->str_smap, so->filenametab, env)))
    newSfxInterfaceWithReadersErrRet();

  if (!so->isplain)
    characterdistribution = initcharacterdistribution(iface->alpha,env);

  if (fasta2sequencekeyvalues(&iface->numofsequences, &iface->length,
                              &iface->specialcharinfo, so->filenametab,
                              &iface->filelengthtab,
                              getsymbolmapAlphabet(iface->alpha),
                              so->isplain, characterdistribution,
                              env) != 0)
    newSfxInterfaceWithReadersErrRet();
    
  numofchars = getnumofcharsAlphabet(iface->alpha);

  if(!(encseq = files2encodedsequence(
         true, so->filenametab, so->isplain,
         iface->length, &iface->specialcharinfo, iface->alpha,
         str_length(so->str_sat) > 0 ? str_get(so->str_sat) : NULL,
         env)))
    newSfxInterfaceWithReadersErrRet();

  if (so->outtistab
      && flushencseqfile(so->str_indexname, encseq, env) != 0)
    newSfxInterfaceWithReadersErrRet();

  if(so->outdestab
     && outputsequencedescription(so->str_indexname, so->filenametab, env))
    newSfxInterfaceWithReadersErrRet();

  initOutfileInfo(&iface->outfileinfo);
  if (so->outlcptab)
    iface->outfileinfo.lvi = newLcpvalueiterator(encseq,so->readmode,env);
  else
    iface->outfileinfo.lvi = NULL;

  iface->outfileinfo.encseq = encseq;

  if(!INITOUTFILEPTR(iface->outfileinfo.outfpsuftab,
                     so->outsuftab, SUFTABSUFFIX))
    newSfxInterfaceWithReadersErrRet();
  if(!INITOUTFILEPTR(iface->outfileinfo.outfplcptab,
                     so->outlcptab, LCPTABSUFFIX))
    newSfxInterfaceWithReadersErrRet();
  if(!INITOUTFILEPTR(iface->outfileinfo.outfpllvtab,
                     so->outlcptab, LARGELCPTABSUFFIX))
    newSfxInterfaceWithReadersErrRet();
  if(!INITOUTFILEPTR(iface->outfileinfo.outfpbwttab,
                     so->outbwttab, BWTTABSUFFIX))
    newSfxInterfaceWithReadersErrRet();

  printf("# specialcharacters=" FormatSeqpos "\n",
         PRINTSeqposcast(iface->specialcharinfo.specialcharacters));
  printf("# specialranges=" FormatSeqpos "\n",
         PRINTSeqposcast(iface->specialcharinfo.specialranges));
  if (!so->isplain)
  {
    showcharacterdistribution(iface->alpha,characterdistribution);
  }
  env_ma_free(characterdistribution, env);

  if((so->readmode == Complementmode || so->readmode == Reversecomplementmode)
     && !isdnaalphabet(iface->alpha,env))
  {
    env_error_set(env,"option %s only can be used for DNA alphabets",
                  so->readmode == Complementmode ? "-cpl" : "rcl");
    newSfxInterfaceWithReadersErrRet();
  }

  if (so->outsuftab || so->outbwttab || so->outlcptab)
  {
    if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
    {
      iface->so.prefixlength = so->prefixlength =
        recommendedprefixlength(numofchars, iface->length);
      printf("# automatically determined prefixlength = %u\n",
             so->prefixlength);
    }
    else
    {
      unsigned int maxprefixlen;
      
      maxprefixlen
        = whatisthemaximalprefixlength(numofchars,
                                       iface->length,
                                       (unsigned int) PREFIXLENBITS);
      if (checkprefixlength(maxprefixlen,so->prefixlength,env) != 0)
        newSfxInterfaceWithReadersErrRet();
      else
        showmaximalprefixlength(maxprefixlen,
                                recommendedprefixlength(
                                  numofchars,
                                  iface->length));
    }
  }
  else if (so->readmode != Forwardmode)
  {
    env_error_set(env,"option '-dir %s' only makes sense in combination "
                  "with at least one of the options -suf, -lcp, or "
                  "-bwt", showreadmode(so->readmode));
    newSfxInterfaceWithReadersErrRet();        
  }
  /* FIXME: replace with different code in iterator */
  if(!(iface->sfi = newSfxiterator(iface->specialcharinfo.specialcharacters,
                                   iface->specialcharinfo.specialranges,
                                   encseq, so->readmode,
                                   numofchars, so->prefixlength,
                                   so->numofparts, NULL, env)))
    newSfxInterfaceWithReadersErrRet();
  {
    size_t i;
    for(i = 0; i < numReaders; ++i)
      if(!SfxIRegisterReader(iface, ids + i, requestors[i], env))
        newSfxInterfaceWithReadersErrRet();
  }
  iface->specialsuffixes = false;
  iface->allRequests = SFX_REQUEST_NONE;
  iface->lastGeneratedSufTabSegment
    = nextSfxiterator(&iface->lastGeneratedLen, &iface->specialsuffixes,
                      NULL, iface->sfi, env);
  iface->lastGeneratedStart = 0;
  return iface;
}


extern int
deleteSfxInterface(sfxInterface *iface, Env *env)
{
  int had_err = 0;

  env_ma_free(iface->prevGeneratedSufTabSegments, env);
  
  freeSfxiterator(&iface->sfi, env);

  env_fa_fclose(iface->outfileinfo.outfpsuftab,env);
  env_fa_fclose(iface->outfileinfo.outfplcptab,env);
  env_fa_fclose(iface->outfileinfo.outfpllvtab,env);
  env_fa_fclose(iface->outfileinfo.outfpbwttab,env);

  if (iface->outfileinfo.lvi != NULL)
    freeLcpvalueiterator(&iface->outfileinfo.lvi,env);

  if (outprjfile(iface->so.str_indexname,
                 iface->so.filenametab,
                 iface->so.readmode,
                 iface->filelengthtab,
                 iface->length,
                 iface->numofsequences,
                 &iface->specialcharinfo,
                 iface->so.prefixlength,
                 iface->outfileinfo.numoflargelcpvalues,
                 iface->outfileinfo.maxbranchdepth,
                 &iface->outfileinfo.longest,
                 env) != 0)
    had_err = 1;
  env_ma_free(iface->filelengthtab, env);
  if(iface->alpha)
    freeAlphabet(&iface->alpha, env);

  freeEncodedsequence(&iface->outfileinfo.encseq, env);
  return !had_err;
}

static int outlcpvalue(Seqpos lcpvalue,Seqpos pos,Outfileinfo *outfileinfo,
                       Env *env)
{
  Uchar outvalue;
  bool haserr = false;

  if (lcpvalue >= (Seqpos) UCHAR_MAX)
  {
    Largelcpvalue largelcpvalue;

    outfileinfo->numoflargelcpvalues++;
    largelcpvalue.position = outfileinfo->pageoffset + pos;
    largelcpvalue.value = lcpvalue;
    if (fwrite(&largelcpvalue,sizeof (Largelcpvalue),(size_t) 1,
               outfileinfo->outfpllvtab) != (size_t) 1)
    {
      env_error_set(env,"cannot write 1 item of size %lu: "
                        "errormsg=\"%s\"",
                        (unsigned long) sizeof (Largelcpvalue),
                        strerror(errno));
      haserr = true;
    }
    outvalue = (Uchar) UCHAR_MAX;
  } else
  {
    outvalue = (Uchar) lcpvalue;
  }
  if (!haserr && fwrite(&outvalue,sizeof (Uchar),(size_t) 1,
                        outfileinfo->outfplcptab) != (size_t) 1)
  {
    env_error_set(env,"cannot write 1 item of size %lu: "
                      "errormsg=\"%s\"",
                      (unsigned long) sizeof (Uchar),
                      strerror(errno));
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static int suftab2file(Outfileinfo *outfileinfo,
                       const Seqpos *suftab,
                       Readmode readmode,
                       Seqpos numberofsuffixes,
                       Env *env)
{
  bool haserr = false;

  env_error_check(env);
  if (outfileinfo->outfpsuftab != NULL)
  {
    if (fwrite(suftab,
              sizeof (*suftab),
              (size_t) numberofsuffixes,
              outfileinfo->outfpsuftab)
              != (size_t) numberofsuffixes)
    {
      env_error_set(env,"cannot write " FormatSeqpos " items of size %u: "
                    "errormsg=\"%s\"",
           PRINTSeqposcast(numberofsuffixes),
           (unsigned int) sizeof (*suftab),
           strerror(errno));
      haserr = true;
    }
  }
  if (!haserr)
  {
    Seqpos startpos, lcpvalue, pos;
    Uchar cc = 0;

    for (pos=0; pos < numberofsuffixes; pos++)
    {
      startpos = suftab[pos];
      if (startpos == 0)
      {
        cc = (Uchar) UNDEFBWTCHAR;
        if (outfileinfo->longest.defined)
        {
          env_error_set(env,"longest = " FormatSeqpos " is already defined",
                        PRINTSeqposcast(outfileinfo->longest.valueseqpos));
          haserr = true;
          break;
        }
        outfileinfo->longest.defined = true;
        outfileinfo->longest.valueseqpos = outfileinfo->pageoffset + pos;
      } else
      {
        if (outfileinfo->outfpbwttab != NULL)
        {
          cc = getencodedchar(outfileinfo->encseq,startpos - 1,readmode);
        }
      }
      if (outfileinfo->outfpbwttab != NULL)
      {
        if (fwrite(&cc,sizeof (Uchar),(size_t) 1,outfileinfo->outfpbwttab)
                    != (size_t) 1)
        {
          env_error_set(env,"cannot write 1 item of size %lu: "
                            "errormsg=\"%s\"",
                          (unsigned long) sizeof (Uchar),
                          strerror(errno));
          haserr = true;
          break;
        }
      }
      if (outfileinfo->outfplcptab != NULL)
      {
        lcpvalue = nextLcpvalueiterator(outfileinfo->lvi,
                                        (outfileinfo->pageoffset == 0)
                                          ? true : false,
                                        suftab,
                                        numberofsuffixes);
        if (outlcpvalue(lcpvalue,pos,outfileinfo,env) != 0)
        {
          haserr = true;
          break;
        }
        if (outfileinfo->maxbranchdepth < lcpvalue)
        {
          outfileinfo->maxbranchdepth = lcpvalue;
        }
      }
    }
  }
  outfileinfo->pageoffset += numberofsuffixes;
  return haserr ? -1 : 0;
}

extern const Uchar *
SfxIReadESQRange(sfxInterface *iface, Seqpos start, Seqpos len,
                 Uchar *dest)
{
  size_t i;
  assert(dest);
  for(i = 0; i < (size_t)len; ++i)
  {
    dest[i] = getencodedchar(iface->outfileinfo.encseq, start + i,
                            iface->so.readmode);
  }
  return dest;
}


extern int
SfxIRegisterReader(sfxInterface *iface, listenerID *id,
                   enum sfxDataRequest request, Env *env)
{
  size_t numReaders = iface->numReaders++;
  iface->readers = env_ma_realloc(
    env, iface->readers, sizeof(iface->readers[0]) * (numReaders + 1));
  iface->readers[numReaders].nextReadPos = 0;
  iface->readers[numReaders].readFlag = request;
  iface->allRequests |= request;
  return 1;
}

static Seqpos
getSufTabVal(sfxInterface *iface, Seqpos pos, Env *env)
{
  assert(iface && pos < iface->length);
  while(1)
  {
    if(pos >= iface->lastGeneratedStart
       && iface->lastGeneratedLen)
    {
      /* pos is in last block */
      return iface->lastGeneratedSufTabSegment[pos
                                               - iface->lastGeneratedStart];
    }
    else if(pos >= iface->prevGeneratedStart
            && iface->prevGeneratedLen)
    {
      /* pos is in prev block */
      return iface->prevGeneratedSufTabSegments[pos
                                                - iface->prevGeneratedStart];
    }
    else
    {
      while(pos >= iface->lastGeneratedStart + iface->lastGeneratedLen)
        if(!SfxIReadAdvance(iface, pos, env))
          break;
    }    
  }
}

extern size_t
SfxIReadBWTRange(sfxInterface *iface, listenerID id, size_t len,
                 Uchar *dest, Env *env)
{
  size_t i, effLen;
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  assert(iface->readers[id].readFlag == SFX_REQUEST_BWTTAB);
  start = iface->readers[id].nextReadPos;
  effLen = MIN((size_t)len, (size_t)(iface->length - start));
  for(i = 0; i < (size_t)len; ++i)
  {
    dest[i] = getencodedchar(iface->outfileinfo.encseq, 
                             getSufTabVal(iface, start + i, env),
                             iface->so.readmode);
  }
  iface->readers[id].nextReadPos = start + len;
  return effLen;
}

extern size_t
SfxIReadLCPRange(sfxInterface *iface, listenerID id, size_t len,
                 Seqpos *dest, Env *env)
{
  size_t i, effLen;
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  assert(iface->readers[id].readFlag == SFX_REQUEST_LCPTAB);
  start = iface->readers[id].nextReadPos;
  effLen = MIN(len, iface->length - start);
  if(start == 0)
    dest[0] = 0, i = 1;
  else  
    i = 0;
  for(; i < (size_t)len; ++i)
  {
#ifndef NDEBUG
    int cmp =
#endif /* NDEBUG */
      comparetwosuffixes(iface->outfileinfo.encseq,
                         iface->so.readmode,
                         dest + i,
                         false,
                         false,
                         0,
                         getSufTabVal(iface, start + i - 1, env),
                         getSufTabVal(iface, start + i, env));
#ifndef NDEBUG
    if(cmp > 0)
    {
      fprintf(stderr,"pos = " FormatSeqpos
              ": cmp " FormatSeqpos
              " " FormatSeqpos " = %d",
              PRINTSeqposcast(start + (Seqpos)i - 1),
              PRINTSeqposcast(start + (Seqpos)i),
              PRINTSeqposcast(getSufTabVal(iface, start + (Seqpos)i, env)),
              cmp);
      exit(EXIT_FAILURE);
    }
#endif /* NDEBUG */
  }
  iface->readers[id].nextReadPos = start + len - 1;
  return effLen;
}



extern size_t
SfxIReadSufTabRange(sfxInterface *iface, listenerID id, size_t len,
                    Seqpos *dest, Env *env)
{
  Seqpos start;
  assert(iface && id < iface->numReaders && dest);
  start = iface->readers[id].nextReadPos;
  while(1)
  {
    if(start >= iface->lastGeneratedStart
       && len <= iface->lastGeneratedLen)
    {
      /* chunk is completely contained in last block */
      iface->readers[id].nextReadPos = start + len;
      memcpy(dest,
             iface->lastGeneratedSufTabSegment
             + start - iface->lastGeneratedStart,
             sizeof(dest[0]) * len);
      return len;
    }
    else if(start >= iface->prevGeneratedStart
            && len <= iface->prevGeneratedLen)
    {
      /* chunk is completely contained in prev block */
      iface->readers[id].nextReadPos = start + len;
      memcpy(dest,
             iface->prevGeneratedSufTabSegments
             + start - iface->prevGeneratedStart,
             sizeof(dest[0]) * len);
      return len;
    }
    else if(start + len < iface->lastGeneratedStart + iface->lastGeneratedLen)
    {
      /* chunk overlaps both segments, but does not extend beyond
       * last segment */
      memcpy(dest,
             iface->prevGeneratedSufTabSegments
             + start - iface->prevGeneratedStart,
             sizeof(dest[0])
             * (iface->prevGeneratedLen + iface->prevGeneratedStart - start));
      memcpy(dest,
             iface->lastGeneratedSufTabSegment,
             sizeof(dest[0])
             * (start + len - iface->lastGeneratedStart));
      iface->readers[id].nextReadPos = start + len;
      return len;
    }
    else
    {
      while(start + len > iface->prevGeneratedStart + iface->prevGeneratedLen)
        if(!SfxIReadAdvance(iface, start + len, env))
        {
          /* Caution: sneakily updates the length parameter to obtain a
           * request we can fulfill */
          len = iface->lastGeneratedStart + iface->lastGeneratedLen - start;
          break;
        }
    }
  }
}

static Seqpos
findMinOpenRequest(sfxInterface *iface, int reqType)
{
  Seqpos min;
  size_t i;
  assert(iface);
  min = iface->length;
  for(i = 0; i < iface->numReaders; ++i)
    if(iface->readers[i].readFlag & reqType)
      min = MIN(min, iface->readers[i].nextReadPos);
  return min;
}

static Seqpos
SfxIReadAdvance(sfxInterface *iface,
                Seqpos requestMaxPos,
                Env *env)
{
  Seqpos requestMinPos, lengthOfExtension = 0;
  assert(iface && env);
  /* 1. find first position that still needs to be read */
  requestMinPos = findMinOpenRequest(iface, SFX_REQUEST_ANY);
  /* 2. move still unread old values as far as possible to head of copy */
  /* FIXME: try to optimize if requestMinPos is beyond end of prevSegments */
  assert(requestMinPos >= iface->prevGeneratedStart);
  memmove(iface->prevGeneratedSufTabSegments,
          iface->prevGeneratedSufTabSegments + requestMinPos
          - iface->prevGeneratedStart,
          sizeof(iface->prevGeneratedSufTabSegments[0]) *
          (iface->prevGeneratedLen - requestMinPos
           + iface->prevGeneratedStart));
  iface->prevGeneratedStart = requestMinPos;
  iface->prevGeneratedLen =
    iface->prevGeneratedLen - requestMinPos + iface->prevGeneratedStart;
  do
  {
    /* 3. extend cache to also accept all values in last region read */
    /* FIXME: only copy parts actually required */
    if(iface->lastGeneratedLen)
    {
      size_t prevGeneratedSizeLeft = iface->prevGeneratedSize 
        - iface->prevGeneratedLen;
      if(iface->lastGeneratedLen > prevGeneratedSizeLeft)
      {
        iface->prevGeneratedSufTabSegments =
          env_ma_realloc(env, iface->prevGeneratedSufTabSegments,
                         sizeof(iface->lastGeneratedSufTabSegment[0]) *
                         (iface->prevGeneratedLen + iface->lastGeneratedLen));
        iface->prevGeneratedSize =
          iface->prevGeneratedLen + iface->lastGeneratedLen;
      }
      memcpy(iface->prevGeneratedSufTabSegments + iface->prevGeneratedLen,
             iface->lastGeneratedSufTabSegment,
             sizeof(iface->lastGeneratedSufTabSegment[0])
             * iface->lastGeneratedLen);
      iface->prevGeneratedLen += iface->lastGeneratedLen;
    }
    /* 4. read next region of sequence by calling nextSfxIterator */
    if((iface->lastGeneratedSufTabSegment = 
        nextSfxiterator(&iface->lastGeneratedLen, &iface->specialsuffixes,
                        NULL, iface->sfi, env)))
      lengthOfExtension += iface->lastGeneratedLen;
    if(iface->lastGeneratedSufTabSegment == NULL
       || suftab2file(&iface->outfileinfo, iface->lastGeneratedSufTabSegment,
                      iface->so.readmode,
                      iface->lastGeneratedLen, env) != 0)
      break;
    /* 5. if positions in region don't suffice go back to step 3. */
  } while(requestMaxPos >= iface->lastGeneratedStart + iface->lastGeneratedLen);
  return lengthOfExtension;
}


