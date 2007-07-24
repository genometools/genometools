/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "emimergeesa.h"
#include "esafileend.h"
#include "fmindex.h"
#include "divmodmul.h"
#include "chardef.h"

#include "mergeesa.pr"
#include "sfxmap.pr"
#include "encseq2offset.pr"
#include "mkidxfilecpy.pr"
#include "fmkeyval.pr"
#include "mapspecfm.pr"
#include "alphabet.pr"
#include "opensfxfile.pr"

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Seqpos);

static void allocatefmtables(Fmindex *fm,bool storeindexpos,Env *env)
{
  ALLOCASSIGNSPACE (fm->tfreq, NULL, Seqpos,TFREQSIZE(fm->mapsize));
  ALLOCASSIGNSPACE (fm->superbfreq, NULL, Seqpos ,
                    SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks));
  if(storeindexpos)
  {
    ALLOCASSIGNSPACE (fm->markpostable,NULL,Seqpos,
                      MARKPOSTABLELENGTH(fm->bwtlength,fm->markdist));
  } else
  {
    fm->markpostable = NULL;
  }
  fm->specpos.nextfreePairBwtidx = 0;
  fm->specpos.allocatedPairBwtidx
    = (unsigned long) determinenumberofspecialstostore(&fm->specialcharinfo);
  printf("# %lu wildcards in the last " FormatSeqpos 
         " characters (%.2f)\n",
          (unsigned long) (fm->specialcharinfo.specialcharacters - 
                          fm->specpos.allocatedPairBwtidx),
          PRINTSeqposcast(fm->specialcharinfo.specialcharacters),
          (double) (fm->specialcharinfo.specialcharacters - 
                    fm->specpos.allocatedPairBwtidx)/
                        fm->specialcharinfo.specialcharacters);
  ALLOCASSIGNSPACE(fm->specpos.spacePairBwtidx,NULL,PairBwtidx,
                   fm->specpos.allocatedPairBwtidx);
  ALLOCASSIGNSPACE (fm->bfreq, NULL, Uchar,
                    BFREQSIZE(fm->mapsize,fm->nofblocks));
}

static void set0frequencies(Fmindex *fm)
{
  Seqpos i;

  for (i = 0; i < TFREQSIZE(fm->mapsize); i++)
  {
    fm->tfreq[i] = 0;
  }
  for (i = 0; i < BFREQSIZE(fm->mapsize,fm->nofblocks); i++)
  {
    fm->bfreq[i] = 0;
  }
  for (i = 0; i < SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks); i++)
  {
    fm->superbfreq[i] = 0;
  }
}

static void finalizefmfrequencies(Fmindex *fm)
{
  uint32_t j; 
  Seqpos i, *freqptr;

  for (j = (uint32_t) 2; j <= fm->mapsize; j++)
  {
    fm->tfreq[j] += fm->tfreq[j - 1];
  }
  freqptr = fm->superbfreq;
  for (j = 0; j < fm->mapsize; j++)
  {
    for (i = (Seqpos) 2; i < fm->nofsuperblocks; i++)
    {
      freqptr[i] += freqptr[i-1];
    }
    freqptr += fm->nofsuperblocks;
  }
}

static void showconstructionmessage(const Str *indexname,
                                    Seqpos totallength,
                                    unsigned long fmsize,
                                    uint32_t log2bsize,
                                    uint32_t log2markdist,
                                    uint32_t mapsize)
{
  printf("# construct fmindex \"%s\" for bsize=%lu, superbsize=%lu,",
          str_get(indexname),
          (unsigned long) POW2(log2bsize),
          (unsigned long) POW2(log2markdist));
  printf(" len=" FormatSeqpos ", alphasize=%u: size ",
          PRINTSeqposcast(totallength),
          (unsigned int) (mapsize-1));
  printf("%lu bytes, space overhead %.2f\n",
          fmsize,
          (double) fmsize/(double) (totallength+1)); 
}

static int nextesamergedsufbwttabvalues(DefinedSeqpos *longest,
                                        Uchar *bwtvalue,
                                        Seqpos *suftabvalue,
                                        Emissionmergedesa *emmesa,
                                        const Seqpos *sequenceoffsettable,
                                        Seqpos bwtpos,
                                        Env *env)
{
  Indexedsuffix indexedsuffix;

  if(emmesa->buf.nextaccessidx >= emmesa->buf.nextstoreidx)
  {
    if(emmesa->numofentries == 0)
    {
      return 0;
    }
    if(stepdeleteandinsertothersuffixes(emmesa,env) != 0)
    {
      return -1;
    }
    if(emmesa->buf.nextstoreidx == 0)
    {
      return 0;
    }
    emmesa->buf.nextaccessidx = 0;
  }
  indexedsuffix = emmesa->buf.suftabstore[emmesa->buf.nextaccessidx];
  *suftabvalue = sequenceoffsettable[indexedsuffix.idx] + 
                 indexedsuffix.startpos;
  if(indexedsuffix.startpos == 0)
  {
    if(indexedsuffix.idx == 0)
    {
      if(longest->defined)
      {
        env_error_set(env,"longest is already defined as " FormatSeqpos,
                      longest->value);
        return -2;
      }
      longest->defined = true;
      longest->value = bwtpos;
      *bwtvalue = (Uchar) UNDEFBWTCHAR;
    } else
    {
      *bwtvalue = (Uchar) SEPARATOR;
    }
  } else
  {
    *bwtvalue = getencodedchar(
                            emmesa->suffixarraytable[indexedsuffix.idx].encseq,
                            indexedsuffix.startpos-1);
  }
  emmesa->buf.nextaccessidx++;
  return 1;
}

int sufbwt2fmindex(Fmindex *fm,
                   uint32_t log2bsize,
                   uint32_t log2markdist,
                   const Str *outfmindex,
                   const StrArray *indexnametab,
                   bool storeindexpos,
                   Env *env)
{
  Suffixarray suffixarray;
  Emissionmergedesa emmesa;
  Uchar cc;
  Seqpos bwtpos, 
         totallength, 
         suftabvalue, 
         *sequenceoffsettable = NULL,
         firstignorespecial = 0,
         nextmark, 
         *markptr, 
         nextprogress, 
         tmpsuftabvalue, 
         stepprogress;
  uint32_t mapsize = 0, 
           suffixlength = 0,
           numofindexes; 
  int retval;
  DefinedSeqpos longest = { false, 0 };
  PairBwtidx *pairptr;
  FILE *outbwt = NULL;
  Str *tmpfilename = NULL;
  Specialcharinfo specialcharinfo;
  bool haserr = false;
  
  numofindexes = (uint32_t) strarray_size(indexnametab);
  if(numofindexes == (uint32_t) 1)
  {
    Str *indexname = str_new_cstr(strarray_get(indexnametab,0),env);
    if(streamsuffixarray(&suffixarray,
                         SARR_SUFTAB | SARR_BWTTAB,
                         indexname,
                         env) != 0)
    {
      haserr = true;
    }
    if(!haserr)
    {
      totallength = getencseqtotallength(suffixarray.encseq);
      mapsize = getmapsizeAlphabet(suffixarray.alpha);
      specialcharinfo = suffixarray.specialcharinfo;
      firstignorespecial = totallength - specialcharinfo.specialcharacters;
      if(makeindexfilecopy(outfmindex,indexname,ALPHATABSUFFIX,0,env) != 0)
      {
        haserr = true;
      }
    }
    if(!haserr)
    {
      if(makeindexfilecopy(outfmindex,
                           indexname,
                           BWTTABSUFFIX,
                           firstignorespecial,
                           env) != 0)
      {
        haserr = true;
      }
    }
    str_delete(indexname,env);
  } else
  {
    if(initEmissionmergedesa(&emmesa,
                             indexnametab,
                             SARR_ESQTAB | SARR_SUFTAB | SARR_LCPTAB,
                             env) != 0)
    {
      haserr = true;
    }
    if(!haserr)
    {
      Str *indexname = str_new_cstr(strarray_get(indexnametab,0),env);
      suffixlength = 0;
      if(makeindexfilecopy(outfmindex,indexname,ALPHATABSUFFIX,0,env) != 0)
      {
        haserr = true;
      }
      str_delete(indexname,env);
    }
    if(!haserr)
    {
      sequenceoffsettable = encseqtable2seqoffsets(&totallength,
                                                   &specialcharinfo,
                                                   emmesa.suffixarraytable,
                                                   numofindexes,
                                                   env);
      if(sequenceoffsettable == NULL)
      {
        haserr = true;
      }
    }
    if(!haserr)
    {
      longest.defined = false;
      longest.value = 0;
      outbwt = opensfxfile(outfmindex,BWTTABSUFFIX,"wb",env);
      if (outbwt == NULL)
      {
        haserr = true;
      }
    }
    if(!haserr)
    {
      mapsize = getmapsizeAlphabet(emmesa.alpha);
      firstignorespecial = totallength - specialcharinfo.specialcharacters;
    }
  }
  if(!haserr)
  {
    printf("# firstignorespecial=" FormatSeqpos "\n",
              PRINTSeqposcast(firstignorespecial));
    computefmkeyvalues (fm,
                        totallength+1,
                        log2bsize,
                        log2markdist,
                        mapsize,
                        suffixlength,
                        storeindexpos,
                        &specialcharinfo);
    showconstructionmessage(outfmindex,
                            totallength,
                            fm->sizeofindex,
                            log2bsize,
                            log2markdist,
                            mapsize);
    allocatefmtables(fm,storeindexpos,env);
    set0frequencies(fm);
    if(storeindexpos)
    {
      markptr = fm->markpostable;
    } else
    {
      markptr = NULL;
    }
    nextprogress = stepprogress = totallength/78;
    for(bwtpos = 0, nextmark = 0; ; bwtpos++)
    {
      if(numofindexes == (uint32_t) 1)
      {
        retval = readnextSeqposfromstream(&tmpsuftabvalue,
                                          &suffixarray.suftabstream,
                                          env);
        if(retval < 0)
        {
          haserr = true;
          break;
        } 
        if(retval == 0)
        {
          break;
        }
        suftabvalue = (Seqpos) tmpsuftabvalue;
        retval = readnextUcharfromstream(&cc,&suffixarray.bwttabstream,env);
        if(retval < 0)
        {
          haserr = true;
          break;
        }
        if(retval == 0)
        {
          break;
        }
      } else
      {
        retval = nextesamergedsufbwttabvalues(&longest,
                                              &cc,
                                              &suftabvalue,
                                              &emmesa,
                                              sequenceoffsettable,
                                              bwtpos,
                                              env);
        if(retval < 0)
        {
          haserr = true;
          break;
        }
        if(retval == 0)
        {
          break;
        }
        if(fwrite(&cc,
                  sizeof(Uchar),
                  (size_t) 1,
                  outbwt) != (size_t) 1)
        {
          haserr = true;
          break;
        }
      }
      if(bwtpos == nextprogress)
      {
        if(bwtpos == stepprogress)
        {
          (void) putchar('#');
        }
        (void) putchar('.');
        (void) fflush(stdout);
        nextprogress += stepprogress;
      } 
      if(storeindexpos && bwtpos == nextmark)
      {
        *markptr++ = suftabvalue;
        nextmark += fm->markdist;
      }
      if(ISBWTSPECIAL(cc))
      {
        if (bwtpos < firstignorespecial)
        {
          if(storeindexpos)
          {
            pairptr = fm->specpos.spacePairBwtidx +
                      fm->specpos.nextfreePairBwtidx++;
            if(pairptr >= fm->specpos.spacePairBwtidx + 
                          fm->specpos.allocatedPairBwtidx)
            {
              env_error_set(env,"program error: not enough space for specpos");
              haserr = true;
              break;
            }
            pairptr->bwtpos = bwtpos;
            pairptr->suftabvalue = suftabvalue;
          }
        }
      } else
      {
        fm->tfreq[cc+1]++;
        fm->bfreq[(cc * fm->nofblocks) + (bwtpos >> fm->log2bsize)]++;
        fm->superbfreq[(cc * fm->nofsuperblocks) +
                       (bwtpos >> fm->log2superbsize) + 1]++;
      }
    }
  }
  if(!haserr)
  {
    if(storeindexpos && 
       fm->specpos.allocatedPairBwtidx != fm->specpos.nextfreePairBwtidx)
    {
      env_error_set(env,"program error: too much space for specpos: "
             "allocated = %lu != %lu = used",
             fm->specpos.allocatedPairBwtidx,
             fm->specpos.nextfreePairBwtidx);
      haserr = true;
    }
  }
  if(!haserr)
  {
    (void) putchar('\n');
    finalizefmfrequencies(fm);
    if(fm->suffixlength > 0)
    {
      ALLOCASSIGNSPACE(fm->boundarray,NULL,Bwtbound,fm->numofcodes);
    }
    if(numofindexes == (uint32_t) 1)
    {
      fm->longestsuffixpos = suffixarray.longest.value;
      freesuffixarray(&suffixarray,env);
    } else
    {
      if(!longest.defined)
      {
        env_error_set(env,"longest is not defined after merging");
        haserr = true;
      }
      if(!haserr)
      {
        fm->longestsuffixpos = longest.value;
      }
      env_fa_xfclose(outbwt, env);
      wraptEmissionmergedesa(&emmesa,env);
    }
  }
  FREESPACE(sequenceoffsettable);
  FREESPACE(tmpfilename);
  return haserr ? -1 : 0;
}
