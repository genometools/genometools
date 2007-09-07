#include <limits.h>
#include "sarr-def.h"
#include "spacedef.h"
#include "esa-seqread.h"
#include "sfx-lcpval.h"

#include "sfx-map.pr"

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

 struct Sequentialsuffixarrayreader
{
  Suffixarray *suffixarray;
  Seqpos totallength,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
  Sequentialaccesstype seqactype;
  Lcpvalueiterator *lvi;
  const Seqpos *suftab;
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos numberofsuffixes;
};

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                        const Str *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        Env *env)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ALLOCASSIGNSPACE(ssar->suffixarray,NULL,Suffixarray,1);
  assert(seqactype == SEQ_mappedboth || seqactype == SEQ_scan);
  if (((seqactype == SEQ_mappedboth) 
         ? mapsuffixarray : streamsuffixarray)(ssar->suffixarray,
                                               &ssar->totallength,
                                               demand,
                                               indexname,
                                               false,
                                               env) != 0)
  {
    FREESPACE(ssar);
    return NULL;
  }
  ssar->nextsuftabindex = 0;
  ssar->nextlcptabindex = (Seqpos) 1;
  ssar->largelcpindex = 0;
  ssar->seqactype = seqactype;
  ssar->suftab = NULL;
  ssar->numberofsuffixes = 0;
  ssar->encseq = ssar->suffixarray->encseq;
  ssar->readmode = ssar->suffixarray->readmode;
  return ssar;
}

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromRAM(
                                        const Encodedsequence *encseq,
                                        const Seqpos *suftab,
                                        Seqpos numberofsuffixes,
                                        Readmode readmode,
                                        Env *env)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ssar->lvi = newLcpvalueiterator(encseq,readmode,env);
  (void) nextLcpvalueiterator(ssar->lvi,
                              true,
                              suftab,
                              numberofsuffixes);
  ssar->suffixarray = NULL;
  ssar->nextsuftabindex = 0;
  ssar->nextlcptabindex = (Seqpos) 1; /* not required here */
  ssar->largelcpindex = 0; /* not required here */
  ssar->seqactype = SEQ_suftabfrommemory;
  ssar->suftab = suftab;
  ssar->numberofsuffixes = numberofsuffixes;
  ssar->readmode = readmode;
  ssar->encseq = encseq;
  return ssar;
}

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar,
                                     Env *env)
{
  if((*ssar)->suffixarray != NULL)
  {
    freesuffixarray((*ssar)->suffixarray,env);
    FREESPACE((*ssar)->suffixarray);
  }
  FREESPACE(*ssar);
}

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Env *env)
{
  Uchar tmpsmalllcpvalue;
  int retval;

  if(ssar->seqactype == SEQ_suftabfrommemory)
  {
    *currentlcp = nextLcpvalueiterator(ssar->lvi,
                                       true,
                                       ssar->suftab,
                                       ssar->numberofsuffixes);
  } else
  {
    if (ssar->seqactype == SEQ_mappedboth)
    {
      if (ssar->nextlcptabindex > ssar->totallength)
      {
        return 0;
      }
      tmpsmalllcpvalue = ssar->suffixarray->lcptab[ssar->nextlcptabindex++];
    } else
    {
      retval = readnextUcharfromstream(&tmpsmalllcpvalue,
                                       &ssar->suffixarray->lcptabstream,
                                       env);
      if (retval < 0)
      {
        return -1;
      }
      if (retval == 0)
      {
        return 0;
      }
    }
    if (tmpsmalllcpvalue == (Uchar) UCHAR_MAX)
    {
      Largelcpvalue tmpexception;

      if (ssar->seqactype == SEQ_suftabfrommemory)
      {
        assert(ssar->suffixarray->llvtab[ssar->largelcpindex].position ==
               ssar->nextlcptabindex-1);
        *currentlcp = ssar->suffixarray->llvtab[ssar->largelcpindex++].value;
      } else
      {
        retval = readnextLargelcpvaluefromstream(
                                          &tmpexception,
                                          &ssar->suffixarray->llvtabstream,
                                          env);
        if (retval < 0)
        {
          return -1;
        }
        if (retval == 0)
        {
          env_error_set(env,"file %s: line %d: unexpected end of file when "
                        "reading llvtab",__FILE__,__LINE__);
          return -1;
        }
        *currentlcp = tmpexception.value;
      }
    } else
    {
      *currentlcp = (Seqpos) tmpsmalllcpvalue;
    }
  }
  return 1;
}

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Env *env)
{
  if(ssar->seqactype == SEQ_scan)
  {
    return readnextSeqposfromstream(currentsuffix,
                                    &ssar->suffixarray->suftabstream,
                                    env);
  }
  if (ssar->seqactype == SEQ_mappedboth)
  {
    *currentsuffix = ssar->suffixarray->suftab[ssar->nextsuftabindex++];
    return 1;
  }
  *currentsuffix = ssar->suftab[ssar->nextsuftabindex++];
  return 1;
}

const Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr)
{
  return sarr->encseq;
}

Readmode readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr)
{
  return sarr->readmode;
}

const Alphabet *alphabetSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr)
{
  assert(sarr->suffixarray != NULL);
  return sarr->suffixarray->alpha;
}

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr)
{
  assert(sarr->suffixarray != NULL);
  return sarr->suffixarray->numofdbsequences;
}
