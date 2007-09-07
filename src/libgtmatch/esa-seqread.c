#include <limits.h>
#include "sarr-def.h"
#include "spacedef.h"
#include "esa-seqread.h"

#include "sfx-map.pr"

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

 struct Sequentialsuffixarrayreader
{
  Suffixarray suffixarray;
  Seqpos totallength,
         nextsuftabindex, /* only relevant if mapped == true */
         nextlcptabindex, /* only relevant if mapped == true */
         largelcpindex;   /* only relevant if mapped == true */
  bool mapped;
};

Sequentialsuffixarrayreader *newSequentialsuffixarrayreader(
                                        const Str *indexname,
                                        unsigned int demand,
                                        bool mapped,
                                        Env *env)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  if ((mapped ? mapsuffixarray : streamsuffixarray)(&ssar->suffixarray,
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
  ssar->mapped = mapped;
  return ssar;
}

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar,
                                     Env *env)
{
  freesuffixarray(&((*ssar)->suffixarray),env);
  FREESPACE(*ssar);
}

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Env *env)
{
  Uchar tmpsmalllcpvalue;
  int retval;

  if (ssar->mapped)
  {
    if (ssar->nextlcptabindex > ssar->totallength)
    {
      return 0;
    }
    tmpsmalllcpvalue = ssar->suffixarray.lcptab[ssar->nextlcptabindex++];
  } else
  {
    retval = readnextUcharfromstream(&tmpsmalllcpvalue,
				     &ssar->suffixarray.lcptabstream,
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

    if (ssar->mapped)
    {
      assert(ssar->suffixarray.llvtab[ssar->largelcpindex].position ==
             ssar->nextlcptabindex-1);
      *currentlcp = ssar->suffixarray.llvtab[ssar->largelcpindex++].value;
    } else
    {
      retval = readnextLargelcpvaluefromstream(&tmpexception,
					       &ssar->suffixarray.llvtabstream,
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
  return 1;
}

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Env *env)
{
  if(ssar->mapped)
  {
    *currentsuffix = ssar->suffixarray.suftab[ssar->nextsuftabindex++];
    return 1;
  }
  return readnextSeqposfromstream(currentsuffix,
				  &ssar->suffixarray.suftabstream,
				  env);
}

Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr)
{
  return sarr->suffixarray.encseq;
}

Readmode readmodeSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr)
{
  return sarr->suffixarray.readmode;
}

const Alphabet *alphabetSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr)
{
  return sarr->suffixarray.alpha;
}

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr)
{
  return sarr->suffixarray.numofdbsequences;
}
