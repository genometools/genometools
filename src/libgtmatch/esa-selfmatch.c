#include <limits.h>
#include "seqpos-def.h"
#include "measure-time-if.h"
#include "esa-seqread.h"
#include "sfx-suffixer.h"

#include "sfx-apfxlen.pr"
#include "esa-maxpairs.pr"

typedef struct
{
  unsigned int minlength;
  Encodedsequence *encseq;
  int (*processmaxmatch)(void *,Seqpos,Seqpos,Seqpos);
  void *processmaxmatchinfo;
} Substringmatchinfo;

static int constructsarrandrunmaxpairs(
                 Substringmatchinfo *ssi,
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 Readmode readmode,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 Measuretime *mtime,
                 Env *env)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false;
  Sfxiterator *sfi;
  bool specialsuffixes = false;
  Sequentialsuffixarrayreader *ssar = NULL;

  sfi = newsfxiterator(specialcharacters,
                       specialranges,
                       ssi->encseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       mtime,
                       env);
  if(sfi == NULL)
  {
    haserr = true;
  } else
  {
    suftabptr = nextSfxiterator(&numberofsuffixes,&specialsuffixes,
                                mtime,sfi,env);
    assert(suftabptr != NULL);
    ssar = newSequentialsuffixarrayreaderfromRAM(ssi->encseq,
                                                 suftabptr,
                                                 numberofsuffixes,
                                                 readmode,
                                                 env);
    if(enumeratemaxpairs(ssar,
                         numofchars,
                         ssi->encseq,
                         readmode,
                         ssi->minlength,
                         ssi->processmaxmatch,
                         ssi->processmaxmatchinfo,
                         env) != 0)
    {
      haserr = true;
    }
  }
  if(ssar != NULL)
  {
    freeSequentialsuffixarrayreader(&ssar,env);
  }
  if(sfi != NULL)
  {
    freeSfxiterator(&sfi,env);
  }
  return haserr ? -1 : 0;
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
  if (constructsarrandrunmaxpairs(&ssi,
                                  samplespecialcharinfo.specialcharacters,
                                  samplespecialcharinfo.specialranges,
                                  Forwardmode,
                                  numofchars,
                                  recommendedprefixlength(numofchars,
                                                          dblen+querylen+1),
                                  (unsigned int) 1, /* parts */
                                  NULL,
		                  env) != 0)
  {
    haserr = true;
  }
  freeEncodedsequence(&ssi.encseq,env);
  return haserr ? -1 : 0;
}
