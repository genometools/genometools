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

#ifndef INLINEDSequentialsuffixarrayreader
#include <limits.h>
#include "seqpos-def.h"
#include "measure-time-if.h"
#include "esa-seqread.h"
#include "sfx-suffixer.h"
#include "verbose-def.h"

#include "sfx-apfxlen.pr"
#include "esa-maxpairs.pr"

typedef struct
{
  unsigned int minlength;
  Encodedsequence *encseq;
  int (*processmaxmatch)(void *,Seqpos,Seqpos,Seqpos,Env *);
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
                 Verboseinfo *verboseinfo,
                 Env *env)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false;
  Sfxiterator *sfi;
  bool specialsuffixes = false;

  sfi = newSfxiterator(specialcharacters,
                       specialranges,
                       ssi->encseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       mtime,
                       NULL,
                       env);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    Sequentialsuffixarrayreader *ssar = NULL;
    bool firstpage = true;

    ssar = newSequentialsuffixarrayreaderfromRAM(ssi->encseq,
                                                 readmode,
                                                 env);
    while (true)
    {
      suftabptr = nextSfxiterator(&numberofsuffixes,&specialsuffixes,
                                  mtime,sfi,env);
      if (suftabptr == NULL || specialsuffixes)
      {
        break;
      }
      updateSequentialsuffixarrayreaderfromRAM(ssar,
                                               suftabptr,
                                               firstpage,
                                               numberofsuffixes);
      firstpage = false;
      if (enumeratemaxpairs(ssar,
                            numofchars,
                            ssi->encseq,
                            readmode,
                            ssi->minlength,
                            ssi->processmaxmatch,
                            ssi->processmaxmatchinfo,
                            verboseinfo,
                            env) != 0)
      {
        haserr = true;
      }
    }
    if (ssar != NULL)
    {
      freeSequentialsuffixarrayreader(&ssar,env);
    }
  }
  if (sfi != NULL)
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
                                                  Seqpos,Seqpos,Env *),
                           void *processmaxmatchinfo,
                           Verboseinfo *verboseinfo,
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
                                     getmapsizeAlphabet(alpha),
                                     verboseinfo,
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
                                  verboseinfo,
                                  env) != 0)
  {
    haserr = true;
  }
  freeEncodedsequence(&ssi.encseq,env);
  return haserr ? -1 : 0;
}
#endif /* ifndef INLINEDSequentialsuffixarrayreader */
