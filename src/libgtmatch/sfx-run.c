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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include "libgtcore/chardef.h"
#include "libgtcore/filelengthvalues.h"
#include "libgtcore/seqiterator.h"
#include "spacedef.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esafileend.h"
#include "intcode-def.h"
#include "sfx-suffixer.h"
#include "sfx-lcpval.h"
#include "readmode-def.h"
#include "verbose-def.h"
#include "stamp.h"

#include "opensfxfile.pr"
#include "fillsci.pr"
#include "sfx-opt.pr"
#include "sfx-outprj.pr"
#include "sfx-apfxlen.pr"

typedef struct
{
  FILE *outfpsuftab,
       *outfplcptab,
       *outfpllvtab,
       *outfpbwttab;
  Seqpos pageoffset,
         numoflargelcpvalues,
         maxbranchdepth;
  const Encodedsequence *encseq;
  DefinedSeqpos longest;
  Lcpvalueiterator *lvi;
} Outfileinfo;

static void initoutfileinfo(Outfileinfo *outfileinfo)
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

static int outal1file(const Str *indexname,const Alphabet *alpha,Env *env)
{
  FILE *al1fp;
  bool haserr = false;

  env_error_check(env);
  al1fp = opensfxfile(indexname,ALPHABETFILESUFFIX,"wb",env);
  if (al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    outputalphabet(al1fp,alpha);
    env_fa_xfclose(al1fp,env);
  }
  return haserr ? -1 : 0;
}

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)\
        if (!haserr && (FLAG))\
        {\
          outfileinfo.encseq = encseq;\
          PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",env);\
          if ((PTR) == NULL)\
          {\
            haserr = true;\
          }\
        }

static int suffixeratorwithoutput(
                 Outfileinfo *outfileinfo,
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 const Encodedsequence *encseq,
                 Readmode readmode,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 Measuretime *mtime,
                 Verboseinfo *verboseinfo,
                 Env *env)
{
  const Seqpos *suftabptr;
  Seqpos len;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(specialcharacters,
                       specialranges,
                       encseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       mtime,
                       verboseinfo,
                       env);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    while (true)
    {
      suftabptr = nextSfxiterator(&len,&specialsuffixes,mtime,sfi,env);
      if (suftabptr == NULL)
      {
        break;
      }
      if (suftab2file(outfileinfo,suftabptr,readmode,len,env) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi,env);
  }
  return haserr ? -1 : 0;
}

static unsigned long *initcharacterdistribution(const Alphabet *alpha,Env *env)
{
  unsigned long *characterdistribution;
  unsigned int mapsize, idx;

  mapsize = getmapsizeAlphabet(alpha);
  ALLOCASSIGNSPACE(characterdistribution,NULL,unsigned long,mapsize-1);
  for (idx=0; idx<mapsize-1; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
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

static int runsuffixerator(Suffixeratoroptions *so,Verboseinfo *verboseinfo,
                           Env *env)
{
  unsigned int numofchars = 0;
  unsigned long numofsequences;
  Seqpos totallength;
  Alphabet *alpha;
  Specialcharinfo specialcharinfo;
  Filelengthvalues *filelengthtab = NULL;
  bool haserr = false;
  Encodedsequence *encseq = NULL;
  Measuretime *mtime;
  Outfileinfo outfileinfo;
  unsigned long *characterdistribution = NULL;

  env_error_check(env);
  inittheclock(&mtime,
               "determining sequence length and number of special symbols",
               env);
  alpha = assigninputalphabet(so->isdna,
                              so->isprotein,
                              so->str_smap,
                              so->filenametab,
                              env);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (!so->isplain)
    {
      characterdistribution = initcharacterdistribution(alpha,env);
    }
    if (fasta2sequencekeyvalues(so->str_indexname,
                                &numofsequences,
                                &totallength,
                                &specialcharinfo,
                                so->filenametab,
                                &filelengthtab,
                                getsymbolmapAlphabet(alpha),
                                so->isplain,
                                so->outdestab,
                                characterdistribution,
                                verboseinfo,
                                env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    numofchars = getnumofcharsAlphabet(alpha);
    if (outal1file(so->str_indexname,alpha,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    deliverthetime(stdout,mtime,"computing sequence encoding",env);
    encseq = files2encodedsequence(true,
                                   so->filenametab,
                                   so->isplain,
                                   totallength,
                                   &specialcharinfo,
                                   alpha,
                                   str_length(so->str_sat) > 0
                                         ? str_get(so->str_sat)
                                         : NULL,
                                   verboseinfo,
                                   env);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,encseq,env) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  initoutfileinfo(&outfileinfo);
  if (so->outlcptab)
  {
    outfileinfo.lvi = newLcpvalueiterator(encseq,so->readmode,env);
  } else
  {
    outfileinfo.lvi = NULL;
  }
  INITOUTFILEPTR(outfileinfo.outfpsuftab,so->outsuftab,SUFTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfplcptab,so->outlcptab,LCPTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfpllvtab,so->outlcptab,LARGELCPTABSUFFIX);
  INITOUTFILEPTR(outfileinfo.outfpbwttab,so->outbwttab,BWTTABSUFFIX);
  if (!haserr)
  {
    printf("# specialcharacters=" FormatSeqpos "\n",
           PRINTSeqposcast(specialcharinfo.specialcharacters));
    printf("# specialranges=" FormatSeqpos "\n",
           PRINTSeqposcast(specialcharinfo.specialranges));
    if (!so->isplain)
    {
      showcharacterdistribution(alpha,characterdistribution);
    }
    if (so->readmode == Complementmode || so->readmode == Reversecomplementmode)
    {
      if (!isdnaalphabet(alpha,env))
      {
        env_error_set(env,"option %s only can be used for DNA alphabets",
                          so->readmode == Complementmode ? "-cpl" : "rcl");
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab)
    {
      if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
      {
        so->prefixlength = recommendedprefixlength(numofchars,totallength);
        printf("# automatically determined prefixlength = %u\n",
                so->prefixlength);
      } else
      {
        unsigned int maxprefixlen;

        maxprefixlen
          = whatisthemaximalprefixlength(numofchars,
                                         totallength,
                                         (unsigned int) PREFIXLENBITS);
        if (checkprefixlength(maxprefixlen,so->prefixlength,env) != 0)
        {
          haserr = true;
        } else
        {
          showmaximalprefixlength(maxprefixlen,
                                  recommendedprefixlength(
                                  numofchars,
                                  totallength));
        }
      }
      if (!haserr)
      {
        if (suffixeratorwithoutput(
                         &outfileinfo,
                         specialcharinfo.specialcharacters,
                         specialcharinfo.specialranges,
                         encseq,
                         so->readmode,
                         numofchars,
                         so->prefixlength,
                         so->numofparts,
                         mtime,
                         verboseinfo,
                         env) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      if (so->readmode != Forwardmode)
      {
        env_error_set(env,"option '-dir %s' only makes sense in combination "
                          "with at least one of the options -suf, -lcp, or "
                          "-bwt",
                          showreadmode(so->readmode));
        haserr = true;
      }
    }
  }
  env_fa_fclose(outfileinfo.outfpsuftab,env);
  env_fa_fclose(outfileinfo.outfplcptab,env);
  env_fa_fclose(outfileinfo.outfpllvtab,env);
  env_fa_fclose(outfileinfo.outfpbwttab,env);
  if (outfileinfo.lvi != NULL)
  {
    freeLcpvalueiterator(&outfileinfo.lvi,env);
  }
  if (!haserr)
  {
    if (outprjfile(so->str_indexname,
                   so->filenametab,
                   so->readmode,
                   filelengthtab,
                   totallength,
                   numofsequences,
                   &specialcharinfo,
                   so->prefixlength,
                   outfileinfo.numoflargelcpvalues,
                   outfileinfo.maxbranchdepth,
                   &outfileinfo.longest,
                   env) != 0)
    {
      haserr = true;
    }
  }
  FREESPACE(filelengthtab);
  if (alpha != NULL)
  {
    freeAlphabet(&alpha,env);
  }
  freeEncodedsequence(&encseq,env);
  FREESPACE(characterdistribution);
  deliverthetime(stdout,mtime,NULL,env);
  return haserr ? -1 : 0;
}

int parseargsandcallsuffixerator(int argc,const char **argv,Env *env)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  env_error_check(env);
  retval = suffixeratoroptions(&so,argc,argv,env);
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false,env);
    showverbose(verboseinfo,"# sizeof (Seqpos)=%lu\n",
                (unsigned long) (sizeof (Seqpos) * CHAR_BIT));
#ifdef INLINEDENCSEQ
    showverbose(verboseinfo,"# inlined encodeded sequence\n");
#endif
    freeverboseinfo(&verboseinfo,env);
  }
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(false,env);
    if (runsuffixerator(&so,verboseinfo,env) < 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo,env);
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  wrapsfxoptions(&so,env);
  return haserr ? -1 : 0;
}
