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
#include "libgtcore/fa.h"
#include "libgtcore/filelengthvalues.h"
#include "libgtcore/seqiterator.h"
#include "spacedef.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esafileend.h"
#include "readmode-def.h"
#include "verbose-def.h"
#include "intcode-def.h"
#include "sfx-suffixer.h"
#include "sfx-lcpval.h"
#include "sfx-outlcp.h"
#include "stamp.h"

#include "opensfxfile.pr"
#include "fillsci.pr"
#include "sfx-opt.pr"
#include "sfx-outprj.pr"
#include "sfx-apfxlen.pr"

#include "eis-encidxseq.h"
#include "eis-bwtseqconstruct.h"
#include "eis-suffixerator-interface.h"
#include "eis-bwtconstruct_params.h"

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)\
        if (!haserr && (FLAG))\
        {\
          PTR = opensfxfile(so->str_indexname,SUFFIX,"wb",err);\
          if ((PTR) == NULL)\
          {\
            haserr = true;\
          }\
        }

typedef struct
{
  FILE *outfpsuftab,
       *outfpbwttab;
  Seqpos pageoffset;
  const Encodedsequence *encseq;
  DefinedSeqpos longest;
  Outlcpinfo *outlcpinfo;
  Lcpvalueiterator *lvi;
} Outfileinfo;

static int initoutfileinfo(Outfileinfo *outfileinfo,
                           const Encodedsequence *encseq,
                           const Suffixeratoroptions *so,
                           Error *err)
{
  bool haserr = false;

  outfileinfo->outfpsuftab = NULL;
  outfileinfo->outfpbwttab = NULL;
  outfileinfo->pageoffset = 0;
  outfileinfo->longest.defined = false;
  outfileinfo->longest.valueseqpos = 0;
  if (so->outlcptab)
  {
    outfileinfo->lvi = newLcpvalueiterator(encseq,so->readmode);
    outfileinfo->outlcpinfo
      = newlcpoutfileinfo(so->outlcptab ? so->str_indexname : NULL,err,true);
    if (outfileinfo->outlcpinfo == NULL)
    {
      haserr = true;
    }
  } else
  {
    outfileinfo->lvi = NULL;
    outfileinfo->outlcpinfo = NULL;
  }
  INITOUTFILEPTR(outfileinfo->outfpsuftab,so->outsuftab,SUFTABSUFFIX);
  INITOUTFILEPTR(outfileinfo->outfpbwttab,so->outbwttab,BWTTABSUFFIX);
  if (so->outsuftab || so->outbwttab || so->outlcptab)
  {
    outfileinfo->encseq = encseq;
  } else
  {
    outfileinfo->encseq = NULL;
  }
  return haserr  ? -1 : 0;
}

static int suftab2file(Outfileinfo *outfileinfo,
                       const Seqpos *suftab,
                       Readmode readmode,
                       Seqpos numberofsuffixes,
                       Error *err)
{
  bool haserr = false;

  error_check(err);
  if (outfileinfo->outfpsuftab != NULL)
  {
    if (fwrite(suftab,
              sizeof (*suftab),
              (size_t) numberofsuffixes,
              outfileinfo->outfpsuftab)
              != (size_t) numberofsuffixes)
    {
      error_set(err,"cannot write " FormatSeqpos " items of size %u: "
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
          error_set(err,"longest = " FormatSeqpos " is already defined",
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
          cc = getencodedchar(outfileinfo->encseq, /* Random access */
                              startpos - 1,readmode);
        }
      }
      if (outfileinfo->outfpbwttab != NULL)
      {
        if (fwrite(&cc,sizeof (Uchar),(size_t) 1,outfileinfo->outfpbwttab)
                    != (size_t) 1)
        {
          error_set(err,"cannot write 1 item of size %lu: "
                            "errormsg=\"%s\"",
                          (unsigned long) sizeof (Uchar),
                          strerror(errno));
          haserr = true;
          break;
        }
      }
      if (outfileinfo->outlcpinfo != NULL)
      {
        lcpvalue = nextLcpvalueiterator(outfileinfo->lvi,
                                        (outfileinfo->pageoffset == 0)
                                          ? true : false,
                                        suftab,
                                        numberofsuffixes);
        outlcpvalue(lcpvalue,outfileinfo->pageoffset+pos,
                    outfileinfo->outlcpinfo);
      }
    }
  }
  outfileinfo->pageoffset += numberofsuffixes;
  return haserr ? -1 : 0;
}

static int outal1file(const Str *indexname,const Alphabet *alpha,Error *err)
{
  FILE *al1fp;
  bool haserr = false;

  error_check(err);
  al1fp = opensfxfile(indexname,ALPHABETFILESUFFIX,"wb",err);
  if (al1fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    outputalphabet(al1fp,alpha);
    fa_xfclose(al1fp);
  }
  return haserr ? -1 : 0;
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
                 const Str *indexname,
                 Measuretime *mtime,
                 Verboseinfo *verboseinfo,
                 Error *err)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(specialcharacters,
                       specialranges,
                       encseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       indexname,
                       mtime,
                       verboseinfo,
                       err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    while (true)
    {
      suftabptr = nextSfxiterator(&numberofsuffixes,&specialsuffixes,
                                  mtime,sfi,err);
      if (suftabptr == NULL)
      {
        break;
      }
      if (suftab2file(outfileinfo,suftabptr,readmode,numberofsuffixes,err) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi);
  }
  return haserr ? -1 : 0;
}

static unsigned long *initcharacterdistribution(const Alphabet *alpha,
                                                Error *err)
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
                   const unsigned long *characterdistribution,
                   Verboseinfo *verboseinfo)
{
  unsigned int mapsize, idx;

  mapsize = getmapsizeAlphabet(alpha);
  assert(characterdistribution != NULL);
  for (idx=0; idx<mapsize-1; idx++)
  {
    showverbose(verboseinfo,"occurrences(%c)=%lu",
                (int) getprettysymbol(alpha,idx),
                characterdistribution[idx]);
  }
}

static int runsuffixerator(bool doesa,
                           Suffixeratoroptions *so,
                           Verboseinfo *verboseinfo,
                           Error *err)
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

  error_check(err);
  inittheclock(&mtime,
               "determining sequence length and number of special symbols");
  alpha = assigninputalphabet(so->isdna,
                              so->isprotein,
                              so->str_smap,
                              so->filenametab,
                              err);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (!so->isplain)
    {
      characterdistribution = initcharacterdistribution(alpha,err);
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
                                err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    numofchars = getnumofcharsAlphabet(alpha);
    if (outal1file(so->str_indexname,alpha,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    deliverthetime(stdout,mtime,"computing sequence encoding");
    encseq = files2encodedsequence(true,
                                   so->filenametab,
                                   so->isplain,
                                   totallength,
                                   specialcharinfo.specialranges,
                                   alpha,
                                   str_length(so->str_sat) > 0
                                         ? str_get(so->str_sat)
                                         : NULL,
                                   verboseinfo,
                                   err);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      if (so->outtistab)
      {
        if (flushencseqfile(so->str_indexname,encseq,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (!haserr)
  {
    showverbose(verboseinfo,"specialcharacters=" FormatSeqpos,
                PRINTSeqposcast(specialcharinfo.specialcharacters));
    showverbose(verboseinfo,"specialranges=" FormatSeqpos,
                PRINTSeqposcast(specialcharinfo.specialranges));
    if (!so->isplain)
    {
      showcharacterdistribution(alpha,characterdistribution,verboseinfo);
    }
    if (so->readmode == Complementmode || so->readmode == Reversecomplementmode)
    {
      if (!isdnaalphabet(alpha))
      {
        error_set(err,"option %s only can be used for DNA alphabets",
                          so->readmode == Complementmode ? "-cpl" : "rcl");
        haserr = true;
      }
    }
  }
  if (initoutfileinfo(&outfileinfo,encseq,so,err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab || !doesa)
    {
      if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
      {
        so->prefixlength = recommendedprefixlength(numofchars,totallength);
        showverbose(verboseinfo,
                    "automatically determined prefixlength = %u",
                    so->prefixlength);
      } else
      {
        unsigned int maxprefixlen;

        maxprefixlen
          = whatisthemaximalprefixlength(numofchars,
                                         totallength,
                                         (unsigned int) PREFIXLENBITS);
        if (checkprefixlength(maxprefixlen,so->prefixlength,err) != 0)
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
        if (doesa)
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
                             so->outlcptab ? so->str_indexname : NULL,
                             mtime,
                             verboseinfo,
                             err) != 0)
          {
            haserr = true;
          }
        } else
        {
          sfxInterface *si;
          BWTSeq *bwtSeq;
          showverbose(verboseinfo, "run construction of packed index for:\n"
                      "blocksize=%u\nblocks-per-bucket=%u\nlocfreq=%u",
                      so->bwtIdxParams.final.seqParams.blockEnc.blockSize,
                      so->bwtIdxParams.final.seqParams.blockEnc.bucketBlocks,
                      so->bwtIdxParams.final.locateInterval);
          if (!(si = newSfxInterface(so, encseq, &specialcharinfo,
                                     numofsequences, mtime, totallength + 1,
                                     alpha, characterdistribution,
                                     verboseinfo, err)))
          {
            fputs("Index creation failed.\n", stderr);
            haserr = true;
          }
          else if (
            !(bwtSeq = createBWTSeqFromSfxI(&so->bwtIdxParams.final, si,
                                            getSfxILength(si), err)))
          {
            fputs("Index creation failed.\n", stderr);
            deleteSfxInterface(si);
            haserr = true;
          }
          else
            deleteBWTSeq(bwtSeq); /**< the actual object is not * used here */
          outfileinfo.longest = getSfxILongestPos(si); /* XXX Thomas, is
                                                          this necessary */

          deleteSfxInterface(si);
        }
      }
    } else
    {
      if (so->readmode != Forwardmode)
      {
        error_set(err,"option '-dir %s' only makes sense in combination "
                          "with at least one of the options -suf, -lcp, or "
                          "-bwt",
                          showreadmode(so->readmode));
        haserr = true;
      }
    }
  }
  fa_fclose(outfileinfo.outfpsuftab);
  fa_fclose(outfileinfo.outfpbwttab);
  if (outfileinfo.lvi != NULL)
  {
    freeLcpvalueiterator(&outfileinfo.lvi);
  }
  if (!haserr)
  {
    Seqpos numoflargelcpvalues,
           maxbranchdepth;

    if (outfileinfo.outlcpinfo == NULL)
    {
      numoflargelcpvalues = maxbranchdepth = 0;
    } else
    {
      numoflargelcpvalues = getnumoflargelcpvalues(outfileinfo.outlcpinfo);
      maxbranchdepth = getmaxbranchdepth(outfileinfo.outlcpinfo);
    }
    if (outprjfile(so->str_indexname,
                   so->filenametab,
                   so->readmode,
                   filelengthtab,
                   totallength,
                   numofsequences,
                   &specialcharinfo,
                   so->prefixlength,
                   numoflargelcpvalues,
                   maxbranchdepth,
                   &outfileinfo.longest,
                   err) != 0)
    {
      haserr = true;
    }
  }
  if (outfileinfo.outlcpinfo != NULL)
  {
    freeoutlcptab(&outfileinfo.outlcpinfo);
  }
  FREESPACE(filelengthtab);
  if (alpha != NULL)
  {
    freeAlphabet(&alpha);
  }
  freeEncodedsequence(&encseq);
  FREESPACE(characterdistribution);
  deliverthetime(stdout,mtime,NULL);
  return haserr ? -1 : 0;
}

int parseargsandcallsuffixerator(bool doesa,int argc,
                                 const char **argv,Error *err)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  error_check(err);
  retval = suffixeratoroptions(&so,doesa,argc,argv,err);
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(so.beverbose);

    showverbose(verboseinfo,"sizeof (Seqpos)=%lu",
                (unsigned long) (sizeof (Seqpos) * CHAR_BIT));
#ifdef INLINEDENCSEQ
    showverbose(verboseinfo,"inlined encodeded sequence");
#endif
    freeverboseinfo(&verboseinfo);
  }
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(so.beverbose);
    if (runsuffixerator(doesa,&so,verboseinfo,err) < 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  wrapsfxoptions(&so,err);
  return haserr ? -1 : 0;
}
