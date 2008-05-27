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
#include "libgtcore/unused.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esafileend.h"
#include "readmode-def.h"
#include "verbose-def.h"
#include "intcode-def.h"
#include "stamp.h"
#include "sfx-suffixer.h"
#include "sfx-outlcp.h"
#include "sfx-input.h"
#include "sfx-run.h"

#include "opensfxfile.pr"
#include "sfx-opt.pr"
#include "sfx-outprj.pr"
#include "sfx-apfxlen.pr"

#include "eis-encidxseq.h"
#include "eis-bwtseq-construct.h"
#include "eis-bwtseq-param.h"
#include "eis-suffixerator-interface.h"

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
       *outfpbwttab,
       *outfpbcktab;
  Seqpos pageoffset;
  const Encodedsequence *encseq;
  DefinedSeqpos longest;
  Outlcpinfo *outlcpinfo;
} Outfileinfo;

static int initoutfileinfo(Outfileinfo *outfileinfo,
                           unsigned int prefixlength,
                           unsigned int numofchars,
                           const Encodedsequence *encseq,
                           const Suffixeratoroptions *so,
                           Error *err)
{
  bool haserr = false;

  outfileinfo->outfpsuftab = NULL;
  outfileinfo->outfpbwttab = NULL;
  outfileinfo->outfpbcktab = NULL;
  outfileinfo->pageoffset = 0;
  outfileinfo->longest.defined = false;
  outfileinfo->longest.valueseqpos = 0;
  if (so->outlcptab)
  {
    outfileinfo->outlcpinfo
      = newlcpoutinfo(so->outlcptab ? so->str_indexname : NULL,
                      prefixlength,
                      numofchars,
                      getencseqtotallength(encseq),err);
    if (outfileinfo->outlcpinfo == NULL)
    {
      haserr = true;
    }
  } else
  {
    outfileinfo->outlcpinfo = NULL;
  }
  INITOUTFILEPTR(outfileinfo->outfpsuftab,so->outsuftab,SUFTABSUFFIX);
  INITOUTFILEPTR(outfileinfo->outfpbwttab,so->outbwttab,BWTTABSUFFIX);
  INITOUTFILEPTR(outfileinfo->outfpbcktab,so->outbcktab,BCKTABSUFFIX);
  if (so->outsuftab || so->outbwttab || so->outlcptab || so->outbcktab)
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
    Seqpos startpos, pos;
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
    }
  }
  outfileinfo->pageoffset += numberofsuffixes;
  return haserr ? -1 : 0;
}

static int suffixeratorwithoutput(
                 Outfileinfo *outfileinfo,
                 Seqpos specialcharacters,
                 Seqpos realspecialranges,
                 const Encodedsequence *encseq,
                 Readmode readmode,
                 unsigned int numofchars,
                 const Uchar *characters,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 const Sfxstrategy *sfxstrategy,
                 Measuretime *mtime,
                 Verboseinfo *verboseinfo,
                 Error *err)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(specialcharacters,
                       realspecialranges,
                       encseq,
                       readmode,
                       numofchars,
                       characters,
                       prefixlength,
                       numofparts,
                       outfileinfo->outlcpinfo,
                       sfxstrategy,
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
                                  mtime,sfi);
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
  if (outfileinfo->outfpbcktab != NULL)
  {
    if (sfibcktab2file(outfileinfo->outfpbcktab,sfi,err) != 0)
    {
      haserr = true;
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi);
  }
  return haserr ? -1 : 0;
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

static void showsequencefeatures(Verboseinfo *verboseinfo,
                                 const Specialcharinfo *specialcharinfo,
                                 const Alphabet *alpha,
                                 const unsigned long *characterdistribution)
{
  showverbose(verboseinfo,"specialcharacters=" FormatSeqpos,
              PRINTSeqposcast(specialcharinfo->specialcharacters));
  showverbose(verboseinfo,"specialranges=" FormatSeqpos,
              PRINTSeqposcast(specialcharinfo->specialranges));
  showverbose(verboseinfo,"realspecialranges=" FormatSeqpos,
              PRINTSeqposcast(specialcharinfo->realspecialranges));
  if (characterdistribution != NULL)
  {
    showcharacterdistribution(alpha,characterdistribution,verboseinfo);
  }
}

static int detpfxlenandmaxdepth(unsigned int *prefixlength,
                                Definedunsignedint *maxdepth,
                                const Suffixeratoroptions *so,
                                unsigned int numofchars,
                                Seqpos totallength,
                                Verboseinfo *verboseinfo,
                                Error *err)
{
  bool haserr = false;

  if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
  {
    *prefixlength = recommendedprefixlength(numofchars,totallength);
    showverbose(verboseinfo,
                "automatically determined prefixlength = %u",
                *prefixlength);
  } else
  {
    unsigned int maxprefixlen;

    *prefixlength = so->prefixlength;
    maxprefixlen
      = whatisthemaximalprefixlength(numofchars,
                                     totallength,
                                     so->sfxstrategy.storespecialcodes
                                     ? getprefixlenbits()
                                     : 0);
    if (checkprefixlength(maxprefixlen,*prefixlength,err) != 0)
    {
      haserr = true;
    } else
    {
      showmaximalprefixlength(verboseinfo,
                              maxprefixlen,
                              recommendedprefixlength(
                              numofchars,
                              totallength));
    }
  }
  if (!haserr && so->maxdepth.defined)
  {
    if (so->maxdepth.valueunsignedint == MAXDEPTH_AUTOMATIC)
    {
      maxdepth->defined = true;
      maxdepth->valueunsignedint = *prefixlength;
      showverbose(verboseinfo,
                  "automatically determined maxdepth = %u",
                  maxdepth->valueunsignedint);
    } else
    {
      if (so->maxdepth.valueunsignedint < *prefixlength)
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint = *prefixlength;
        showverbose(verboseinfo,
                    "set maxdepth = %u",maxdepth->valueunsignedint);
      } else
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint = so->maxdepth.valueunsignedint;
        showverbose(verboseinfo,
                    "use maxdepth = %u",maxdepth->valueunsignedint);
      }
    }
  }
  return haserr ? -1 : 0;
}

static int run_packedindexconstruction(Verboseinfo *verboseinfo,
                                       Measuretime *mtime,
                                       FILE *outfpbcktab,
                                       const Suffixeratoroptions *so,
                                       unsigned int prefixlength,
                                       Sfxseqinfo *sfxseqinfo,
                                       const Sfxstrategy *sfxstrategy,
                                       Error *err)
{
  sfxInterface *si;
  BWTSeq *bwtSeq;
  const Sfxiterator *sfi;
  bool haserr = false;

  showverbose(verboseinfo, "run construction of packed index for:\n"
              "blocksize=%u\nblocks-per-bucket=%u\nlocfreq=%u",
              so->bwtIdxParams.final.seqParams.encParams.blockEnc.blockSize,
              so->bwtIdxParams.final.seqParams.encParams.blockEnc.bucketBlocks,
              so->bwtIdxParams.final.locateInterval);
  si = newSfxInterface(so->readmode,
                       prefixlength,
                       so->numofparts,
                       sfxstrategy,
                       sfxseqinfo->encseq,
                       &sfxseqinfo->specialcharinfo,
                       sfxseqinfo->numofsequences,
                       mtime,
                       getencseqtotallength(sfxseqinfo->encseq) + 1,
                       sfxseqinfo->alpha,
                       sfxseqinfo->characterdistribution,
                       verboseinfo, err);
  if (si == NULL)
  {
    haserr = true;
  } else
  {
    bwtSeq = createBWTSeqFromSfxI(&so->bwtIdxParams.final, si, err);
    if (bwtSeq == NULL)
    {
      deleteSfxInterface(si);
      haserr = true;
    } else
    {
      deleteBWTSeq(bwtSeq); /**< the actual object is not * used here */
      /*
        outfileinfo.longest = SfxIGetRot0Pos(si);
      */
      sfi = SfxInterface2Sfxiterator(si);
      assert(sfi != NULL);
      if (outfpbcktab != NULL)
      {
        if (sfibcktab2file(outfpbcktab,sfi,err) != 0)
        {
          haserr = true;
        }
      }
      deleteSfxInterface(si);
    }
  }
  return haserr ? -1 : 0;
}

static int runsuffixerator(bool doesa,
                           const Suffixeratoroptions *so,
                           Verboseinfo *verboseinfo,
                           Error *err)
{
  Measuretime *mtime = NULL;
  Outfileinfo outfileinfo;
  bool haserr = false;
  Sfxseqinfo sfxseqinfo;
  unsigned int prefixlength;
  Seqpos totallength = 0;
  Sfxstrategy sfxstrategy;

  error_check(err);
  if (so->showtime)
  {
    mtime = inittheclock("determining sequence length and number of "
                         "special symbols");
  }
  if (str_length(so->str_inputindex) > 0)
  {
    if (fromsarr2Sfxseqinfo(&sfxseqinfo,
                            so->str_inputindex,
                            verboseinfo,
                            err) != 0)
    {
      haserr = true;
    }
    if (so->outtistab && strcmp(str_get(so->str_inputindex),
                                str_get(so->str_indexname)) != 0)
    {
      if (flushencseqfile(so->str_indexname,sfxseqinfo.encseq,err) != 0)
      {
        haserr = true;
      }
    }
  } else
  {
    if (fromfiles2Sfxseqinfo(&sfxseqinfo,
                             mtime,
                             so,
                             verboseinfo,
                             err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    totallength = getencseqtotallength(sfxseqinfo.encseq);
  }
  if (!haserr)
  {
    showsequencefeatures(verboseinfo,
                         &sfxseqinfo.specialcharinfo,
                         sfxseqinfo.alpha,
                         sfxseqinfo.characterdistribution);
    if (sfxseqinfo.characterdistribution != NULL)
    {
      if (so->readmode == Complementmode ||
          so->readmode == Reversecomplementmode)
      {
        if (!isdnaalphabet(sfxseqinfo.alpha))
        {
          error_set(err,"option -%s only can be used for DNA alphabets",
                            so->readmode == Complementmode ? "cpl" : "rcl");
          haserr = true;
        }
      }
    }
  }
  prefixlength = so->prefixlength;
  sfxstrategy = so->sfxstrategy;
  sfxstrategy.maxdepth.defined = false;
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab || so->outbcktab ||
        !doesa)
    {
      unsigned int numofchars = getnumofcharsAlphabet(sfxseqinfo.alpha);

      if (detpfxlenandmaxdepth(&prefixlength,
                               &sfxstrategy.maxdepth,
                               so,
                               numofchars,
                               totallength,
                               verboseinfo,
                               err) != 0)
      {
        haserr = true;
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
  outfileinfo.outfpsuftab = NULL;
  outfileinfo.outfpbwttab = NULL;
  outfileinfo.outlcpinfo = NULL;
  outfileinfo.outfpbcktab = NULL;
  if (!haserr)
  {
    if (initoutfileinfo(&outfileinfo,prefixlength,
                        getnumofcharsAlphabet(sfxseqinfo.alpha),
                        sfxseqinfo.encseq,so,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab || !doesa)
    {
      if (doesa)
      {
        if (suffixeratorwithoutput(
                           &outfileinfo,
                           sfxseqinfo.specialcharinfo.specialcharacters,
                           sfxseqinfo.specialcharinfo.realspecialranges,
                           sfxseqinfo.encseq,
                           so->readmode,
                           getnumofcharsAlphabet(sfxseqinfo.alpha),
                           getcharactersAlphabet(sfxseqinfo.alpha),
                           prefixlength,
                           so->numofparts,
                           &sfxstrategy,
                           mtime,
                           verboseinfo,
                           err) != 0)
        {
          haserr = true;
        }
      } else
      {
        if (run_packedindexconstruction(verboseinfo,
                                        mtime,
                                        outfileinfo.outfpbcktab,
                                        so,
                                        prefixlength,
                                        &sfxseqinfo,
                                        &sfxstrategy,
                                        err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  fa_fclose(outfileinfo.outfpsuftab);
  fa_fclose(outfileinfo.outfpbwttab);
  fa_fclose(outfileinfo.outfpbcktab);
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
    assert(sfxseqinfo.numofsequences > 0);
    assert(sfxseqinfo.filelengthtab != NULL);
    if (outprjfile(so->str_indexname,
                   sfxseqinfo.filenametab,
                   sfxseqinfo.readmode,
                   sfxseqinfo.filelengthtab,
                   totallength,
                   sfxseqinfo.numofsequences,
                   &sfxseqinfo.specialcharinfo,
                   prefixlength,
                   &sfxstrategy.maxdepth,
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
  freeSfxseqinfo(&sfxseqinfo,
                 (str_length(so->str_inputindex) > 0) ? true : false);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,NULL);
  }
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
  wrapsfxoptions(&so);
  return haserr ? -1 : 0;
}
