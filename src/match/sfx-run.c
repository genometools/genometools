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
#include "core/chardef.h"
#include "core/fa.h"
#include "core/unused_api.h"
#include "alphadef.h"
#include "sfx-optdef.h"
#include "encseq-def.h"
#include "measure-time-if.h"
#include "esa-fileend.h"
#include "readmode-def.h"
#include "verbose-def.h"
#include "intcode-def.h"
#include "spacedef.h"
#include "stamp.h"
#include "sfx-suffixer.h"
#include "sfx-bentsedg.h"
#include "sfx-input.h"
#include "sfx-run.h"
#include "opensfxfile.h"
#include "stamp.h"

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
                           const Encodedsequence *encseq,
                           const Suffixeratoroptions *so,
                           GtError *err)
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
      = newOutlcpinfo(so->outlcptab ? so->str_indexname : NULL,
                      prefixlength,
                      getencseqAlphabetnumofchars(encseq),
                      getencseqtotallength(encseq),
                      so->sfxstrategy.ssortmaxdepth.defined ? false : true,
                      err);
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
                       GtError *err)
{
  bool haserr = false;

  gt_error_check(err);
  if (outfileinfo->outfpsuftab != NULL)
  {
    if (fwrite(suftab,
              sizeof (*suftab),
              (size_t) numberofsuffixes,
              outfileinfo->outfpsuftab)
              != (size_t) numberofsuffixes)
    {
      gt_error_set(err,"cannot write " FormatSeqpos " items of size %u: "
                    "errormsg=\"%s\"",
           PRINTSeqposcast(numberofsuffixes),
           (unsigned int) sizeof (*suftab),
           strerror(errno));
      haserr = true;
    }
  }
  printf("longest.defined = %s\n",
         outfileinfo->longest.defined ? "true" : "false");
  if (!haserr &&
      (!outfileinfo->longest.defined || outfileinfo->outfpbwttab != NULL))
  {
    Seqpos startpos, pos;
    Uchar cc = 0;

    STAMP;
    for (pos=0; pos < numberofsuffixes; pos++)
    {
      startpos = suftab[pos];
      if (startpos == 0)
      {
        cc = (Uchar) UNDEFBWTCHAR;
        if (!outfileinfo->longest.defined)
        {
          STAMP;
          outfileinfo->longest.defined = true;
          outfileinfo->longest.valueseqpos = outfileinfo->pageoffset + pos;
          printf("set longest = %lu\n",(unsigned long) 
                     outfileinfo->longest.valueseqpos);
        }
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
          gt_error_set(err,"cannot write 1 item of size %lu: "
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
                 const Encodedsequence *encseq,
                 Readmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 const Sfxstrategy *sfxstrategy,
                 Measuretime *mtime,
                 Verboseinfo *verboseinfo,
                 GtError *err)
{
  const Seqpos *suftabptr;
  Seqpos numberofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(encseq,
                       readmode,
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
      Seqpos longest;

      suftabptr = nextSfxiterator(&numberofsuffixes,&specialsuffixes,
                                  mtime,sfi);
      if (suftabptr == NULL)
      {
        break;
      }
      if (numofparts == 1U && sfi2longestsuffixpos(&longest,sfi))
      {
        outfileinfo->longest.defined = true;
        outfileinfo->longest.valueseqpos = longest;
        printf("set longest = %lu\n",(unsigned long) 
                     outfileinfo->longest.valueseqpos);

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
                   const SfxAlphabet *alpha,
                   const unsigned long *characterdistribution,
                   Verboseinfo *verboseinfo)
{
  unsigned int numofchars, idx;

  numofchars = getnumofcharsAlphabet(alpha);
  gt_assert(characterdistribution != NULL);
  for (idx=0; idx<numofchars; idx++)
  {
    showverbose(verboseinfo,"occurrences(%c)=%lu",
                (int) getprettysymbol(alpha,idx),
                characterdistribution[idx]);
  }
}

static void showsequencefeatures(Verboseinfo *verboseinfo,
                                 const Encodedsequence *encseq,
                                 const unsigned long *characterdistribution)
{
  showverbose(verboseinfo,"specialcharacters=" FormatSeqpos,
              PRINTSeqposcast(getencseqspecialcharacters(encseq)));
  showverbose(verboseinfo,"specialranges=" FormatSeqpos,
              PRINTSeqposcast(getencseqspecialranges(encseq)));
  showverbose(verboseinfo,"realspecialranges=" FormatSeqpos,
              PRINTSeqposcast(getencseqrealspecialranges(encseq)));
  if (characterdistribution != NULL)
  {
    const SfxAlphabet *alpha = getencseqAlphabet(encseq);
    showcharacterdistribution(alpha,characterdistribution,verboseinfo);
  }
}

static int detpfxlenandmaxdepth(unsigned int *prefixlength,
                                Definedunsignedint *maxdepth,
                                const Suffixeratoroptions *so,
                                unsigned int numofchars,
                                Seqpos totallength,
                                Verboseinfo *verboseinfo,
                                GtError *err)
{
  bool haserr = false;

  if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
  {
    *prefixlength = recommendedprefixlength(numofchars,totallength);
    showverbose(verboseinfo,
                "automatically determined prefixlength=%u",
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
  if (!haserr && so->sfxstrategy.ssortmaxdepth.defined)
  {
    if (so->sfxstrategy.ssortmaxdepth.valueunsignedint == MAXDEPTH_AUTOMATIC)
    {
      maxdepth->defined = true;
      maxdepth->valueunsignedint = *prefixlength;
      showverbose(verboseinfo,
                  "automatically determined maxdepth=%u",
                  maxdepth->valueunsignedint);
    } else
    {
      if (so->sfxstrategy.ssortmaxdepth.valueunsignedint < *prefixlength)
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint = *prefixlength;
        showverbose(verboseinfo,
                    "set maxdepth=%u",maxdepth->valueunsignedint);
      } else
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint
          = so->sfxstrategy.ssortmaxdepth.valueunsignedint;
        showverbose(verboseinfo,
                    "use maxdepth=%u",maxdepth->valueunsignedint);
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
                                       GtError *err)
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
                       mtime,
                       getencseqtotallength(sfxseqinfo->encseq) + 1,
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
      gt_assert(sfi != NULL);
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
                           GtError *err)
{
  Measuretime *mtime = NULL;
  Outfileinfo outfileinfo;
  bool haserr = false;
  Sfxseqinfo sfxseqinfo;
  unsigned int prefixlength;
  Sfxstrategy sfxstrategy;

  gt_error_check(err);
  if (so->showtime)
  {
    mtime = inittheclock("determining sequence length and number of "
                         "special symbols");
  }
  if (gt_str_length(so->str_inputindex) > 0)
  {
    if (fromsarr2Sfxseqinfo(&sfxseqinfo,
                            so->str_inputindex,
                            so->readmode,
                            verboseinfo,
                            err) != 0)
    {
      haserr = true;
    }
    if (so->outtistab && strcmp(gt_str_get(so->str_inputindex),
                                gt_str_get(so->str_indexname)) != 0)
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
    if (!haserr && so->outssptab)
    {
      FILE *outfp;

      outfp = opensfxfile(so->str_indexname,SSPTABSUFFIX,"wb",err);
      if (outfp == NULL)
      {
        haserr = true;
      } else
      {
        if (fwrite(sfxseqinfo.sequenceseppos.spaceSeqpos,
                   sizeof (*sfxseqinfo.sequenceseppos.spaceSeqpos),
                   (size_t) sfxseqinfo.sequenceseppos.nextfreeSeqpos,
                   outfp)
                   != (size_t) sfxseqinfo.sequenceseppos.nextfreeSeqpos)
        {
          gt_error_set(err,"cannot write %lu items of size %u: "
                           "errormsg=\"%s\"",
                            sfxseqinfo.sequenceseppos.nextfreeSeqpos,
                            (unsigned int) sizeof (*sfxseqinfo.sequenceseppos.
                                                   spaceSeqpos),
                            strerror(errno));
          haserr = true;
        }
      }
      FREESPACE(sfxseqinfo.sequenceseppos.spaceSeqpos);
      gt_fa_fclose(outfp);
    }
  }
  if (!haserr)
  {
    showsequencefeatures(verboseinfo,
                         sfxseqinfo.encseq,
                         sfxseqinfo.characterdistribution);
    if (sfxseqinfo.characterdistribution != NULL)
    {
      if (so->readmode == Complementmode ||
          so->readmode == Reversecomplementmode)
      {
        if (!isdnaalphabet(getencseqAlphabet(sfxseqinfo.encseq)))
        {
          gt_error_set(err,"option -%s only can be used for DNA alphabets",
                            so->readmode == Complementmode ? "cpl" : "rcl");
          haserr = true;
        }
      }
    }
  }
  prefixlength = so->prefixlength;
  sfxstrategy = so->sfxstrategy;
  sfxstrategy.ssortmaxdepth.defined = false;
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab || so->outbcktab ||
        !doesa)
    {
      unsigned int numofchars = getencseqAlphabetnumofchars(sfxseqinfo.encseq);

      if (detpfxlenandmaxdepth(&prefixlength,
                               &sfxstrategy.ssortmaxdepth,
                               so,
                               numofchars,
                               getencseqtotallength(sfxseqinfo.encseq),
                               verboseinfo,
                               err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (so->readmode != Forwardmode)
      {
        gt_error_set(err,"option '-dir %s' only makes sense in combination "
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
                           sfxseqinfo.encseq,
                           so->readmode,
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
  gt_fa_fclose(outfileinfo.outfpsuftab);
  gt_fa_fclose(outfileinfo.outfpbwttab);
  gt_fa_fclose(outfileinfo.outfpbcktab);
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
                   sfxseqinfo.readmode,
                   sfxseqinfo.encseq,
                   prefixlength,
                   &sfxstrategy.ssortmaxdepth,
                   numoflargelcpvalues,
                   maxbranchdepth,
                   &outfileinfo.longest,
                   err) != 0)
    {
      haserr = true;
    }
  }
  if (gt_str_length(so->str_inputindex) == 0 && sfxseqinfo.encseq != NULL)
  {
    removefilenametabref(sfxseqinfo.encseq);
  }
  if (outfileinfo.outlcpinfo != NULL)
  {
    freeOutlcptab(&outfileinfo.outlcpinfo);
  }
  freeSfxseqinfo(&sfxseqinfo);
  if (mtime != NULL)
  {
    deliverthetime(stdout,mtime,NULL);
  }
  return haserr ? -1 : 0;
}

int parseargsandcallsuffixerator(bool doesa,int argc,
                                 const char **argv,GtError *err)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  gt_error_check(err);
  retval = suffixeratoroptions(&so,doesa,argc,argv,err);
  if (retval == 0)
  {
    Verboseinfo *verboseinfo = newverboseinfo(so.beverbose);

    showverbose(verboseinfo,"sizeof (Seqpos)=%lu",
                (unsigned long) (sizeof (Seqpos) * CHAR_BIT));
    if (runsuffixerator(doesa,&so,verboseinfo,err) < 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
    /*showgetencodedcharcounters(); */
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
