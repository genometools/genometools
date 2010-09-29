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
#include "core/alphabet.h"
#include "core/chardef.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "core/fa.h"
#include "core/logger.h"
#include "core/progress_timer_api.h"
#include "core/readmode.h"
#include "core/showtime.h"
#include "core/unused_api.h"
#include "esa-fileend.h"
#include "giextract.h"
#include "stamp.h"
#include "intcode-def.h"
#include "sfx-optdef.h"
#include "sfx-suffixer.h"
#include "sfx-run.h"
#include "sfx-opt.pr"
#include "sfx-outprj.pr"
#include "sfx-apfxlen.h"
#include "sfx-bentsedg.h"
#include "sfx-suffixgetset.h"

#ifndef S_SPLINT_S
#include "eis-encidxseq.h"
#include "eis-bwtseq-construct.h"
#include "eis-bwtseq-param.h"
#include "eis-suffixerator-interface.h"
#endif

#define INITOUTFILEPTR(PTR,FLAG,SUFFIX)\
        if (!haserr && (FLAG))\
        {\
          PTR = gt_fa_fopen_with_suffix(gt_str_get(so->fn2encopt.indexname),\
                                        SUFFIX,"wb",err);\
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
  unsigned long pageoffset;
  const GtEncseq *encseq;
  Definedunsignedlong longest;
  Outlcpinfo *outlcpinfo;
} Outfileinfo;

static int initoutfileinfo(Outfileinfo *outfileinfo,
                          unsigned int prefixlength,
                          const GtEncseq *encseq,
                          const Suffixeratoroptions *so,
                          GtError *err)
{
  bool haserr = false;

  outfileinfo->outfpsuftab = NULL;
  outfileinfo->outfpbwttab = NULL;
  outfileinfo->outfpbcktab = NULL;
  outfileinfo->pageoffset = 0;
  outfileinfo->longest.defined = false;
  outfileinfo->longest.valueunsignedlong = 0;
  if (so->outlcptab)
  {
    outfileinfo->outlcpinfo
      = gt_newOutlcpinfo(so->outlcptab ?
                         gt_str_get(so->fn2encopt.indexname) : NULL,
                         gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq)),
                         prefixlength,
                         gt_encseq_total_length(encseq),
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

static int bwttab2file(Outfileinfo *outfileinfo,
                       const GtSuffixsortspace *suffixsortspace,
                       GtReadmode readmode,
                       unsigned long numberofsuffixes,
                       GtError *err)
{
  bool haserr = false;

  gt_error_check(err);
  if (outfileinfo->outfpbwttab != NULL)
  {
    unsigned long startpos, pos;
    GtUchar cc = 0;

    for (pos=0; pos < numberofsuffixes; pos++)
    {
      startpos = gt_suffixsortspace_getdirect(suffixsortspace,pos);
      if (startpos == 0)
      {
        cc = (GtUchar) UNDEFBWTCHAR;
      } else
      {
        if (outfileinfo->outfpbwttab != NULL)
        {
          /* Random access */
          cc = gt_encseq_get_encoded_char(outfileinfo->encseq,
                                          startpos - 1,
                                          readmode);
        }
      }
      if (outfileinfo->outfpbwttab != NULL)
      {
        if (fwrite(&cc,sizeof (GtUchar),(size_t) 1,outfileinfo->outfpbwttab)
                    != (size_t) 1)
        {
          gt_error_set(err,"cannot write 1 item of size %u: errormsg=\"%s\"",
                          (unsigned int) sizeof (GtUchar),
                          strerror(errno));
          haserr = true;
          break;
        }
      }
    }
  }
  return haserr ? -1 : 0;
}

static int suffixeratorwithoutput(const GtStr *indexname,
                                  Outfileinfo *outfileinfo,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  unsigned int prefixlength,
                                  unsigned int numofparts,
                                  unsigned long maximumspace,
                                  const Sfxstrategy *sfxstrategy,
                                  GtProgressTimer *sfxprogress,
                                  bool withprogressbar,
                                  GtLogger *logger,
                                  GtError *err)
{
  unsigned long numberofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi = NULL;

  sfi = gt_Sfxiterator_new(encseq,
                           readmode,
                           prefixlength,
                           numofparts,
                           maximumspace,
                           outfileinfo->outlcpinfo,
                           sfxstrategy,
                           sfxprogress,
                           withprogressbar,
                           logger,
                           err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    const GtSuffixsortspace *suffixsortspace;
    while (true)
    {
      suffixsortspace = gt_Sfxiterator_next(&numberofsuffixes,&specialsuffixes,
                                            sfi);
      if (suffixsortspace == NULL)
      {
        break;
      }
      if (outfileinfo->outfpsuftab != NULL)
      {
        if (gt_suffixsortspace_to_file (outfileinfo->outfpsuftab,
                                        suffixsortspace,numberofsuffixes,
                                        err) != 0)
        {
          haserr = true;
          break;
        }
      }
      /* XXX be careful: postpone not output bwtab if streamsuftab is true */
      if (bwttab2file(outfileinfo,suffixsortspace,readmode,numberofsuffixes,
                      err) != 0)
      {
        haserr = true;
        break;
      }
      outfileinfo->pageoffset += numberofsuffixes;
    }
  }
  outfileinfo->longest.defined = true;
  outfileinfo->longest.valueunsignedlong = gt_Sfxiterator_longest(sfi);
  if (!haserr && sfxstrategy->streamsuftab)
  {
    gt_fa_fclose(outfileinfo->outfpsuftab);
    outfileinfo->outfpsuftab = NULL;
    if (gt_Sfxiterator_postsortfromstream(sfi,indexname,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && outfileinfo->outfpbcktab != NULL)
  {
    if (gt_Sfxiterator_bcktab2file(outfileinfo->outfpbcktab,sfi,err) != 0)
    {
      haserr = true;
    }
  }
  gt_Sfxiterator_delete(sfi);
  return haserr ? -1 : 0;
}

static int detpfxlenandmaxdepth(unsigned int *prefixlength,
                                Definedunsignedint *maxdepth,
                                const Suffixeratoroptions *so,
                                unsigned int numofchars,
                                unsigned long totallength,
                                GtLogger *logger,
                                GtError *err)
{
  bool haserr = false;

  if (so->prefixlength == PREFIXLENGTH_AUTOMATIC)
  {
    *prefixlength = gt_recommendedprefixlength(numofchars,totallength);
    gt_logger_log(logger,
                "automatically determined prefixlength=%u",
                *prefixlength);
  } else
  {
    unsigned int maxprefixlen;

    *prefixlength = so->prefixlength;
    maxprefixlen
      = gt_whatisthemaximalprefixlength(numofchars,
                                        totallength,
                                        so->sfxstrategy.storespecialcodes
                                        ? (unsigned int) PREFIXLENBITS
                                        : 0);
    if (gt_checkprefixlength(maxprefixlen,*prefixlength,err) != 0)
    {
      haserr = true;
    } else
    {
      gt_showmaximalprefixlength(logger,
                              maxprefixlen,
                              gt_recommendedprefixlength(
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
      gt_logger_log(logger,
                  "automatically determined maxdepth=%u",
                  maxdepth->valueunsignedint);
    } else
    {
      if (so->sfxstrategy.ssortmaxdepth.valueunsignedint < *prefixlength)
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint = *prefixlength;
        gt_logger_log(logger,
                    "set maxdepth=%u",maxdepth->valueunsignedint);
      } else
      {
        maxdepth->defined = true;
        maxdepth->valueunsignedint
          = so->sfxstrategy.ssortmaxdepth.valueunsignedint;
        gt_logger_log(logger,
                    "use maxdepth=%u",maxdepth->valueunsignedint);
      }
    }
  }
  return haserr ? -1 : 0;
}

#ifndef S_SPLINT_S
static int run_packedindexconstruction(GtLogger *logger,
                                       GtProgressTimer *sfxprogress,
                                       bool withprogressbar,
                                       FILE *outfpbcktab,
                                       const Suffixeratoroptions *so,
                                       unsigned int prefixlength,
                                       const GtEncseq *encseq,
                                       const Sfxstrategy *sfxstrategy,
                                       GtError *err)
{
  sfxInterface *si;
  BWTSeq *bwtSeq;
  const Sfxiterator *sfi;
  bool haserr = false;

  gt_logger_log(logger, "run construction of packed index for:\n"
              "blocksize=%u\nblocks-per-bucket=%u\nlocfreq=%u",
              so->bwtIdxParams.final.seqParams.encParams.blockEnc.blockSize,
              so->bwtIdxParams.final.seqParams.encParams.blockEnc.bucketBlocks,
              so->bwtIdxParams.final.locateInterval);
  si = gt_newSfxInterface(so->readmode,
                          prefixlength,
                          so->numofparts,
                          so->maximumspace,
                          sfxstrategy,
                          encseq,
                          sfxprogress,
                          withprogressbar,
                          gt_encseq_total_length(encseq) + 1,
                          logger,
                          err);
  if (si == NULL)
  {
    haserr = true;
  } else
  {
    bwtSeq = gt_createBWTSeqFromSfxI(&so->bwtIdxParams.final, si, err);
    if (bwtSeq == NULL)
    {
      gt_deleteSfxInterface(si);
      haserr = true;
    } else
    {
      gt_deleteBWTSeq(bwtSeq); /**< the actual object is not * used here */
      /*
        outfileinfo.longest = gt_SfxIGetRot0Pos(si);
      */
      sfi = gt_SfxInterface2Sfxiterator(si);
      gt_assert(sfi != NULL);
      if (outfpbcktab != NULL)
      {
        if (gt_Sfxiterator_bcktab2file(outfpbcktab,sfi,err) != 0)
        {
          haserr = true;
        }
      }
      gt_deleteSfxInterface(si);
    }
  }
  return haserr ? -1 : 0;
}
#else

#endif

static int runsuffixerator(bool doesa,
                           const Suffixeratoroptions *so,
                           GtLogger *logger,
                           GtError *err)
{
  GtProgressTimer *sfxprogress = NULL;
  Outfileinfo outfileinfo;
  bool haserr = false;
  unsigned int prefixlength;
  Sfxstrategy sfxstrategy;
  GtEncseq *encseq = NULL;
  GtEncseqEncoder *ee = NULL;

  gt_error_check(err);

  ee = gt_encseq_encoder_new();
  if (gt_showtime_enabled())
  {
    sfxprogress = gt_progress_timer_new("determining sequence length and "
                                        "number of special symbols");
  }
  if (gt_str_length(so->inputindex) > 0)
  {
    GtEncseqLoader *el;
    el = gt_encseq_loader_new();
    if (!so->fn2encopt.outdestab)
      gt_encseq_loader_do_not_require_des_tab(el);
    if (!so->fn2encopt.outsdstab)
      gt_encseq_loader_do_not_require_sds_tab(el);
    if (!so->fn2encopt.outssptab)
      gt_encseq_loader_do_not_require_ssp_tab(el);
    gt_encseq_loader_set_logger(el, logger);
    encseq = gt_encseq_loader_load(el, gt_str_get(so->inputindex), err);
    gt_encseq_loader_delete(el);
    if (encseq == NULL)
    {
      haserr = true;
    }
  } else
  {
    /* -tis is always implied */
    if (!so->fn2encopt.outdestab)
      gt_encseq_encoder_do_not_create_des_tab(ee);
    if (!so->fn2encopt.outsdstab)
      gt_encseq_encoder_do_not_create_sds_tab(ee);
    if (!so->fn2encopt.outssptab)
      gt_encseq_encoder_do_not_create_ssp_tab(ee);
    if (so->fn2encopt.isdna)
      gt_encseq_encoder_set_input_dna(ee);
    if (so->fn2encopt.isprotein)
      gt_encseq_encoder_set_input_protein(ee);
    if (so->fn2encopt.isplain)
      gt_encseq_encoder_set_input_preencoded(ee);
    gt_encseq_encoder_set_progresstimer(ee, sfxprogress);
    if (gt_encseq_encoder_use_symbolmap_file(ee,
                                             gt_str_get(so->fn2encopt.smap),
                                             err) != 0)
    {
      haserr = true;
    }
    if (!haserr) {
      if (gt_encseq_encoder_use_representation(ee,
                                               gt_str_get(so->fn2encopt.sat),
                                               err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr) {
      int rval;
      gt_encseq_encoder_set_logger(ee, logger);
      rval = gt_encseq_encoder_encode(ee, so->fn2encopt.filenametab,
                                      gt_str_get(so->fn2encopt.indexname),
                                      err);
      if (rval != 0)
        haserr = true;
      if (!haserr) {
        GtEncseqLoader *el = gt_encseq_loader_new();
        /* as we only construct the des, sds, and ssptable, but do not
           need it later during the construction of the enhanced suffix
           array, we do not load it again: we just load the information
           stored in the .esq file. It would be better to use the
           in memory constructed encseq rather than the one mapped into
           memory. This would require to free the des, sds, and ssptable. */
        gt_encseq_loader_do_not_require_des_tab(el);
        gt_encseq_loader_do_not_require_sds_tab(el);
        gt_encseq_loader_do_not_require_ssp_tab(el);
        encseq = gt_encseq_loader_load(el, gt_str_get(so->fn2encopt.indexname),
                                       err);
        if (!encseq)
          haserr = true;
        gt_encseq_loader_delete(el);
      }
    }
    if (encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_encseq_show_features(encseq,logger,false);
    if (so->readmode == GT_READMODE_COMPL ||
        so->readmode == GT_READMODE_REVCOMPL)
    {
      if (!gt_alphabet_is_dna(gt_encseq_alphabet(encseq)))
      {
        gt_error_set(err,"option -%s only can be used for DNA alphabets",
                          so->readmode == GT_READMODE_COMPL ? "cpl" : "rcl");
        haserr = true;
      }
    }
  }
  if (!haserr && so->outkystab && !so->outkyssort)
  {
    if (gt_extractkeysfromdesfile(gt_str_get(so->fn2encopt.indexname),
                                  false, logger, err) != 0)
    {
      haserr = true;
    }
  }
  prefixlength = so->prefixlength;
  sfxstrategy = so->sfxstrategy;
  sfxstrategy.ssortmaxdepth.defined = false;
  sfxstrategy.iteratorbasedkmerscanning = so->iteratorbasedkmerscanning;
  sfxstrategy.suftabasulongarray = so->suftabasulongarray;
  if (!haserr)
  {
    if (so->outsuftab || so->outbwttab || so->outlcptab || so->outbcktab ||
        !doesa)
    {
      unsigned int numofchars = gt_alphabet_num_of_chars(
                                          gt_encseq_alphabet(encseq));

      if (detpfxlenandmaxdepth(&prefixlength,
                              &sfxstrategy.ssortmaxdepth,
                              so,
                              numofchars,
                              gt_encseq_total_length(encseq),
                              logger,
                              err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (so->readmode != GT_READMODE_FORWARD)
      {
        gt_error_set(err,"option '-dir %s' only makes sense in combination "
                          "with at least one of the options -suf, -lcp, or "
                          "-bwt",
                          gt_readmode_show(so->readmode));
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
    if (initoutfileinfo(&outfileinfo,prefixlength,encseq,so,err) != 0)
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
        if (suffixeratorwithoutput(so->fn2encopt.indexname,
                                   &outfileinfo,
                                   encseq,
                                   so->readmode,
                                   prefixlength,
                                   so->numofparts,
                                   so->maximumspace,
                                   &sfxstrategy,
                                   sfxprogress,
                                   so->showprogress,
                                   logger,
                                   err) != 0)
        {
          haserr = true;
        }
      } else
      {
#ifndef S_SPLINT_S
        if (run_packedindexconstruction(logger,
                                        sfxprogress,
                                        so->showprogress,
                                        outfileinfo.outfpbcktab,
                                        so,
                                        prefixlength,
                                        encseq,
                                        &sfxstrategy,
                                        err) != 0)
        {
          haserr = true;
        }
#endif
      }
    }
  }
  gt_fa_fclose(outfileinfo.outfpsuftab);
  gt_fa_fclose(outfileinfo.outfpbwttab);
  gt_fa_fclose(outfileinfo.outfpbcktab);
  if (!haserr)
  {
    unsigned long numoflargelcpvalues,
          maxbranchdepth;

    if (outfileinfo.outlcpinfo == NULL)
    {
      numoflargelcpvalues = maxbranchdepth = 0;
    } else
    {
      numoflargelcpvalues = getnumoflargelcpvalues(outfileinfo.outlcpinfo);
      maxbranchdepth = getmaxbranchdepth(outfileinfo.outlcpinfo);
    }
    if (gt_outprjfile(gt_str_get(so->fn2encopt.indexname),
                      so->readmode,
                      encseq,
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
  gt_freeOutlcptab(outfileinfo.outlcpinfo);
  gt_encseq_delete(encseq);
  encseq = NULL;
  if (!haserr && so->outkystab && so->outkyssort)
  {
    if (gt_extractkeysfromdesfile(gt_str_get(so->fn2encopt.indexname), true,
                                  logger, err) != 0)
    {
      haserr = true;
    }
  }
  if (sfxprogress != NULL)
  {
    gt_progress_timer_start_new_state(sfxprogress,NULL,stdout);
    gt_progress_timer_delete(sfxprogress);
  }
  gt_encseq_encoder_delete(ee);
  return haserr ? -1 : 0;
}

void showgetencodedcharcounters(void);

int gt_parseargsandcallsuffixerator(bool doesa,int argc,
                                const char **argv,GtError *err)
{
  Suffixeratoroptions so;
  int retval;
  bool haserr = false;

  gt_error_check(err);
  retval = gt_suffixeratoroptions(&so,doesa,argc,argv,err);
  if (retval == 0)
  {
    GtLogger *logger = gt_logger_new(so.beverbose,
                                     GT_LOGGER_DEFLT_PREFIX, stdout);

    gt_logger_log(logger,"sizeof (unsigned long)=%lu",
                  (unsigned long) (sizeof (unsigned long) * CHAR_BIT));
    if (runsuffixerator(doesa,&so,logger,err) < 0)
    {
      haserr = true;
    }
    gt_logger_delete(logger);
    logger = NULL;
  } else
  {
    if (retval < 0)
    {
      haserr = true;
    }
  }
  gt_wrapsfxoptions(&so);
  /*
  showgetencodedcharcounters();
  printf("countgt_encseq_compare_viatwobitencoding=%lu\n",
          countgt_encseq_compare_viatwobitencoding_get());
  printf("# inverse_order = %lu\n",inverse_order);
  printf("# no_order = %lu\n",no_order);
  printf("# correct_order = %lu\n",correct_order);
  */
  return haserr ? -1 : 0;
}
