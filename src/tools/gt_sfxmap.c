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

#include "core/encseq_metadata.h"
#include "core/error.h"
#include "core/logger.h"
#include "core/option_api.h"
#include "core/str_array_api.h"
#include "match/echoseq.h"
#include "match/esa-seqread.h"
#include "match/esa-map.h"
#include "match/eis-voiditf.h"
#include "match/index_options.h"
#include "match/optionargmode.h"
#include "match/pckdfs.h"
#include "match/sfx-bentsedg.h"
#include "match/sfx-diffcov.h"
#include "match/sfx-strategy.h"
#include "match/sfx-suffixgetset.h"
#include "match/sfx-suftaborder.h"
#include "match/test-mappedstr.pr"
#include "match/twobits2kmers.h"
#include "match/esa-lcpintervals.h"
#include "tools/gt_sfxmap.h"

typedef struct
{
  bool usestream,
       verbose,
       cmpsuf,
       cmplcp,
       inputtis,
       inputsuf,
       inputdes,
       inputsds,
       inputbwt,
       inputlcp,
       inputbck,
       inputssp,
       enumlcpitvs,
       enumlcpitvtree,
       enumlcpitvtreeBU,
       diffcovercheck;
  unsigned long delspranges;
  GtStr *esaindexname,
        *pckindexname;
  unsigned int sortmaxdepth;
  GtStrArray *algbounds,
             *streamesq;
} Sfxmapoptions;

static void deletethespranges(const GtEncseq *encseq,
                              unsigned long delspranges)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  unsigned long rangewidth, nextpos = 0, totallength;
  const unsigned long fastawidth = 70UL;

  sri = gt_specialrangeiterator_new(encseq,true);
  printf(">\n");
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_assert(range.end > range.start);
    rangewidth = range.end - range.start;
    if (rangewidth > (unsigned long) delspranges)
    {
      if (range.start == 0)
      {
        nextpos = range.end;
      } else
      {
        if (range.start > nextpos)
        {
          gt_encseq2symbolstring(stdout,
                              encseq,
                              GT_READMODE_FORWARD,
                              nextpos,
                              range.start + delspranges - nextpos,
                              fastawidth);
          nextpos = range.end;
        }
      }
    }
  }
  totallength = gt_encseq_total_length(encseq);
  if (nextpos < totallength-1)
  {
    gt_encseq2symbolstring(stdout,
                        encseq,
                        GT_READMODE_FORWARD,
                        nextpos,
                        totallength - nextpos,
                        fastawidth);
  }
  gt_specialrangeiterator_delete(sri);
}

static void *gt_sfxmap_arguments_new(void)
{
  Sfxmapoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->esaindexname = gt_str_new();
  arguments->algbounds = gt_str_array_new();
  arguments->pckindexname = gt_str_new();
  arguments->streamesq = gt_str_array_new();
  return arguments;
}

static void gt_sfxmap_arguments_delete(void *tool_arguments)
{
  Sfxmapoptions *arguments = tool_arguments;

  if (arguments != NULL)
  {
    gt_str_delete(arguments->esaindexname);
    gt_str_delete(arguments->pckindexname);
    gt_str_array_delete(arguments->streamesq);
    gt_str_array_delete(arguments->algbounds);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sfxmap_option_parser_new(void *tool_arguments)
{
  Sfxmapoptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionstream, *optionverbose, *optionbck, *optionsuf,
         *optiondes, *optionsds, *optionbwt, *optionlcp, *optiontis, *optionssp,
         *optiondelspranges, *optionpckindex, *optionesaindex,
         *optioncmpsuf, *optioncmplcp, *optionstreamesq,
         *optionsortmaxdepth, *optionalgbounds, *optiondiffcov,
         *optionenumlcpitvs, *optionenumlcpitvtree, *optionenumlcpitvtreeBU;

  gt_assert(arguments != NULL);
  op = gt_option_parser_new("[options]",
                            "Map or Stream <indexname> and check consistency.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  optionesaindex = gt_option_new_string("esa",
                                        "Specify index (enhanced suffix array)",
                                        arguments->esaindexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                        "Specify index (packed index)",
                                        arguments->pckindexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);

  optionstreamesq = gt_option_new_stringarray("stream-esq",
                                              "Stream the encoded sequence",
                                              arguments->streamesq);
  gt_option_parser_add_option(op, optionstreamesq);

  optionsortmaxdepth = gt_option_new_uint("sortmaxdepth",
                                          "sort suffixes up to some depth",
                                          &arguments->sortmaxdepth,0);
  gt_option_parser_add_option(op, optionsortmaxdepth);

  optionalgbounds = gt_option_new_stringarray("algbds",
                                "length boundaries for the different "
                                "algorithms to sort buckets of suffixes\n"
                                "first number: maxbound for insertion sort\n"
                                "second number: maxbound for blindtrie sort\n"
                                "third number: maxbound for counting sort\n",
                                arguments->algbounds);
  gt_option_parser_add_option(op, optionalgbounds);

  gt_option_is_mandatory_either_3(optionesaindex,optionpckindex,
                                  optionstreamesq);

  optionstream = gt_option_new_bool("stream","stream the index",
                                    &arguments->usestream,false);
  gt_option_parser_add_option(op, optionstream);

  optiondelspranges = gt_option_new_ulong("delspranges",
                                          "delete ranges of special values",
                                           &arguments->delspranges,
                                           0);
  gt_option_parser_add_option(op, optiondelspranges);

  optiontis = gt_option_new_bool("tis","input the transformed input sequence",
                                 &arguments->inputtis,
                                 false);
  gt_option_parser_add_option(op, optiontis);

  optiondes = gt_option_new_bool("des","input the descriptions",
                                 &arguments->inputdes,
                                 false);
  gt_option_parser_add_option(op, optiondes);

  optionsds = gt_option_new_bool("sds","input the description end positions",
                                 &arguments->inputsds,
                                 false);
  gt_option_parser_add_option(op, optionsds);

  optionsuf = gt_option_new_bool("suf","input the suffix array",
                                 &arguments->inputsuf,
                                 false);
  gt_option_parser_add_option(op, optionsuf);

  optionlcp = gt_option_new_bool("lcp","input the lcp-table",
                                 &arguments->inputlcp,
                                 false);
  gt_option_parser_add_option(op, optionlcp);

  optionbwt = gt_option_new_bool("bwt",
                                 "input the Burrows-Wheeler Transformation",
                                 &arguments->inputbwt,
                                 false);
  gt_option_parser_add_option(op, optionbwt);

  optionbck = gt_option_new_bool("bck","input the bucket table",
                                 &arguments->inputbck,
                                 false);
  gt_option_parser_add_option(op, optionbck);

  optioncmpsuf = gt_option_new_bool("cmpsuf","compare pck derived suftab with "
                                    "esa-suftab",
                                    &arguments->cmpsuf,
                                    false);
  gt_option_parser_add_option(op, optioncmpsuf);

  optioncmplcp = gt_option_new_bool("cmplcp","compare pck derived lcptab with "
                                    "esa-lcptab",
                                    &arguments->cmplcp,
                                    false);
  gt_option_parser_add_option(op, optioncmplcp);
  gt_option_exclude(optioncmpsuf,optioncmplcp);

  optionssp = gt_option_new_bool("ssp","input the sequence separator table",
                                 &arguments->inputssp,
                                 false);
  gt_option_parser_add_option(op, optionssp);

  optiondiffcov = gt_option_new_bool("diffcover","check difference covers",
                                     &arguments->diffcovercheck,false);
  gt_option_parser_add_option(op, optiondiffcov);

  optionenumlcpitvs = gt_option_new_bool("enumlcpitvs",
                                         "enumerate the lcp-intervals",
                                         &arguments->enumlcpitvs,false);
  gt_option_parser_add_option(op, optionenumlcpitvs);

  optionenumlcpitvtree = gt_option_new_bool("enumlcpitvtree",
                                            "enumerate the lcp-interval tree",
                                            &arguments->enumlcpitvtree,false);
  gt_option_parser_add_option(op, optionenumlcpitvtree);

  optionenumlcpitvtreeBU = gt_option_new_bool("enumlcpitvtreeBU",
                                      "enumerate the lcp-interval tree "
                                      "(using a bottom-up strategy)",
                                       &arguments->enumlcpitvtreeBU,false);
  gt_option_parser_add_option(op, optionenumlcpitvtreeBU);

  optionverbose = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, optionverbose);

  gt_option_exclude(optionenumlcpitvs,optionenumlcpitvtree);
  gt_option_exclude(optionenumlcpitvs,optionenumlcpitvtreeBU);
  gt_option_exclude(optionenumlcpitvtree,optionenumlcpitvtreeBU);
  gt_option_imply(optionlcp,optionsuf);
  gt_option_imply(optionenumlcpitvs,optionesaindex);
  gt_option_imply(optionsortmaxdepth,optionesaindex);
  gt_option_imply(optionalgbounds,optionsortmaxdepth);
  return op;
}

static int gt_sfxmap_arguments_check(GT_UNUSED int rest_argc,
                                     GT_UNUSED void *tool_arguments,
                                     GT_UNUSED GtError *err)
{
  return 0;
}

static void showcomparisonfailureESA(const char *filename,
                                     int line,
                                     const char *where,
                                     const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     const ESASuffixptr *suftab,
                                     unsigned long depth,
                                     unsigned long idx1,
                                     unsigned long idx2,
                                     int cmp,
                                     unsigned long maxlcp)
{
  fprintf(stderr,"ERROR: file \"%s\", line %d: ",filename,line);
  fprintf(stderr,"%s(%lu vs %lu"
                 " %lu=\"",
                       where,
                       idx1,
                       idx2,
                       ESASUFFIXPTRGET(suftab,idx1));
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    ESASUFFIXPTRGET(suftab,idx1),depth);
  fprintf(stderr,"\",\"");
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    ESASUFFIXPTRGET(suftab,idx2),depth);
  fprintf(stderr,"\"=%lu)=%d with maxlcp %lu\n",ESASUFFIXPTRGET(suftab,idx2),
                                                cmp,
                                                maxlcp);
}

static void gt_checkentiresuftab(const char *filename,
                                 int line,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 const ESASuffixptr *suftab,
                                 unsigned long numberofsuffixes,
                                 Sequentialsuffixarrayreader *ssar,
                                 bool specialsareequal,
                                 bool specialsareequalatdepth0,
                                 unsigned long depth,
                                 GT_UNUSED GtError *err)
{
  unsigned long idx, maxlcp,
                currentlcp = 0,
                totallength = gt_encseq_total_length(encseq);
  int cmp;
  GtEncseqReader *esr1, *esr2;
  bool haserr = false;

#ifdef INLINEDSequentialsuffixarrayreader
  GtUchar tmpsmalllcpvalue;
#else
  int retval;
#endif

  gt_error_check(err);
  gt_assert(!specialsareequal || specialsareequalatdepth0);
  if (numberofsuffixes == totallength+1)
  {
    GtBitsequence *startposoccurs;
    unsigned long countbitsset = 0;

    GT_INITBITTAB(startposoccurs,totallength+1);
    for (idx = 0; idx <= totallength; idx++)
    {
      if (GT_ISIBITSET(startposoccurs,ESASUFFIXPTRGET(suftab,idx)))
      {
        fprintf(stderr,"ERROR: suffix with startpos %lu"
                       " already occurs\n",
                        ESASUFFIXPTRGET(suftab,idx));
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      GT_SETIBIT(startposoccurs,ESASUFFIXPTRGET(suftab,idx));
      countbitsset++;
    }
    if (countbitsset != totallength+1)
    {
      fprintf(stderr,"ERROR: not all bits are set\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_free(startposoccurs);
  }
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  gt_assert(ESASUFFIXPTRGET(suftab,0) < totallength);
  for (idx = 1UL; !haserr && idx < numberofsuffixes; idx++)
  {
    if (idx < numberofsuffixes - 1)
    {
      gt_assert(ESASUFFIXPTRGET(suftab,idx) < totallength);
      cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                               readmode,
                                               &maxlcp,
                                               specialsareequal,
                                               specialsareequalatdepth0,
                                               depth,
                                               ESASUFFIXPTRGET(suftab,idx-1),
                                               ESASUFFIXPTRGET(suftab,idx),
                                               esr1,
                                               esr2);
      if (cmp > 0)
      {
        showcomparisonfailureESA(filename,
                              line,
                              "checkentiresuftab",
                              encseq,
                              readmode,
                              suftab,
                              depth,
                              idx-1,
                              idx,
                              cmp,
                              maxlcp);
        haserr = true;
        break;
      }
    } else
    {
      maxlcp = 0;
      if (numberofsuffixes == totallength+1)
      {
        gt_assert(ESASUFFIXPTRGET(suftab,idx) == totallength);
      }
    }
    if (ssar != NULL)
    {
#ifdef INLINEDSequentialsuffixarrayreader
      NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
#else
      retval = gt_nextSequentiallcpvalue(&currentlcp,ssar,err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
#endif
      if (maxlcp != currentlcp)
      {
        fprintf(stderr,"%lu: startpos=%lu, firstchar=%u, "
                "startpos=%lu,firstchar=%u",
                idx,
                ESASUFFIXPTRGET(suftab,idx-1),
                (unsigned int)
                gt_encseq_get_encoded_char(encseq,
                                           ESASUFFIXPTRGET(suftab,idx-1),
                                           readmode),
                ESASUFFIXPTRGET(suftab,idx),
                (ESASUFFIXPTRGET(suftab,idx) < totallength)
                   ? (unsigned int) gt_encseq_get_encoded_char(encseq,
                                                     ESASUFFIXPTRGET(suftab,
                                                                     idx),
                                                     readmode)
                   : SEPARATOR);
        fprintf(stderr,", maxlcp(bruteforce) = %lu != %lu(fast)\n",
                          maxlcp, currentlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  if (haserr)
  {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  /*
  printf("# gt_checkentiresuftab with mode 'specials are %s'\n",
               specialsareequal ? "equal" : "different");
  */
}

static int sfxmap_esa(const Sfxmapoptions *arguments, GtLogger *logger,
                      GtError *err)
{
  bool haserr = false;
  Suffixarray suffixarray;
  unsigned int demand = 0;

  if (arguments->inputtis || arguments->delspranges > 0 || arguments->inputsuf)
  {
    demand |= SARR_ESQTAB;
  }
  if (arguments->inputdes)
  {
    demand |= SARR_DESTAB;
  }
  if (arguments->inputsds)
  {
    demand |= SARR_SDSTAB;
  }
  if (arguments->inputsuf)
  {
    demand |= SARR_SUFTAB;
  }
  if (arguments->inputlcp)
  {
    demand |= SARR_LCPTAB;
  }
  if (arguments->inputbwt)
  {
    demand |= SARR_BWTTAB;
  }
  if (arguments->inputbck)
  {
    demand |= SARR_BCKTAB;
  }
  if (arguments->inputssp)
  {
    demand |= SARR_SSPTAB;
  }
  if ((arguments->usestream
         ? streamsuffixarray
         : gt_mapsuffixarray)(&suffixarray,
                              demand,
                              gt_str_get(arguments->esaindexname),
                              logger,
                              err) != 0)
  {
    haserr = true;
  }
  if (!haserr && suffixarray.encseq != NULL)
  {
    if (arguments->delspranges > 0)
    {
      deletethespranges(suffixarray.encseq,arguments->delspranges);
    } else
    {
      if (!haserr && arguments->inputsuf && !arguments->usestream)
      {
        Sequentialsuffixarrayreader *ssar;

        if (arguments->inputlcp)
        {
          ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                        gt_str_get(arguments->esaindexname),
                                        SARR_LCPTAB | SARR_ESQTAB,
                                        SEQ_scan,
                                        err);
        } else
        {
          ssar = NULL;
        }
        gt_logger_log(logger, "checkentiresuftab");
        gt_checkentiresuftab(__FILE__,
                             __LINE__,
                             suffixarray.encseq,
                             suffixarray.readmode,
                             suffixarray.suftab,
                             gt_encseq_total_length(suffixarray.encseq)+1,
                             ssar,
                             false, /* specialsareequal  */
                             false,  /* specialsareequalatdepth0 */
                             0,
                             err);
        if (ssar != NULL)
        {
          gt_freeSequentialsuffixarrayreader(&ssar);
        }
        gt_logger_log(logger, "okay");
      }
      if (!haserr && arguments->inputbwt)
      {
        unsigned long totallength, bwtdifferentconsecutive = 0, idx, longest;

        gt_assert(suffixarray.longest.defined);
        longest = suffixarray.longest.valueunsignedlong;
        printf("longest=%lu\n",(unsigned long) longest);
        totallength = gt_encseq_total_length(suffixarray.encseq);
        printf("totallength=%lu\n",(unsigned long) totallength);
        if (!arguments->usestream)
        {
          for (idx = (unsigned long) 1; idx<totallength; idx++)
          {
            if (suffixarray.bwttab[idx-1] != suffixarray.bwttab[idx] ||
                ISSPECIAL(suffixarray.bwttab[idx]))
            {
              bwtdifferentconsecutive++;
            }
          }
        } else
        {
          GtUchar prevcc;

          if (readnextGtUcharfromstream(&prevcc,&suffixarray.bwttabstream) == 1)
          {
            GtUchar cc;
            while (readnextGtUcharfromstream(&cc,&suffixarray.bwttabstream)
                   == 1)
            {
              if (prevcc != cc || ISSPECIAL(cc))
              {
                bwtdifferentconsecutive++;
              }
            }
          }
        }
        printf("bwtdifferentconsecutive=%lu (%.4f)\n",
               bwtdifferentconsecutive,
               (double) bwtdifferentconsecutive/totallength);
      }
    }
  }
  if (!haserr && arguments->inputdes && arguments->inputsds)
  {
    gt_logger_log(logger, "checkallsequencedescriptions");
    gt_encseq_check_descriptions(suffixarray.encseq);
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int comparelcpvalue(void *info,unsigned long lcp,GtError *err)
{
  int retval;
  unsigned long currentlcpvalue;
  Sequentialsuffixarrayreader *ssar = (Sequentialsuffixarrayreader *) info;

  retval = gt_nextSequentiallcpvalue(&currentlcpvalue,ssar,err);
  if (retval < 0)
  {
    return -1;
  }
  if (retval == 0)
  {
    gt_error_set(err,"missing lcpvalue value");
    return -1;
  }
  if (lcp != currentlcpvalue)
  {
    fprintf(stderr,"lcp = %lu != %lu = currentlcpvalue\n",lcp,currentlcpvalue);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  return 0;
}

static int sfxmap_pck(const Sfxmapoptions *arguments,GtLogger *logger,
                      GtError *err)
{
  bool haserr = false;
  FMindex *fmindex;
  unsigned long totallength = 0;
  unsigned int numofchars = 0;
  GtEncseqMetadata *encseqmetadata = NULL;
  Sequentialsuffixarrayreader *ssar;

  gt_assert(gt_str_length(arguments->pckindexname) > 0);
  fmindex = gt_loadvoidBWTSeqForSA(gt_str_get(arguments->pckindexname),false,
                                   err);
  if (fmindex == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    encseqmetadata = gt_encseq_metadata_new(gt_str_get(arguments->pckindexname),
                                            err);
    if (encseqmetadata == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_str_length(arguments->esaindexname) > 0 &&
        (arguments->cmpsuf || arguments->cmplcp))
    {
      ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                          gt_str_get(arguments->esaindexname),
                                          arguments->cmpsuf ? SARR_SUFTAB
                                                            : SARR_LCPTAB,
                                          SEQ_scan,
                                          err);
      if (ssar == NULL)
      {
        haserr = true;
      }
    } else
    {
      ssar = NULL;
    }
  }
  if (!haserr)
  {
    unsigned long idx, pos, numofnonspecials, currentsuffix;
    GtSpecialcharinfo specialcharinfo;
    Bwtseqpositioniterator *bspi;
    int retval;

    gt_assert(encseqmetadata != NULL);
    totallength = gt_encseq_metadata_total_length(encseqmetadata);
    specialcharinfo = gt_encseq_metadata_specialcharinfo(encseqmetadata);
    gt_assert(totallength >= specialcharinfo.specialcharacters);
    numofnonspecials = totallength - specialcharinfo.specialcharacters;
    bspi = gt_Bwtseqpositioniterator_new(fmindex,0,totallength);
    gt_logger_log(logger, "iterate over all suftab values");
    for (idx = 0; idx < numofnonspecials; idx++)
    {
      if (!gt_Bwtseqpositioniterator_next(&pos,bspi))
      {
        gt_error_set(err,"cannot decode enough symbols");
        haserr = true;
        break;
      }
      if (arguments->cmpsuf && ssar != NULL)
      {
        retval = gt_nextSequentialsuftabvalue(&currentsuffix,ssar);
        gt_assert(retval >= 0);
        if (retval == 0)
        {
          gt_error_set(err,"missing suftab value");
          haserr = true;
          break;
        }
        gt_assert(pos == currentsuffix);
      }
      /*printf("%lu: pos = %lu\n",idx,pos);*/
    }
    gt_assert(idx == numofnonspecials);
    gt_Bwtseqpositioniterator_delete(bspi);
  }
  if (!haserr)
  {
    GtAlphabet *alphabet;

    alphabet = gt_alphabet_new_from_file(gt_str_get(arguments->pckindexname),
                                         err);
    if (alphabet == NULL)
    {
      haserr = true;
    }
    numofchars = gt_alphabet_num_of_chars(alphabet);
    gt_alphabet_delete(alphabet);
  }
  if (!haserr)
  {
    if (arguments->cmplcp && ssar != NULL)
    {
      gt_logger_log(logger, "perform dfs traversal");
      if (gt_fmindex_dfstraverse(fmindex,numofchars,totallength,comparelcpvalue,
                                  (void *) ssar,err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_deletevoidBWTSeq(fmindex);
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_encseq_metadata_delete(encseqmetadata);
  return haserr ? -1 : 0;
}

static const Optionargmodedesc stream_esq_operation[] =
{
  {"stream_words","read words from stream",(unsigned int) BSRS_stream_words},
  {"stream_single","read single characters from word stream",
                   (unsigned int) BSRS_stream_single},
  {"reader_single","read single characters with encseq reader",
                   (unsigned int) BSRS_reader_single},
  {"stream_reader_single","read single characters with encseq reader and "
                          "from word stream",
                   (unsigned int) BSRS_stream_reader_single},
  {"reader_multi","read kmers with encseq reader",
                   (unsigned int) BSRS_reader_multi},
  {"stream_reader_multi","read kmers with encseq reader and from word stream",
                   (unsigned int) BSRS_stream_reader_multi},
  {"stream_reader_multi3","read kmers with encseq reader and from word stream",
                   (unsigned int) BSRS_stream_reader_multi3}
};

static int stream_esq(const Sfxmapoptions *arguments,GtError *err)
{
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  bool haserr = false;
  Bitstreamreadmode brsmode = BSRS_stream_single;
  int multiarg = 0;
  unsigned long streamesq_size = gt_str_array_size(arguments->streamesq);

  if (streamesq_size == 2UL || streamesq_size == 3UL)
  {
    int mode;

    mode = gt_optionargsetsingle(stream_esq_operation,
                                 sizeof stream_esq_operation/
                                 sizeof stream_esq_operation[0],
                                 "stream-esq",
                                 gt_str_array_get(arguments->streamesq,1UL),
                                 err);
    if (mode < 0)
    {
      haserr = true;
    } else
    {
      brsmode = (Bitstreamreadmode) mode;
      if ((brsmode == BSRS_stream_words ||
           brsmode == BSRS_stream_single ||
           brsmode == BSRS_reader_single ||
           brsmode == BSRS_stream_reader_single) &&
           streamesq_size != 2UL)
      {
        gt_error_set(err,"if option -streamesq has one of the arguments "
                         "stream_words stream_single reader_single "
                         "stream_reader_single then no other argument "
                         "is allowed");
        haserr = true;
      } else
      {
        if ((brsmode == BSRS_reader_multi ||
             brsmode == BSRS_stream_reader_multi ||
             brsmode == BSRS_stream_reader_multi3) &&
             streamesq_size != 3UL)
        {
          gt_error_set(err,"if option -streamesq has one of the arguments "
                           "stream_multi reader_multi "
                           "stream_reader_multi then one more argument "
                           "is required");
          haserr = true;
        }
      }
    }
    if (!haserr && streamesq_size == 3UL)
    {
      if (sscanf(gt_str_array_get(arguments->streamesq,2UL),"%d",&multiarg)
          != 1 || multiarg < 1)
      {
        gt_error_set(err,"if option -streamesq has three arguments, "
                         "then third argument must be positive integer");
        haserr = true;
      }
    }
  } else
  {
    gt_error_set(err,"option -streamesq must have two or three arguments");
    haserr = true;
  }
  if (!haserr)
  {
    el = gt_encseq_loader_new();
    encseq = gt_encseq_loader_load(el, gt_str_array_get(arguments->streamesq,0),
                                   err);
    if (encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_assert(multiarg >= 0);
    gt_encseq_faststream(encseq,brsmode,(unsigned int) multiarg);
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

typedef struct
{
  GtEncseq *encseq;
  GtReadmode readmode;
  unsigned int sortmaxdepth;
  GtSuffixsortspace *sssp;
} Checkunsortedrangeinfo;

static void sortmaxdepth_processunsortedrange(void *voiddcov,
                                              unsigned long subbucketleft,
                                              unsigned long width,
                                              GT_UNUSED unsigned long depth)
{
  Checkunsortedrangeinfo *curi = voiddcov;

  gt_checkifprefixesareidentical(__FILE__,
                                 __LINE__,
                                 curi->encseq,
                                 curi->readmode,
                                 curi->sssp,
                                 subbucketleft,
                                 width,
                                 (unsigned long) curi->sortmaxdepth);
}

static int performsortmaxdepth(const Sfxmapoptions *arguments,
                               GtLogger *logger,GtError *err)
{
  bool haserr = false;
  GtEncseqLoader *el = NULL;
  Checkunsortedrangeinfo curi;
  Sfxstrategy sfxstrategy;
  const char *indexname;

  indexname = gt_str_get(arguments->esaindexname);
  el = gt_encseq_loader_new();
  curi.encseq = gt_encseq_loader_load(el, indexname, err);
  if (curi.encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_logger_log(logger,"performsortmaxdepth(%s,%u)",
                  indexname,arguments->sortmaxdepth);
    defaultsfxstrategy(&sfxstrategy,
                       gt_encseq_bitwise_cmp_ok(curi.encseq) ? false : true);
    if (gt_str_array_size(arguments->algbounds) > 0 &&
        gt_parse_algbounds(&sfxstrategy,arguments->algbounds,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    unsigned long idx, totallength = gt_encseq_total_length(curi.encseq);

    curi.sssp = gt_suffixsortspace_new(totallength+1, totallength, true);
    for (idx=0; idx<=totallength; idx++)
    {
      gt_suffixsortspace_setdirect(curi.sssp,idx,idx);
    }
    curi.readmode = GT_READMODE_FORWARD;
    curi.sortmaxdepth = arguments->sortmaxdepth;
    gt_sortallsuffixesfromstart(curi.sssp,
                                totallength+1,
                                curi.encseq,
                                curi.readmode,
                                NULL,
                                arguments->sortmaxdepth,
                                &sfxstrategy,
                                sortmaxdepth_processunsortedrange,
                                &curi, /* voiddcov */
                                logger);
    gt_checksortedsuffixes(__FILE__,
                           __LINE__,
                           curi.encseq,
                           curi.readmode,
                           curi.sssp,
                           0,
                           totallength+1,
                           false, /* specialsareequal  */
                           false,  /* specialsareequalatdepth0 */
                           (unsigned long) curi.sortmaxdepth);
    gt_suffixsortspace_delete(curi.sssp,false);
  }
  gt_encseq_delete(curi.encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

static int run_diffcover_check(const Sfxmapoptions *arguments, GtError *err)
{
  int had_err = 0;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  const char *indexname;

  indexname = gt_str_get(arguments->esaindexname);
  el = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(el, indexname, err);
  if (!encseq)
    had_err = -1;

  if (!had_err) {
    GtReadmode readmode;
    for (readmode = GT_READMODE_FORWARD;
         readmode <= GT_READMODE_REVCOMPL; readmode++) {
      gt_differencecover_check(encseq, readmode);
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return had_err;
}

static int gt_sfxmap_runner(GT_UNUSED int argc,
                            GT_UNUSED const char **argv,
                            GT_UNUSED int parsed_args,
                            void *tool_arguments, GtError *err)
{
  bool haserr = false;
  Sfxmapoptions *arguments = tool_arguments;
  GtLogger *logger;

  gt_error_check(err);
  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (gt_str_length(arguments->esaindexname) > 0)
  {
    if (sfxmap_esa(arguments,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && gt_str_length(arguments->pckindexname) > 0)
  {
    if (sfxmap_pck(arguments,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && gt_str_array_size(arguments->streamesq) > 0)
  {
    if (stream_esq(arguments,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->sortmaxdepth > 0)
  {
    if (performsortmaxdepth(arguments,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->diffcovercheck)
  {
    if (run_diffcover_check(arguments, err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && (arguments->enumlcpitvs || arguments->enumlcpitvtree ||
                  arguments->enumlcpitvtreeBU))
  {
    if (gt_runenumlcpvalues(gt_str_get(arguments->esaindexname),
                            arguments->enumlcpitvs ? false : true,
                            arguments->enumlcpitvtreeBU,
                            logger, err) != 0)
    {
      haserr = true;
    }
  }
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}

GtTool* gt_sfxmap(void)
{
  return gt_tool_new(gt_sfxmap_arguments_new,
                     gt_sfxmap_arguments_delete,
                     gt_sfxmap_option_parser_new,
                     gt_sfxmap_arguments_check,
                     gt_sfxmap_runner);
}
