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

#include "core/error.h"
#include "core/logger.h"
#include "core/option.h"
#include "core/encseq_metadata.h"
#include "match/echoseq.h"
#include "match/sarr-def.h"
#include "match/esa-seqread.h"
#include "match/esa-map.h"
#include "match/eis-voiditf.h"
#include "match/pckdfs.h"
#include "match/test-mappedstr.pr"
#include "tools/gt_sfxmap.h"

typedef struct
{
  bool usestream,
       verbose,
       inputtis,
       inputsuf,
       inputdes,
       inputsds,
       inputbwt,
       inputlcp,
       inputbck,
       inputssp;
  unsigned long scantrials,
                multicharcmptrials,
                delspranges;
  GtStr *pckindexname, *esaindexname;
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
  arguments->pckindexname = gt_str_new();
  arguments->esaindexname = gt_str_new();
  return arguments;
}

static void gt_sfxmap_arguments_delete(void *tool_arguments)
{
  Sfxmapoptions *arguments = tool_arguments;

  if (arguments != NULL)
  {
    gt_str_delete(arguments->pckindexname);
    gt_str_delete(arguments->esaindexname);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sfxmap_option_parser_new(void *tool_arguments)
{
  Sfxmapoptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionstream, *optionverbose, *optionscantrials,
         *optionmulticharcmptrials, *optionbck, *optionsuf,
         *optiondes, *optionsds, *optionbwt, *optionlcp, *optiontis, *optionssp,
         *optiondelspranges, *optionpckindex, *optionesaindex;

  gt_assert(arguments != NULL);
  op = gt_option_parser_new("[options]",
                            "Map or Stream <indexname> and check consistency.");
  gt_option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");

  optionesaindex = gt_option_new_string("esa",
                                        "Specify index (enhanced suffix array)",
                                        arguments->esaindexname, NULL);
  gt_option_parser_add_option(op, optionesaindex);

  optionpckindex = gt_option_new_string("pck",
                                        "Specify index (packed index)",
                                        arguments->pckindexname, NULL);
  gt_option_parser_add_option(op, optionpckindex);
  gt_option_is_mandatory_either(optionesaindex,optionpckindex);

  optionstream = gt_option_new_bool("stream","stream the index",
                                 &arguments->usestream,false);
  gt_option_parser_add_option(op, optionstream);

  optionscantrials = gt_option_new_ulong("scantrials",
                                         "specify number of scan trials",
                                         &arguments->scantrials,0);
  gt_option_parser_add_option(op, optionscantrials);

  optionmulticharcmptrials
    = gt_option_new_ulong("multicharcmptrials",
                          "specify number of multichar cmp trials",
                          &arguments->multicharcmptrials,0);
  gt_option_parser_add_option(op, optionmulticharcmptrials);

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

  optionssp = gt_option_new_bool("ssp","input the sequence separator table",
                                 &arguments->inputssp,
                                 false);
  gt_option_parser_add_option(op, optionssp);

  optionverbose = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, optionverbose);

  gt_option_imply(optionlcp,optionsuf);
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
                                 GtError *err)
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

static int sfxmap_esa(Sfxmapoptions *arguments,GtError *err)
{
  bool haserr = false;
  Suffixarray suffixarray;
  GtLogger *logger;
  unsigned int demand = 0;

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
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
      if (!haserr && arguments->inputtis)
      {
        int readmode;

        for (readmode = 0; readmode < 4; readmode++)
        {
          if (gt_alphabet_is_dna(
            gt_encseq_alphabet(suffixarray.encseq)) ||
               ((GtReadmode) readmode) == GT_READMODE_FORWARD ||
               ((GtReadmode) readmode) == GT_READMODE_REVERSE)
          {
            gt_logger_log(logger, "testencseq(readmode=%s)",
                                   gt_readmode_show((GtReadmode) readmode));
            if (gt_encseq_check_consistency(suffixarray.encseq,
                               gt_encseq_filenames(suffixarray.encseq),
                               (GtReadmode) readmode,
                               arguments->scantrials,
                               arguments->multicharcmptrials,
                               arguments->inputssp,
                               err) != 0)
            {
              haserr = true;
              break;
            }
          }
        }
      }
      if (!haserr && arguments->inputtis)
      {
        gt_logger_log(logger, "checkspecialrangesfast");
        if (gt_encseq_check_specialranges(suffixarray.encseq) != 0)
        {
          haserr = true;
        }
      }
      if (!haserr && arguments->inputtis)
      {
        gt_logger_log(logger, "gt_encseq_check_markpos");
        gt_encseq_check_markpos(suffixarray.encseq);
      }
      if (!haserr && arguments->inputtis &&
          suffixarray.readmode == GT_READMODE_FORWARD &&
          suffixarray.prefixlength > 0)
      {
        gt_logger_log(logger, "verifymappedstr");
        if (gt_verifymappedstr(suffixarray.encseq,suffixarray.prefixlength,
                               err) != 0)
        {
          haserr = true;
        }
      }
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
  if (!haserr && arguments->inputdes)
  {
    gt_logger_log(logger, "checkallsequencedescriptions");
    gt_encseq_check_descriptions(suffixarray.encseq);
  }
  gt_freesuffixarray(&suffixarray);
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}

static int sfxmap_pck(Sfxmapoptions *arguments,GtError *err)
{
  bool haserr = false;
  FMindex *fmindex;
  unsigned long totallength = 0;
  unsigned int numofchars = 0;
  GtEncseqMetadata *encseqmetadata = NULL;
  Sequentialsuffixarrayreader *ssar;

  gt_assert(gt_str_length(arguments->pckindexname));
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
    if (gt_str_length(arguments->esaindexname) > 0)
    {
      ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                          gt_str_get(arguments->esaindexname),
                                          SARR_SUFTAB,
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
    for (idx = 0; idx < numofnonspecials; idx++)
    {
      if (!gt_Bwtseqpositioniterator_next(&pos,bspi))
      {
        gt_error_set(err,"cannot decode enough symbols");
        haserr = true;
        break;
      }
      if (ssar != NULL)
      {
        retval = gt_nextSequentialsuftabvalue(&currentsuffix,ssar);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          gt_error_set(err,"missing suftab values");
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
  gt_fmindex_dfstraverse(fmindex,numofchars,totallength);
  gt_deletevoidBWTSeq(fmindex);
  if (ssar != NULL)
  {
    gt_freeSequentialsuffixarrayreader(&ssar);
  }
  gt_encseq_metadata_delete(encseqmetadata);
  return haserr ? -1 : 0;
}

static int gt_sfxmap_runner(GT_UNUSED int argc,
                            GT_UNUSED const char **argv,
                            GT_UNUSED int parsed_args,
                            void *tool_arguments, GtError *err)
{
  Sfxmapoptions *arguments = tool_arguments;

  gt_error_check(err);
  if (gt_str_length(arguments->esaindexname) > 0)
  {
    return sfxmap_esa(arguments,err);
  }
  if (gt_str_length(arguments->pckindexname) > 0)
  {
    return sfxmap_pck(arguments,err);
  }
  return 0;
}

GtTool* gt_sfxmap(void)
{
  return gt_tool_new(gt_sfxmap_arguments_new,
                     gt_sfxmap_arguments_delete,
                     gt_sfxmap_option_parser_new,
                     gt_sfxmap_arguments_check,
                     gt_sfxmap_runner);
}
