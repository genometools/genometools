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
#include "core/unused_api.h"
#include "core/format64.h"
#include "core/fa.h"
#include "core/mathsupport.h"
#include "match/echoseq.h"
#include "match/eis-voiditf.h"
#include "match/esa-lcpintervals.h"
#include "match/esa-map.h"
#include "match/esa-seqread.h"
#include "match/esa-spmitvs.h"
#include "match/index_options.h"
#include "match/optionargmode.h"
#include "match/pckdfs.h"
#include "match/sfx-bentsedg.h"
#include "match/sfx-diffcov.h"
#include "match/sfx-strategy.h"
#include "match/sfx-suffixgetset.h"
#include "match/sfx-suftaborder.h"
#include "match/sfx-linlcp.h"
#include "match/sfx-lwcheck.h"
#include "match/twobits2kmers.h"
#include "match/esa-fileend.h"
#include "tools/gt_sfxmap.h"

typedef struct
{
  bool usestream,
       bfcheck,
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
       diffcovercheck,
       wholeleafcheck,
       compressedesa,
       compresslcp,
       spmitv,
       ownencseq2file;
  GtUword delspranges;
  GtStr *esaindexname,
        *pckindexname;
  unsigned int sortmaxdepth,
               scanesa;
  GtStrArray *algbounds,
             *streamesq;
} Sfxmapoptions;

static void gt_sfxmap_deletethespranges(const GtEncseq *encseq,
                                        GtUword delspranges)
{
  GtSpecialrangeiterator *sri;
  GtRange range;
  GtUword rangewidth, nextpos = 0, totallength;
  const GtUword fastawidth = 70UL;

  sri = gt_specialrangeiterator_new(encseq,true);
  printf(">\n");
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_assert(range.end > range.start);
    rangewidth = range.end - range.start;
    if (rangewidth > (GtUword) delspranges)
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
  arguments->pckindexname = gt_str_new();
  arguments->streamesq = gt_str_array_new();
  arguments->algbounds = gt_str_array_new();
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
         *optionwholeleafcheck,
         *optionenumlcpitvs, *optionenumlcpitvtree, *optionenumlcpitvtreeBU,
         *optionscanesa, *optionspmitv, *optionownencseq2file,
         *optionbfcheck, *optioncompressedesa,
         *optioncompresslcp;

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

  optionstreamesq = gt_option_new_string_array("stream-esq",
                                               "Stream the encoded sequence",
                                               arguments->streamesq);
  gt_option_parser_add_option(op, optionstreamesq);

  optionsortmaxdepth = gt_option_new_uint("sortmaxdepth",
                                          "sort suffixes up to some depth",
                                          &arguments->sortmaxdepth,0);
  gt_option_parser_add_option(op, optionsortmaxdepth);

  optionalgbounds = gt_option_new_string_array("algbds",
                                "length boundaries for the different "
                                "algorithms to sort buckets of suffixes\n"
                                "first number: maxbound for insertion sort\n"
                                "second number: maxbound for blindtrie sort\n"
                                "third number: maxbound for counting sort",
                                arguments->algbounds);
  gt_option_parser_add_option(op, optionalgbounds);

  gt_option_is_mandatory_either_3(optionesaindex,optionpckindex,
                                  optionstreamesq);

  optionstream = gt_option_new_bool("stream","stream the index",
                                    &arguments->usestream,false);
  gt_option_parser_add_option(op, optionstream);

  optionbfcheck = gt_option_new_bool("bfcheck",
                                     "perform check by brute force algorithm "
                                     "(this can be slow if lcps are long)",
                                     &arguments->bfcheck,false);
  gt_option_parser_add_option(op, optionbfcheck);

  optiondelspranges = gt_option_new_uword("delspranges",
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

  optionwholeleafcheck = gt_option_new_bool("wholeleafcheck",
                                            "check if all whole leaves are "
                                            "present",
                                            &arguments->wholeleafcheck,false);
  gt_option_parser_add_option(op, optionwholeleafcheck);

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

  optionscanesa = gt_option_new_uint("scanesa",
                                     "scan suftab and lcptab",
                                     &arguments->scanesa,0);
  gt_option_parser_add_option(op, optionscanesa);

  optionspmitv = gt_option_new_bool("spmitv",
                                    "determine distribution of intervals "
                                    "with whole leaves",
                                    &arguments->spmitv,false);
  gt_option_parser_add_option(op, optionspmitv);

  optionownencseq2file = gt_option_new_bool("ownencseq2file",
                                            "write own encseq to file",
                                            &arguments->ownencseq2file,false);
  gt_option_parser_add_option(op, optionownencseq2file);
  gt_option_imply(optionownencseq2file, optionesaindex);

  optioncompressedesa = gt_option_new_bool("compressedesa",
                                           "input the compressed tables of the "
                                           "enhanced suffix array",
                                           &arguments->compressedesa,
                                           false);
  gt_option_parser_add_option(op, optioncompressedesa);
  gt_option_imply(optioncompressedesa, optionesaindex);

  optioncompresslcp = gt_option_new_bool("compresslcp",
                                          "compress the lcp value into byte "
                                          "compressed form",
                                          &arguments->compresslcp,
                                          false);
  gt_option_parser_add_option(op, optioncompresslcp);
  gt_option_imply(optioncompresslcp, optionesaindex);

  optionverbose = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, optionverbose);

  gt_option_exclude(optionenumlcpitvs,optionenumlcpitvtree);
  gt_option_exclude(optioncompresslcp,optioncompressedesa);
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

static void gt_sfxmap_showcomparisonfailureESA(const char *filename,
                                     int line,
                                     const char *where,
                                     const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     const ESASuffixptr *suftab,
                                     GtUword depth,
                                     GtUword idx1,
                                     GtUword idx2,
                                     int cmp,
                                     GtUword maxlcp)
{
  fprintf(stderr,"ERROR: file \"%s\", line %d: ",filename,line);
  fprintf(stderr,"%s("GT_WU" vs "GT_WU""
                 " "GT_WU"=\"",
                       where,
                       idx1,
                       idx2,
                       ESASUFFIXPTRGET(suftab,idx1));
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    ESASUFFIXPTRGET(suftab,idx1),depth);
  fprintf(stderr,"\",\"");
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    ESASUFFIXPTRGET(suftab,idx2),depth);
  fprintf(stderr,"\"="GT_WU")=%d with maxlcp "GT_WU"\n",
          ESASUFFIXPTRGET(suftab,idx2),cmp,maxlcp);
}

static GtUword gt_sfxmap_determinenumberofwholeleaves(
                                                  const GtEncseq *encseq,
                                                  GtReadmode readmode)
{
  GtUword idx, wholeleafcount = 0, totallength;
  GtEncseqReader* esr;
  GtUchar cc;
  bool sequencestart = true;

  esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
  totallength = gt_encseq_total_length(encseq);
  for (idx = 0; idx < totallength; idx++)
  {
    cc = gt_encseq_reader_next_encoded_char(esr);
    if (cc == (GtUchar) SEPARATOR)
    {
      sequencestart = true;
    } else
    {
      if (sequencestart)
      {
        if (cc != (GtUchar) WILDCARD)
        {
          wholeleafcount++;
        }
        sequencestart = false;
      }
    }
  }
  gt_encseq_reader_delete(esr);
  return wholeleafcount;
}

static int gt_sfxmap_comparefullsuffixes(const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         GtUword *maxlcp,
                                         GtUword start1,
                                         GtUword start2,
                                         GtEncseqReader *esr1,
                                         GtEncseqReader *esr2)
{
  GtUchar cc1, cc2;
  GtUword pos1, pos2, totallength;
  int retval;

  totallength = gt_encseq_total_length(encseq);
  if (esr1 != NULL && esr2 != NULL)
  {
    gt_encseq_reader_reinit_with_readmode(esr1,encseq,readmode,start1);
    gt_encseq_reader_reinit_with_readmode(esr2,encseq,readmode,start2);
  } else
  {
    gt_assert(esr1 == NULL && esr2 == NULL);
  }
  for (pos1=start1, pos2=start2; /* Nothing */; pos1++, pos2++)
  {
    if (pos1 >= totallength)
    {
      cc1 = (GtUchar) SEPARATOR;
    } else
    {
      cc1 = (esr1 != NULL) ? gt_encseq_reader_next_encoded_char(esr1)
                           : gt_encseq_get_encoded_char(encseq,pos1,readmode);
    }
    if (pos2 >= totallength)
    {
      cc2 = (GtUchar) SEPARATOR;
    } else
    {
      cc2 = (esr2 != NULL) ? gt_encseq_reader_next_encoded_char(esr2)
                           : gt_encseq_get_encoded_char(encseq,pos2,readmode);
    }
    if (ISSPECIAL(cc1))
    {
      if (ISSPECIAL(cc2))
      {
        if (pos1 < pos2)
        {
          *maxlcp = pos1  - start1;
          retval = -1; /* a < b */
          break;
        }
        if (pos1 > pos2)
        {
          *maxlcp = pos1 - start1;
          retval = 1; /* a > b */
          break;
        }
        *maxlcp = pos1 - start1 + 1;
        retval = 0; /* a = b */
        break;
      }
      *maxlcp = pos1 - start1;
      retval = 1; /* a > b */
      break;
    } else
    {
      if (ISSPECIAL(cc2))
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 < cc2)
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 > cc2)
      {
        *maxlcp = pos1 - start1;
        retval = 1; /* a > b */
        break;
      }
    }
  }
  return retval;
}

static int gt_sfxmap_checkentiresuftab(const char *filename,
                                int line,
                                const GtEncseq *encseq,
                                GtReadmode readmode,
                                const ESASuffixptr *suftab,
                                GtUword numberofsuffixes,
                                bool wholeleafcheck,
                                Sequentialsuffixarrayreader *ssar,
                                GT_UNUSED bool specialsareequal,
                                GT_UNUSED bool specialsareequalatdepth0,
                                GtUword depth,
                                GT_UNUSED GtError *err)
{
  GtUword idx, maxlcp,
                currentlcp = 0,
                countbitsset = 0,
                wholeleafcount = 0,
                totallength = gt_encseq_total_length(encseq);
  int cmp;
  GtEncseqReader *esr1, *esr2;
  bool haserr = false;
  GtBitsequence *startposoccurs;
  GtEncseqReader *esr;
  /*
#define MAXDIST 100
  GtUword countdist[MAXDIST+1] = {0};
  */

  gt_error_check(err);
  gt_assert(!specialsareequal || specialsareequalatdepth0);
  if (numberofsuffixes == 0)
  {
    return 0;
  }
  GT_INITBITTAB(startposoccurs,totallength+1);
  esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    GtUword position = ESASUFFIXPTRGET(suftab,idx);
    if (GT_ISIBITSET(startposoccurs,position))
    {
      fprintf(stderr,"ERROR: suffix with startpos "GT_WU" already occurs\n",
              ESASUFFIXPTRGET(suftab,idx));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    GT_SETIBIT(startposoccurs,position);
    countbitsset++;
    if (wholeleafcheck)
    {
      if (position == 0 || gt_encseq_position_is_separator(encseq,
                                                           position - 1,
                                                           readmode))
      {
        /*printf("whole "GT_WU"\n",position);*/
        wholeleafcount++;
      }
    }
  }
  gt_encseq_reader_delete(esr);
  gt_free(startposoccurs);
  if (wholeleafcheck)
  {
    GtUword expectednumofwholeleaves
      = gt_sfxmap_determinenumberofwholeleaves(encseq,readmode);
    if (wholeleafcount != expectednumofwholeleaves)
    {
      fprintf(stderr,
              "wholeleafcount="GT_WU" != "GT_WU"=expectednumofwholeleaves\n",
              wholeleafcount,expectednumofwholeleaves);
      exit(EXIT_FAILURE);
    }
  }
  if (numberofsuffixes == totallength + 1 && countbitsset != numberofsuffixes)
  {
    fprintf(stderr,"ERROR: not all bits are set\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  gt_assert(ESASUFFIXPTRGET(suftab,0) < totallength);
  for (idx = 1UL; !haserr && idx < numberofsuffixes; idx++)
  {
    if (idx < totallength)
    {
      gt_assert(ESASUFFIXPTRGET(suftab,idx) < totallength);
      cmp = gt_sfxmap_comparefullsuffixes(encseq,
                                          readmode,
                                          &maxlcp,
                                          ESASUFFIXPTRGET(suftab,idx-1),
                                          ESASUFFIXPTRGET(suftab,idx),
                                          esr1,
                                          esr2);
      if (cmp >= 0)
      {
        gt_sfxmap_showcomparisonfailureESA(filename,
                                 line,
                                 __func__,
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
      SSAR_NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
      if (maxlcp != currentlcp)
      {
        fprintf(stderr,""GT_WU": startpos="GT_WU", firstchar=%u, "
                "startpos="GT_WU",firstchar=%u",
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
        fprintf(stderr,", maxlcp(bruteforce) = "GT_WU" != "GT_WU"(fast)\n",
                          maxlcp, currentlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  return haserr ? -1 : 0;
}

static GtUchar sfx_accesschar_encseq(const void *encseq,GtUword position,
                                     GtReadmode readmode)
{
  return gt_encseq_get_encoded_char((const GtEncseq *) encseq,
                                    position,
                                    readmode);
}

GtUword sfx_charcount_encseq(const void *encseq,GtUchar idx)
{
  return gt_encseq_charcount((const GtEncseq *) encseq, idx);
}

static int gt_sfxmap_esa(const Sfxmapoptions *arguments, GtLogger *logger,
                         GtError *err)
{
  bool haserr = false;
  Suffixarray suffixarray;
  unsigned int demand = 0;

  gt_error_check(err);
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
      gt_sfxmap_deletethespranges(suffixarray.encseq,arguments->delspranges);
    } else
    {
      GtUword totallength = gt_encseq_total_length(suffixarray.encseq);
      if (!haserr && arguments->inputsuf && !arguments->usestream)
      {
        if (suffixarray.numberofallsortedsuffixes != totallength + 1 ||
            arguments->bfcheck)
        {
          Sequentialsuffixarrayreader *ssar;

          if (arguments->inputlcp)
          {
            ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                          gt_str_get(arguments->esaindexname),
                                          SARR_LCPTAB | SARR_ESQTAB,
                                          true,
                                          logger,
                                          err);
          } else
          {
            ssar = NULL;
          }
          gt_logger_log(logger, __func__);
          if (gt_sfxmap_checkentiresuftab(__FILE__,
                                   __LINE__,
                                   suffixarray.encseq,
                                   suffixarray.readmode,
                                   suffixarray.suftab,
                                   suffixarray.numberofallsortedsuffixes,
                                   arguments->wholeleafcheck,
                                   ssar,
                                   false, /* specialsareequal  */
                                   false,  /* specialsareequalatdepth0 */
                                   0,
                                   err) != 0)
          {
            fprintf(stderr,"%s\n",__func__);
            exit(GT_EXIT_PROGRAMMING_ERROR);
          }
          if (ssar != NULL)
          {
            gt_freeSequentialsuffixarrayreader(&ssar);
          }
        } else
        {
          if (suffixarray.numberofallsortedsuffixes == totallength + 1)
          {
            gt_suftab_lightweightcheck(sfx_accesschar_encseq,
                                       sfx_charcount_encseq,
                                       suffixarray.encseq,
                                       suffixarray.readmode,
                                       totallength,
                                       gt_encseq_alphabetnumofchars(
                                                 suffixarray.encseq),
                                       suffixarray.suftab,
                                       sizeof suffixarray.suftab,
                                       logger);
            if (arguments->inputlcp)
            {
              if (gt_lcptab_lightweightcheck(
                                      gt_str_get(arguments->esaindexname),
                                      suffixarray.encseq,
                                      suffixarray.readmode,
                                      suffixarray.suftab,
                                      logger,err) != 0)
              {
                haserr = true;
              }
            }
          }
        }
        if (!haserr)
        {
          gt_logger_log(logger, "okay");
        }
      }
      if (!haserr && arguments->inputbwt)
      {
        GtUword bwtdifferentconsecutive = 0, idx;
        GT_UNUSED GtUword longest;

        gt_assert(suffixarray.longest.defined);
        longest = suffixarray.longest.valueunsignedlong;
        if (!haserr && arguments->inputsuf && !arguments->usestream)
        {
          gt_assert(suffixarray.numberofallsortedsuffixes == totallength+1);
          gt_assert(longest < totallength);
          gt_assert(suffixarray.suftab != NULL);
          gt_assert(ESASUFFIXPTRGET(suffixarray.suftab,longest) == 0);
        }
        if (!arguments->usestream)
        {
          for (idx = (GtUword) 1; idx<totallength; idx++)
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

          if (gt_readnextfromstream_GtUchar(&prevcc,&suffixarray.bwttabstream)
              == 1)
          {
            GtUchar cc;
            while (gt_readnextfromstream_GtUchar(&cc,&suffixarray.bwttabstream)
                   == 1)
            {
              if (prevcc != cc || ISSPECIAL(cc))
              {
                bwtdifferentconsecutive++;
              }
            }
          }
        }
        gt_logger_log(logger,"bwtdifferentconsecutive="GT_WU" (%.4f)",
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

static int gt_sfxmap_compressedesa(const char *indexname,GtError *err)
{
  bool haserr = false;
  GtEncseq *encseq;
  GtEncseqLoader *encseq_loader = gt_encseq_loader_new();

  gt_error_check(err);
  encseq = gt_encseq_loader_load(encseq_loader,indexname, err);
  gt_encseq_loader_delete(encseq_loader);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    uint64_t writtenbits;
    uint8_t bitsperentry;
    GtUword numberofentries = 0;
    GtUword totallength = gt_encseq_total_length(encseq);
    FILE *fp;

    fp = gt_fa_fopen_with_suffix(indexname,GT_SUFTABSUFFIX_BYTECOMPRESSED,
                                 "rb",err);
    if (fp == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      if (fread(&writtenbits,sizeof writtenbits,(size_t) 1,fp)
          != (size_t) 1)
      {
        gt_error_set(err,"cannot read %d bytes from %s%s at offset 0",
                        (int) sizeof writtenbits,
                        indexname,GT_SUFTABSUFFIX_BYTECOMPRESSED);
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (fread(&bitsperentry,sizeof bitsperentry,(size_t) 1,fp)
          != (size_t) 1)
      {
        gt_error_set(err,"cannot read %d bytes from %s%s at offset %d",
                        (int) sizeof bitsperentry,
                        indexname,GT_SUFTABSUFFIX_BYTECOMPRESSED,
                        (int) sizeof writtenbits);
        haserr = true;
      }
      gt_assert(writtenbits % bitsperentry == 0);
      numberofentries = (GtUword) (writtenbits/bitsperentry);
      gt_assert(numberofentries == totallength + 1);
    }
    if (!haserr)
    {
      uint64_t readvalue = 0;
      const unsigned int bitspervalue = 64U;
      unsigned int bits2add = (unsigned int) bitsperentry,
                   remainingbits = 0;
      GtUword bitbuffer = 0,
                    countentries = 0,
                    *suftabptr,
                    *suftab = gt_malloc(sizeof *suftab *
                                        (size_t) numberofentries);

      suftabptr = suftab;
      while (true)
      {
        if (remainingbits == 0)
        {
          if (fread(&readvalue,sizeof readvalue,(size_t) 1,fp) != (size_t) 1)
          {
            break;
          }
          remainingbits = bitspervalue;
        } else
        {
          if (remainingbits >= (unsigned int) bits2add)
          {
            if (suftabptr >= suftab + numberofentries)
            {
              break;
            }
            bitbuffer |= ((GtUword) (readvalue) &
                         ((1UL << bits2add) - 1)) << (bitsperentry-bits2add);
            readvalue >>= bits2add;
            remainingbits -= bits2add;
            gt_assert(bitbuffer < (1UL << bitsperentry));
            *suftabptr++ = bitbuffer;
            bitbuffer = 0;
            bits2add = (unsigned int) bitsperentry;
            countentries++;
          } else
          {
            gt_assert(remainingbits < 64U);
            bitbuffer |= (readvalue & ((((uint64_t) 1) << remainingbits) - 1))
                         << (bitsperentry-bits2add);
            readvalue >>= remainingbits;
            bits2add -= remainingbits;
            remainingbits = 0;
          }
        }
      }
      gt_assert(countentries == numberofentries);
      gt_suftab_lightweightcheck(sfx_accesschar_encseq,
                                 sfx_charcount_encseq,
                                 encseq,
                                 GT_READMODE_FORWARD,
                                 totallength,
                                 gt_encseq_alphabetnumofchars(encseq),
                                 suftab,
                                 sizeof *suftab,
                                 NULL);
      gt_free(suftab);
    }
    gt_fa_fclose(fp);
  }
  gt_encseq_delete(encseq);
  return haserr ? -1 : 0;
}

static int gt_sfxmap_compresslcp(const char *indexname,
                                 GtLogger *logger,GtError *err)
{
  Sequentialsuffixarrayreader *ssar = NULL;
  bool haserr = false;

  gt_error_check(err);
  ssar = gt_newSequentialsuffixarrayreaderfromfile(
                                          indexname,
                                          SARR_LCPTAB | SARR_ESQTAB,
                                          true,
                                          logger,
                                          err);
  if (ssar == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    GtUword elems = 0;
    GT_UNUSED GtUword totallength
      = gt_Sequentialsuffixarrayreader_totallength(ssar);
    GtUword maxbranchdepth
      = gt_Sequentialsuffixarrayreader_maxbranchdepth(ssar);
    FILE *fpcompressedlcp = gt_fa_fopen_with_suffix(indexname,
                                     GT_LCPTABSUFFIX_BYTECOMPRESSED,"wb",err);
    unsigned int bitsperentry = gt_determinebitspervalue(maxbranchdepth);
    GtBitbuffer *bitbuffer = NULL;

    gt_assert(ssar != NULL);
    if (fpcompressedlcp == NULL)
    {
       haserr = true;
    } else
    {
      bitbuffer = gt_bitbuffer_new(fpcompressedlcp,bitsperentry);
    }
    if (bitsperentry > 0)
    {
      if (!haserr)
      {
        gt_bitbuffer_next_fixed_bits_value (bitbuffer,0UL);
        while (true)
        {
          GtUword currentlcp;

          SSAR_NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
          gt_bitbuffer_next_fixed_bits_value (bitbuffer,currentlcp);
          elems++;
        }
      }
      if (!haserr)
      {
        gt_assert(elems == totallength);
      }
    }
    gt_freeSequentialsuffixarrayreader(&ssar);
    gt_bitbuffer_delete(bitbuffer);
    if (fpcompressedlcp != NULL)
    {
      gt_fa_fclose(fpcompressedlcp);
    }
  }
  return haserr ? -1 : 0;
}

static int gt_sfxmap_comparelcpvalue(void *info,GtUword lcp,GtError *err)
{
  GtUword currentlcpvalue;
  Sequentialsuffixarrayreader *ssar = (Sequentialsuffixarrayreader *) info;
  bool haserr = false;

  gt_error_check(err);
  do /* fake loop to allow for the use of a break statement */
  {
    SSAR_NEXTSEQUENTIALLCPTABVALUE(currentlcpvalue,ssar);
    if (lcp != currentlcpvalue)
    {
      gt_error_set(err,"lcp="GT_WU" != "GT_WU"=currentlcpvalue",
                   lcp,currentlcpvalue);
      haserr = true;
      break;
    }
    break;
  } while (true);
  return haserr ? -1 : 0;
}

static int gt_sfxmap_pck(const Sfxmapoptions *arguments,GtLogger *logger,
                         GtError *err)
{
  bool haserr = false;
  FMindex *fmindex;
  GtUword totallength = 0;
  unsigned int numofchars = 0;
  GtEncseqMetadata *encseqmetadata = NULL;
  Sequentialsuffixarrayreader *ssar = NULL;

  gt_error_check(err);
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
                                          true,
                                          logger,
                                          err);
      if (ssar == NULL)
      {
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    GtUword idx, pos, numofnonspecials;
    GT_UNUSED GtUword currentsuffix = 0;
    GtSpecialcharinfo specialcharinfo;
    Bwtseqpositioniterator *bspi;

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
        SSAR_NEXTSEQUENTIALSUFTABVALUE(currentsuffix,ssar);
        gt_assert(pos == currentsuffix);
      }
      /*printf(""GT_WU": pos = "GT_WU"\n",idx,pos);*/
    }
    gt_assert(idx == numofnonspecials);
    gt_Bwtseqpositioniterator_delete(bspi);
  }
  if (!haserr)
  {
    GtAlphabet *alphabet;

    alphabet = gt_encseq_metadata_alphabet(encseqmetadata);
    if (alphabet == NULL)
    {
      haserr = true;
    }
    numofchars = gt_alphabet_num_of_chars(alphabet);
  }
  if (!haserr)
  {
    if (arguments->cmplcp && ssar != NULL)
    {
      gt_logger_log(logger,"perform dfs traversal");
      if (gt_fmindex_dfstraverse(fmindex,
                                 numofchars,
                                 totallength,
                                 gt_sfxmap_comparelcpvalue,
                                 (void *) ssar,
                                 err) != 0)
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
                   (unsigned int) BSRS_stream_reader_multi3},
  {"hashfirstcodes","hash first codes of each sequence in the encseq",
                   (unsigned int) BSRS_hashfirstcodes}
};

static int gt_sfxmap_stream_esq(const Sfxmapoptions *arguments,GtError *err)
{
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  bool haserr = false;
  Bitstreamreadmode brsmode = BSRS_stream_single;
  int multiarg = 0;
  GtUword streamesq_size = gt_str_array_size(arguments->streamesq);

  gt_error_check(err);
  if (streamesq_size == 2UL || streamesq_size == 3UL)
  {
    int mode;

    mode = gt_optionargsetsingle(stream_esq_operation,
                                 sizeof stream_esq_operation/
                                 sizeof stream_esq_operation[0],
                                 __func__,
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
             brsmode == BSRS_stream_reader_multi3 ||
             brsmode == BSRS_hashfirstcodes) &&
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
  if (!haserr && brsmode == BSRS_hashfirstcodes)
  {
    if (gt_encseq_mirror(encseq, err) != 0)
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

static void gt_sfxmap_sortmaxdepth_processunsortedrange(void *voiddcov,
                                              GT_UNUSED GtSuffixsortspace *sssp,
                                              GtUword subbucketleft,
                                              GtUword width,
                                              GT_UNUSED GtUword depth)
{
  Checkunsortedrangeinfo *curi = voiddcov;

  gt_checkifprefixesareidentical(__FILE__,
                                 __LINE__,
                                 curi->encseq,
                                 curi->readmode,
                                 curi->sssp,
                                 subbucketleft,
                                 width,
                                 (GtUword) curi->sortmaxdepth);
}

static int gt_sfxmap_performsortmaxdepth(const Sfxmapoptions *arguments,
                                         GtLogger *logger,GtError *err)
{
  bool haserr = false;
  GtEncseqLoader *el = NULL;
  Checkunsortedrangeinfo curi;
  Sfxstrategy sfxstrategy;
  const char *indexname;

  gt_error_check(err);
  indexname = gt_str_get(arguments->esaindexname);
  el = gt_encseq_loader_new();
  curi.encseq = gt_encseq_loader_load(el, indexname, err);
  if (curi.encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_logger_log(logger,"%s(%s,%u)",__func__,
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
    GtUword totallength = gt_encseq_total_length(curi.encseq);

    curi.sssp = gt_suffixsortspace_new(totallength+1,totallength, true, logger);
    gt_suffixsortspace_init_identity(curi.sssp,totallength+1);
    curi.readmode = GT_READMODE_FORWARD;
    curi.sortmaxdepth = arguments->sortmaxdepth;
    gt_sortallsuffixesfromstart(curi.sssp,
                                totallength+1,
                                curi.encseq,
                                curi.readmode,
                                NULL,
                                arguments->sortmaxdepth,
                                &sfxstrategy,
                                gt_sfxmap_sortmaxdepth_processunsortedrange,
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
                           (GtUword) curi.sortmaxdepth);
    gt_suffixsortspace_delete(curi.sssp,false);
  }
  gt_encseq_delete(curi.encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

static int gt_sfxmap_run_diffcover_check(const Sfxmapoptions *arguments,
                                         GtError *err)
{
  int had_err = 0;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  const char *indexname;

  gt_error_check(err);
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

static int showlcpinterval(GT_UNUSED void *data,const Lcpinterval *lcpinterval)
{
  printf("N "GT_WU" "GT_WU" "GT_WU"\n",lcpinterval->offset,
                           lcpinterval->left,
                           lcpinterval->right);
  return 0;
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
    if (arguments->compressedesa)
    {
      if (gt_sfxmap_compressedesa(gt_str_get(arguments->esaindexname),
                                   err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (arguments->compresslcp)
      {
        if (gt_sfxmap_compresslcp(gt_str_get(arguments->esaindexname),
                                  logger,err) != 0)
        {
          haserr = true;
        }
      } else
      {
        if (gt_sfxmap_esa(arguments,logger,err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  if (!haserr && gt_str_length(arguments->pckindexname) > 0)
  {
    if (gt_sfxmap_pck(arguments,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && gt_str_array_size(arguments->streamesq) > 0)
  {
    if (gt_sfxmap_stream_esq(arguments,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->sortmaxdepth > 0)
  {
    if (gt_sfxmap_performsortmaxdepth(arguments,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->diffcovercheck)
  {
    if (gt_sfxmap_run_diffcover_check(arguments, err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && (arguments->enumlcpitvs || arguments->enumlcpitvtree ||
                  arguments->enumlcpitvtreeBU))
  {
    if (arguments->enumlcpitvtreeBU)
    {
      if (gt_runenumlcpvalues_bottomup(gt_str_get(arguments->esaindexname),
                                       logger, err) != 0)
      {
        haserr = true;
      }
    } else
    {
      if (gt_runenumlcpvalues_process(gt_str_get(arguments->esaindexname),
                                      arguments->enumlcpitvs ? showlcpinterval
                                                             : NULL,
                                      NULL,
                                     logger, err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (!haserr && arguments->scanesa > 0)
  {
    if (gt_runscanesa(gt_str_get(arguments->esaindexname),arguments->scanesa,
                      logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->spmitv)
  {
    if (gt_process_spmitv(gt_str_get(arguments->esaindexname),logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && arguments->ownencseq2file)
  {
    if (gt_encseq_check_external_twobitencoding_to_file(
                                 gt_str_get(arguments->esaindexname), err) != 0)
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
