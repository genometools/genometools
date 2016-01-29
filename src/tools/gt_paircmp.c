/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
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

#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_rec.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "core/versionfunc.h"
#include "extended/diagonalbandalign.h"
#include "extended/diagonalbandalign_affinegapcost.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "extended/squarealign.h"
#include "match/test-pairwise.h"
#include "tools/gt_paircmp.h"

typedef struct
{
  GtStr *charlist;
  GtUword len;
} Charlistlen;

typedef struct
{
  GtStrArray *strings,
             *files,
             *fastasequences0,
             *fastasequences1;
  Charlistlen *charlistlen;
  GtStr *text;
  bool fasta;
  bool showedist;
  bool print;
} Cmppairwiseopt;

static void showsimpleoptions(const Cmppairwiseopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    if (!opt->showedist)
      printf("# two strings \"%s\" \"%s\"\n", gt_str_array_get(opt->strings,0),
             gt_str_array_get(opt->strings,1UL));
    return;
  }
  if (gt_str_array_size(opt->files) > 0)
  {
    if (opt->fasta)
    {
      printf("# two files fasta \"%s\" \"%s\"\n",gt_str_array_get(opt->files,1),
           gt_str_array_get(opt->files,2UL));
      return;
    }
    printf("# two files \"%s\" \"%s\"\n", gt_str_array_get(opt->files,0),
           gt_str_array_get(opt->files,1UL));
    return;
  }
  if (opt->charlistlen != NULL)
  {
    printf("# alphalen \"%s\" " GT_WU "\n",
           gt_str_get(opt->charlistlen->charlist),
           opt->charlistlen->len);
    return;
  }
  if (gt_str_length(opt->text) > 0)
  {
    printf("# text \"%s\"\n", gt_str_get(opt->text));
    return;
  }
}

static GtOPrval parse_options(int *parsed_args,
                              Cmppairwiseopt *pw,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstrings,
         *optionfiles,
         *optioncharlistlen,
         *optiontext,
         *optionshowedist,
         *optionprint;
  GtStrArray *charlistlen;
  GtOPrval oprval;

  gt_error_check(err);
  charlistlen = gt_str_array_new();
  pw->strings = gt_str_array_new();
  pw->files = gt_str_array_new();
  pw->text = gt_str_new();
  pw->charlistlen = NULL;
  pw->fastasequences0 = NULL;
  pw->fastasequences1 = NULL;
  pw->showedist = false;
  pw->print = false;
  pw->fasta = false;
  op = gt_option_parser_new("options", "Apply function to pairs of strings.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             pw->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionfiles = gt_option_new_filename_array("ff", "use two files",
                                             pw->files);
  gt_option_parser_add_option(op, optionfiles);

  optioncharlistlen = gt_option_new_string_array("a",
                                             "use character list and length",
                                             charlistlen);
  gt_option_parser_add_option(op, optioncharlistlen);

  optiontext = gt_option_new_string("t", "use text", pw->text, NULL);
  gt_option_parser_add_option(op, optiontext);

  optionshowedist = gt_option_new_bool("e", "output unit edit distance",
                      &pw->showedist, false);
  gt_option_parser_add_option(op, optionshowedist);

  optionprint = gt_option_new_bool("p", "print edist alignment",
                      &pw->print, false);
  gt_option_parser_add_option(op, optionprint);

  gt_option_exclude(optionstrings, optionfiles);
  gt_option_exclude(optionstrings, optioncharlistlen);
  gt_option_exclude(optionstrings, optiontext);
  gt_option_exclude(optionfiles, optioncharlistlen);
  gt_option_exclude(optionfiles, optiontext);
  gt_option_exclude(optioncharlistlen, optiontext);
  gt_option_imply(optionshowedist, optionstrings);
  gt_option_imply(optionprint, optionstrings);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionstrings))
    {
      if (gt_str_array_size(pw->strings) != 2UL)
      {
        gt_error_set(err, "option -ss requires two string arguments");
        oprval = GT_OPTION_PARSER_ERROR;
      }
    } else
    {
      if (gt_option_is_set(optionfiles))
      {
        if (gt_str_array_size(pw->files) != 2UL)
        {
          if (gt_str_array_size(pw->files) == 3UL &&
              !strcmp(gt_str_array_get(pw->files,0),"fasta"))
          {
            pw->fasta = true;
          }
          if (!pw->fasta)
          {
            gt_error_set(err, "option -ff requires two filename arguments or "
                              "keyword fasta and two filename arguments in "
                              "FASTA format");
            oprval = GT_OPTION_PARSER_ERROR;
          }
        }
      } else
      {
        if (gt_option_is_set(optioncharlistlen))
        {
          GtWord readint;
          if (gt_str_array_size(charlistlen) != 2UL)
          {
            gt_error_set(err,
                         "option -a requires charlist and length argument");
            oprval = GT_OPTION_PARSER_ERROR;
          }else
          {
            pw->charlistlen = gt_malloc(sizeof *pw->charlistlen);
            pw->charlistlen->charlist =
              gt_str_ref(gt_str_array_get_str(charlistlen,
                                                                    0));
            if (sscanf(gt_str_array_get(charlistlen,1UL), GT_WD, &readint) != 1
                || readint < 1L)
            {
              gt_error_set(err,
                           "option -a requires charlist and length argument");
              oprval = GT_OPTION_PARSER_ERROR;
            }
            pw->charlistlen->len = (GtUword) readint;
          }
        } else
        {
          if (!gt_option_is_set(optiontext))
          {
            gt_error_set(err,
                         "use exactly one of the options -ss, -ff, -a, -t");
            oprval = GT_OPTION_PARSER_ERROR;
          }
        }
      }
    }
  }
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && *parsed_args != argc)
  {
    gt_error_set(err, "superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  gt_str_array_delete(charlistlen);
  return oprval;
}

static void freesimpleoption(Cmppairwiseopt *cmppairwise)
{
  gt_str_array_delete(cmppairwise->strings);
  gt_str_array_delete(cmppairwise->files);
  gt_str_delete(cmppairwise->text);
  if (cmppairwise->charlistlen != NULL)
  {
    gt_str_delete(cmppairwise->charlistlen->charlist);
    gt_free(cmppairwise->charlistlen);
  }
  if (cmppairwise->fasta)
  {
    gt_str_array_delete(cmppairwise->fastasequences0);
    gt_str_array_delete(cmppairwise->fastasequences1);
  }
}

static GtUword applycheckfunctiontosimpleoptions(
                                  Checkcmppairfuntype checkfunction,
                                  const Cmppairwiseopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    bool forward = true;
    while (true)
    {
      checkfunction(forward,
                    (const GtUchar *) gt_str_array_get(opt->strings,0),
                    (GtUword) strlen(gt_str_array_get(opt->strings,0)),
                    (const GtUchar *) gt_str_array_get(opt->strings,1UL),
                    (GtUword) strlen(gt_str_array_get(opt->strings,1UL)));
      if (!forward)
      {
        break;
      }
      forward = false;
    }
    return 2UL; /* number of testcases */
  }
  if (gt_str_array_size(opt->files) > 0)
  {
    if (opt->fasta)
    {
      GtUword  i, j;
      for (i = 0; i < gt_str_array_size(opt->fastasequences0); i++)
      {
        for (j = 0; j < gt_str_array_size(opt->fastasequences1); j++)
        {
          checkfunction(true,
                    (const GtUchar *) gt_str_array_get(opt->fastasequences0,i),
                    (GtUword) strlen(gt_str_array_get(opt->fastasequences0,i)),
                    (const GtUchar *) gt_str_array_get(opt->fastasequences1,j),
                    (GtUword) strlen(gt_str_array_get(opt->fastasequences1,j)));
        }
      }
    }
    else
    {
      gt_runcheckfunctionontwofiles(checkfunction,
                                    gt_str_array_get(opt->files,0),
                                    gt_str_array_get(opt->files,1UL));
    }
    return 2UL;
  }
  if (opt->charlistlen != NULL)
  {
    return gt_runcheckfunctiononalphalen(checkfunction,
                                      gt_str_get(opt->charlistlen->charlist),
                                      opt->charlistlen->len);
  }
  if (gt_str_length(opt->text) > 0)
  {
    return gt_runcheckfunctionontext(checkfunction, gt_str_get(opt->text));
  }
  gt_assert(false);
  return 0;
}

static int save_fastaentry(const char *seqpart, GT_UNUSED GtUword length,
                           void *data, GT_UNUSED GtError* err)
{
  gt_error_check(err);
  GtStrArray *fasta_sequences = (GtStrArray*) data;
  gt_str_array_add_cstr(fasta_sequences, seqpart);
  return 0;
}

typedef struct {
  Checkcmppairfuntype function;
  const char *name;
} Checkfunctiontabentry;

#define MAKECheckfunctiontabentry(X) {X, #X}

int gt_paircmp(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  Cmppairwiseopt cmppairwise;
  GtOPrval oprval;
  GtFastaReader *reader0 = NULL,
                *reader1 = NULL;

  gt_error_check(err);

  oprval = parse_options(&parsed_args, &cmppairwise, argc, argv, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    gt_assert(parsed_args == argc);
    showsimpleoptions(&cmppairwise);
    if (cmppairwise.showedist)
    {
      GtUword edist, len1, len2;
      GtStr *s1, *s2;

      gt_assert(gt_str_array_size(cmppairwise.strings) >= 2);
      s1 = gt_str_array_get_str(cmppairwise.strings,0);
      s2 = gt_str_array_get_str(cmppairwise.strings,1UL);
      len1 = gt_str_length(s1);
      len2 = gt_str_length(s2);
      edist = gt_computegreedyunitedist((const GtUchar *) gt_str_get(s1),
                                        len1,
                                        (const GtUchar *) gt_str_get(s2),
                                        len2);
      printf(GT_WU " " GT_WU " " GT_WU " " GT_WU "%% errors\n",
             edist, len1,len2,(200 * edist)/(len1+len2));
    }
    else if (cmppairwise.print)
    {
      const GtStr *str0 = gt_str_array_get_str(cmppairwise.strings,0),
                  *str1 = gt_str_array_get_str(cmppairwise.strings,1);

      gt_squarealign_print_edit_alignment((const GtUchar *) gt_str_get(str0),0,
                                          gt_str_length(str0),
                                          (const GtUchar *) gt_str_get(str1),0,
                                          gt_str_length(str1));
    } else
    {
      size_t idx;
      Checkfunctiontabentry checkfunction_tab[] = {
        MAKECheckfunctiontabentry(gt_checkgreedyunitedist),
        MAKECheckfunctiontabentry(gt_linearalign_check),
        MAKECheckfunctiontabentry(gt_linearalign_check_local),
        MAKECheckfunctiontabentry(gt_linearalign_affinegapcost_check),
        MAKECheckfunctiontabentry(gt_linearalign_affinegapcost_check_local),
        MAKECheckfunctiontabentry(gt_diagonalbandalign_check),
        MAKECheckfunctiontabentry(gt_diagonalbandalign_affinegapcost_check)
      };

      if (cmppairwise.fasta)
      {
        gt_assert(gt_str_array_size(cmppairwise.files) == 3);
        cmppairwise.fastasequences0 = gt_str_array_new();
        cmppairwise.fastasequences1 = gt_str_array_new();

        reader0 = gt_fasta_reader_rec_new(gt_str_array_get_str(
                                                        cmppairwise.files,1UL));
        gt_fasta_reader_run(reader0, NULL, save_fastaentry,
                            NULL, cmppairwise.fastasequences0, err);
        reader1 = gt_fasta_reader_rec_new (gt_str_array_get_str(
                                                        cmppairwise.files,2UL));
        gt_fasta_reader_run(reader1, NULL, save_fastaentry,
                            NULL, cmppairwise.fastasequences1, err);
        gt_error_check(err);
      }
      for (idx = 0; idx < sizeof checkfunction_tab/sizeof checkfunction_tab[0];
           idx++)
      {
        GtUword testcases;

        printf("run %s\n",checkfunction_tab[idx].name);
        testcases
          = applycheckfunctiontosimpleoptions(checkfunction_tab[idx].function,
                                              &cmppairwise);
        printf("# number of testcases for %s: " GT_WU "\n",
               checkfunction_tab[idx].name,testcases);
      }
      gt_fasta_reader_delete(reader0);
      gt_fasta_reader_delete(reader1);
    }
  }
  freesimpleoption(&cmppairwise);
  if (oprval == GT_OPTION_PARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == GT_OPTION_PARSER_ERROR)
  {
    return -1;
  }
  return 0;
}
