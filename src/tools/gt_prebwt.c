/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/alphabet_api.h"
#include "core/encseq_metadata.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/tool_api.h"
#include "match/pckbucket.h"
#include "match/eis-voiditf.h"

typedef struct
{
  GtStr *indexname;
  unsigned int maxdepth;
} Prebwtoptions;

static void *gt_prebwt_arguments_new(void)
{
  return gt_malloc(sizeof (Prebwtoptions));
}

static void gt_prebwt_arguments_delete(void *tool_arguments)
{
  Prebwtoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_free(arguments);
}

static GtOptionParser* gt_prebwt_option_parser_new(void *tool_arguments)
{
  Prebwtoptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionpck;

  gt_assert(arguments != NULL);
  arguments->indexname = gt_str_new();
  op = gt_option_parser_new("[options] -pck indexname",
                            "Precompute bwt-bounds for some prefix length.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optionpck = gt_option_new_string("pck","Specify index (packed index)",
                             arguments->indexname, NULL);
  gt_option_parser_add_option(op, optionpck);
  gt_option_is_mandatory(optionpck);

  option = gt_option_new_uint_min("maxdepth",
                                  "specify maximum depth (value > 0)",
                                  &arguments->maxdepth,0,1U);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_prebwt_runner(GT_UNUSED int argc,
                            GT_UNUSED const char **argv,
                            GT_UNUSED int parsed_args,
                            void *tool_arguments, GtError *err)
{
  unsigned long totallength = 0;
  FMindex *fmindex = NULL;
  bool haserr = false;
  GtEncseqMetadata *emd = NULL;
  Prebwtoptions *prebwtoptions = (Prebwtoptions *) tool_arguments;
  GtAlphabet *alphabet;
  const char *indexname = gt_str_get(prebwtoptions->indexname);

  emd = gt_encseq_metadata_new(indexname, err);
  if (emd == NULL) {
    haserr = true;
  }
  if (!haserr) {
    alphabet = gt_encseq_metadata_alphabet(emd);
    if (alphabet == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    fmindex = gt_loadvoidBWTSeqForSA(indexname,false, err);
    if (fmindex == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    unsigned int numofchars = gt_alphabet_num_of_chars(alphabet);
    Pckbuckettable *pckbt;
    totallength = gt_voidpackedindex_totallength_get(fmindex);

    pckbt = gt_pckbuckettable_new((const void *) fmindex,
                                  numofchars,
                                  totallength,
                                  prebwtoptions->maxdepth);
    if (gt_pckbuckettable_2file(indexname,pckbt,err) != 0)
    {
      haserr = true;
    }
    gt_pckbuckettable_delete(pckbt);
  }
  gt_encseq_metadata_delete(emd);
  if (fmindex != NULL)
  {
    gt_deletevoidBWTSeq(fmindex);
  }
  return haserr ? -1 : 0;
}

GtTool* gt_prebwt(void)
{
  return gt_tool_new(gt_prebwt_arguments_new,
                  gt_prebwt_arguments_delete,
                  gt_prebwt_option_parser_new,
                  NULL,
                  gt_prebwt_runner);
}
