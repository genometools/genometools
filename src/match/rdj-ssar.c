/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/fileutils_api.h"
#include "match/rdj-ssar.h"

#define GT_READJOINER_CHECK_FILE(SUFFIX)\
  filename = gt_str_clone(indexname);\
  gt_str_append_cstr(filename, SUFFIX);\
  if (!gt_file_exists(gt_str_get(filename)))\
  {\
    gt_error_set(err, "file \"%s\" does not exist "\
                      "or is not readable", gt_str_get(filename));\
    gt_str_delete(filename);\
    return NULL;\
  }\
  gt_str_delete(filename)

Sequentialsuffixarrayreader* gt_readjoiner_ssar_new(GtStr *indexname,
    bool esqtab, bool mapped, GtLogger *logger, GtError *err)
{
  Sequentialsuffixarrayreader *ssar;
  GtStr *filename;

  /*@i1@*/ gt_error_check(err);

  /* check that necessary files are available */
  GT_READJOINER_CHECK_FILE(".suf");
  GT_READJOINER_CHECK_FILE(".lcp");
  GT_READJOINER_CHECK_FILE(".llv");
  GT_READJOINER_CHECK_FILE(".al1");
  GT_READJOINER_CHECK_FILE(".prj");
  if (esqtab)
  {
    GT_READJOINER_CHECK_FILE(".esq");
  }
  /* ".ssp" is not checked as not present for sat eqlen */
  ssar = gt_newSequentialsuffixarrayreaderfromfile(
           gt_str_get(indexname),
           esqtab
             ? SARR_LCPTAB |SARR_SUFTAB | SARR_SSPTAB | SARR_ESQTAB
             : SARR_LCPTAB |SARR_SUFTAB | SARR_SSPTAB,
           mapped
             ? SEQ_mappedboth
             : SEQ_scan,
           logger,
           err);
  if (gt_error_is_set(err))
  {
    if (ssar != NULL)
      gt_freeSequentialsuffixarrayreader(&ssar);
    return NULL;
  }
  else
  {
    return ssar;
  }
}

void gt_readjoiner_ssar_delete(Sequentialsuffixarrayreader *ssar)
{
  if (ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&ssar);
}
