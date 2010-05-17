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

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include "core/error.h"
#include "core/fa.h"
#include "core/fileutils.h"
#include "core/array.h"
#include "core/str.h"
#include "core/alphabet.h"
#include "sfx-ri-def.h"
#include "fmindex.h"
#include "core/logger.h"
#include "spacedef.h"

#include "fmi-keyval.pr"
#include "fmi-mapspec.pr"

bool gt_fmindexexists(const GtStr *indexname)
{
  if (!gt_file_with_suffix_exists(indexname,FMASCIIFILESUFFIX))
  {
    return false;
  }
  if (!gt_file_with_suffix_exists(indexname,FMDATAFILESUFFIX))
  {
    return false;
  }
  return true;
}

static int scanfmafileviafileptr(Fmindex *fmindex,
                                 GtSpecialcharinfo *specialcharinfo,
                                 bool *storeindexpos,
                                 const GtStr *indexname,
                                 FILE *fpin,
                                 GtLogger *logger,
                                 GtError *err)
{
  bool haserr = false;
  GtArray *riktab;
  unsigned int intstoreindexpos;

  gt_error_check(err);
  riktab = gt_array_new(gt_sizeofReadintkeys());
  SETREADINTKEYS("bwtlength",&fmindex->bwtlength,NULL);
  SETREADINTKEYS("longest",&fmindex->longestsuffixpos,NULL);
  SETREADINTKEYS("storeindexpos",&intstoreindexpos,NULL);
  SETREADINTKEYS("log2blocksize",&fmindex->log2bsize,NULL);
  SETREADINTKEYS("log2markdist",&fmindex->log2markdist,NULL);
  SETREADINTKEYS("specialcharacters",
                 &specialcharinfo->specialcharacters,NULL);
  SETREADINTKEYS("specialranges",
                 &specialcharinfo->specialranges,NULL);
  SETREADINTKEYS("realspecialranges",
                 &specialcharinfo->realspecialranges,NULL);
  SETREADINTKEYS("lengthofspecialprefix",
                 &specialcharinfo->lengthofspecialprefix,NULL);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &specialcharinfo->lengthofspecialsuffix,NULL);
  SETREADINTKEYS("suffixlength",&fmindex->suffixlength,NULL);
  if (!haserr)
  {
    GtStr *currentline;
    unsigned int linenum;

    currentline = gt_str_new();
    for (linenum = 0; gt_str_read_next_line(currentline, fpin) != EOF;
         linenum++)
    {
      if (gt_analyzeuintline(indexname,
                         FMASCIIFILESUFFIX,
                         linenum,
                         gt_str_get(currentline),
                         gt_str_length(currentline),
                         riktab,
                         err) != 0)
      {
        haserr = true;
        break;
      }
      gt_str_reset(currentline);
    }
    gt_str_delete(currentline);
  }
  if (!haserr && gt_allkeysdefined(indexname,FMASCIIFILESUFFIX,riktab,
                                logger,err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (intstoreindexpos == 1U)
    {
      *storeindexpos = true;
    } else
    {
      if (intstoreindexpos == 0)
      {
        *storeindexpos = false;
      } else
      {
        gt_error_set(err,"illegal value in line matching \"storeindexpos=\"");
        haserr = true;
      }
    }
  }
  gt_array_delete(riktab);
  return haserr ? -1 : 0;
}

void gt_freefmindex(Fmindex *fmindex)
{
  if (fmindex->mappedptr != NULL)
  {
    gt_fa_xmunmap(fmindex->mappedptr);
  }
  if (fmindex->bwtformatching != NULL)
  {
    gt_encseq_delete(fmindex->bwtformatching);
    fmindex->bwtformatching = NULL;
  }
  gt_alphabet_delete((GtAlphabet *) fmindex->alphabet);
}

static GtEncseq *mapbwtencoding(const GtStr *indexname,
                                       GtLogger *logger,
                                       GtError *err)
{
  GtEncseqLoader *el;
  GtEncseq *ret;
  gt_error_check(err);

  el = gt_encseq_loader_new();
  gt_encseq_loader_do_not_require_des_tab(el);
  gt_encseq_loader_do_not_require_ssp_tab(el);
  gt_encseq_loader_do_not_require_sds_tab(el);
  gt_encseq_loader_set_logger(el, logger);
  ret = gt_encseq_loader_load(el, (GtStr*) indexname, err);
  gt_encseq_loader_delete(el);
  return ret;
}

int gt_mapfmindex (Fmindex *fmindex,const GtStr *indexname,
                GtLogger *logger,GtError *err)
{
  FILE *fpin = NULL;
  bool haserr = false, storeindexpos = true;
  GtSpecialcharinfo specialcharinfo;

  gt_error_check(err);
  fmindex->mappedptr = NULL;
  fmindex->bwtformatching = NULL;
  fmindex->alphabet = NULL;
  fpin = gt_fa_fopen_with_suffix(indexname,FMASCIIFILESUFFIX,"rb",err);
  if (fpin == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (scanfmafileviafileptr(fmindex,
                              &specialcharinfo,
                              &storeindexpos,
                              indexname,
                              fpin,
                              logger,
                              err) != 0)
    {
      haserr = true;
    }
  }
  gt_fa_xfclose(fpin);
  if (!haserr)
  {
    fmindex->bwtformatching = mapbwtencoding(indexname,logger,err);
    if (fmindex->bwtformatching == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    fmindex->specpos.nextfreeGtPairBwtidx
      = (unsigned long) gt_determinenumberofspecialstostore(&specialcharinfo);
    fmindex->specpos.spaceGtPairBwtidx = NULL;
    fmindex->specpos.allocatedGtPairBwtidx = 0;
    fmindex->alphabet = gt_alphabet_new_from_file(indexname,err);
    if (fmindex->alphabet == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    GtStr *tmpfilename;

    gt_computefmkeyvalues (fmindex,
                        &specialcharinfo,
                        fmindex->bwtlength,
                        fmindex->log2bsize,
                        fmindex->log2markdist,
                        gt_alphabet_num_of_chars(fmindex->alphabet),
                        fmindex->suffixlength,
                        storeindexpos);
    tmpfilename = gt_str_clone(indexname);
    gt_str_append_cstr(tmpfilename,FMDATAFILESUFFIX);
    if (gt_fillfmmapspecstartptr(fmindex,storeindexpos,tmpfilename,err) != 0)
    {
      haserr = true;
    }
    gt_str_delete(tmpfilename);
  }
  if (haserr)
  {
    gt_freefmindex(fmindex);
  }
  return haserr ? -1 : 0;
}
