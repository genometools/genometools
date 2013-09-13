/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/file.h"
#include "core/undef_api.h"
#include "core/xansi_api.h"
#include "match/reads_libraries_table.h"

typedef struct {
  GtUword first_seqnum;
  GtUword insertlength;
  GtUword stdev;
} GtReadsLibrary;

struct GtReadsLibrariesTable {
  GtUword noflibraries;
  GtUword firstunpaired;
  GtReadsLibrary *library;
  GtUword nextfreelibrary;
};

GtReadsLibrariesTable* gt_reads_libraries_table_new(GtUword noflibraries)
{
  GtReadsLibrariesTable *rlt;
  gt_assert(noflibraries > 0);
  rlt = gt_malloc(sizeof (*rlt));
  rlt->library = gt_malloc(sizeof (*rlt->library) * noflibraries);
  rlt->noflibraries = noflibraries;
  rlt->nextfreelibrary = 0;
  rlt->firstunpaired = GT_UNDEF_UWORD;
  return rlt;
}

void gt_reads_libraries_table_delete(GtReadsLibrariesTable *rlt)
{
  if (rlt != NULL)
  {
    gt_free(rlt->library);
    gt_free(rlt);
  }
}

void gt_reads_libraries_table_add(GtReadsLibrariesTable *rlt,
    GtUword first_seqnum, GtUword insertlength,
    GtUword stdev, bool paired)
{
  gt_assert(rlt != NULL);
  gt_assert(rlt->nextfreelibrary < rlt->noflibraries);
  if (!paired)
  {
    if (rlt->firstunpaired == GT_UNDEF_UWORD)
    {
      rlt->firstunpaired = first_seqnum;
    }
    else
    {
      gt_assert(rlt->firstunpaired < first_seqnum);
    }
  }
  else
  {
    gt_assert(rlt->firstunpaired == GT_UNDEF_UWORD);
  }
  rlt->library[rlt->nextfreelibrary].first_seqnum = first_seqnum;
  rlt->library[rlt->nextfreelibrary].insertlength = insertlength;
  rlt->library[rlt->nextfreelibrary].stdev = stdev;
  rlt->nextfreelibrary++;
}

void gt_reads_libraries_table_get_library(GtReadsLibrariesTable *rlt,
    GtUword libnum, GtUword *first_seqnum,
    GtUword *insertlength, GtUword *stdev)
{
  gt_assert(rlt != NULL);
  gt_assert(libnum < rlt->noflibraries);
  if (first_seqnum != NULL)
    *first_seqnum = rlt->library[libnum].first_seqnum;
  if (insertlength != NULL)
    *insertlength = rlt->library[libnum].insertlength;
  if (stdev != NULL)
    *stdev = rlt->library[libnum].stdev;
}

GtUword gt_reads_libraries_table_noflibraries(GtReadsLibrariesTable *rlt)
{
  gt_assert(rlt != NULL);
  return rlt->noflibraries;
}

GtUword gt_reads_libraries_table_firstunpaired(GtReadsLibrariesTable *rlt)
{
  gt_assert(rlt != NULL);
  return rlt->firstunpaired;
}

void gt_reads_libraries_table_save(GtReadsLibrariesTable *rlt, FILE *rlt_fp)
{
  gt_assert(rlt != NULL);
  gt_assert(rlt_fp != NULL);
  gt_xfwrite(&(rlt->noflibraries), sizeof (rlt->noflibraries), (size_t)1,
      rlt_fp);
  gt_xfwrite(&(rlt->firstunpaired), sizeof (rlt->firstunpaired), (size_t)1,
      rlt_fp);
  gt_assert(rlt->noflibraries > 0);
  gt_xfwrite(&(rlt->library), sizeof (*(rlt->library)),
      (size_t)rlt->noflibraries, rlt_fp);
}

GtReadsLibrariesTable* gt_reads_libraries_table_load(FILE *rlt_fp,
                                                     GtError *err)
{
  GtReadsLibrariesTable *rlt;
  size_t count;
  gt_assert(rlt_fp != NULL);
  gt_error_check(err);
  rlt = gt_malloc(sizeof (*rlt));
  count = fread(&(rlt->noflibraries), sizeof (rlt->noflibraries), (size_t)1,
      rlt_fp);
  if (count != (size_t)1 || rlt->noflibraries == 0)
  {
    gt_error_set(err, "libraries table: error by reading number of libraries");
    gt_free(rlt);
    return NULL;
  }
  else
  {
    rlt->library = gt_malloc(sizeof (*(rlt->library)) * rlt->noflibraries);
    count = fread(rlt->library, sizeof (*(rlt->library)),
        (size_t)rlt->noflibraries, rlt_fp);
    if (count != (size_t)rlt->noflibraries)
    {
      gt_error_set(err, "library table: "GT_WU" libraries expected, "GT_WU
                   " found", rlt->noflibraries, (GtUword)count);
      gt_free(rlt->library);
      gt_free(rlt);
      return NULL;
    }
  }
  return rlt;
}
