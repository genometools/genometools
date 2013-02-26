/*
  Copyright (c) 2010      Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/bittab_api.h"
#include "core/array_api.h"
#include "core/file.h"
#include "core/fileutils_api.h"
#include "extended/match.h"
#include "extended/match_open.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_open.h"
#include "extended/match_iterator_rep.h"

#define READNUMS 5

#define GT_MATCHER_OPEN_CANNOTPARSECOLUMN(S)\
        gt_error_set(err,"file %s, line %lu, column %lu: %s", \
                     mpi->pvt->matchfile, mpi->pvt->curpos, columncount+1, S)

#define GT_MATCHER_OPEN_CANNOTPARSELINE(S)\
        gt_error_set(err,"file %s, line %lu: %s", \
                     mpi->pvt->matchfile, mpi->pvt->curpos, S)

#define gt_match_iterator_open_cast(M)\
        gt_match_iterator_cast(gt_match_iterator_open_class(), M)

typedef struct {
  unsigned long curpos;
  FILE *matchfilep;
  GtFile *gtmatchfilep;
  const char *matchfile;
} GtMatchIteratorOpenMembers;

struct GtMatchIteratorOpen {
  const GtMatchIterator parent_instance;
  GtMatchIteratorOpenMembers *pvt;
};

const GtMatchIteratorClass* gt_match_iterator_open_class(void);

static GtMatchIteratorStatus gt_match_iterator_open_next(GtMatchIterator *gmpi,
                                                         GtMatch **match,
                                                         GtError *err)
{
  unsigned long columncount = 0;
  int readnums;
  long storeinteger[READNUMS];
  int had_err = 0, i = 0;
  char buffer[BUFSIZ], seqid1[BUFSIZ], seqid2[BUFSIZ], matchtype;
  GtMatchIteratorOpen *mpi = gt_match_iterator_open_cast(gmpi);
  gt_assert(mpi);

  if (mpi->pvt->matchfilep) {
    while (true) {
      if (fgetc(mpi->pvt->matchfilep) == '#') {
        GT_UNUSED char *l = fgets(buffer, BUFSIZ, mpi->pvt->matchfilep);
        mpi->pvt->curpos++;
      } else break;
    }
    fseek(mpi->pvt->matchfilep, -1, SEEK_CUR);
    readnums = fscanf(mpi->pvt->matchfilep," %ld %s %ld %c %ld %s %ld %*d %*e "
                      "%ld %*f\n",
                      &storeinteger[0],
                      seqid1,
                      &storeinteger[1],
                      &matchtype,
                      &storeinteger[2],
                      seqid2,
                      &storeinteger[3],
                      &storeinteger[4]);
    if (readnums == EOF)
      return GT_MATCHER_STATUS_END;
    if (readnums != READNUMS + 3)
    {
      GT_MATCHER_OPEN_CANNOTPARSELINE("invalid format");
      had_err = -1;
    }
  } else {
    while (true) {
      while ((buffer[i] = gt_file_xfgetc(mpi->pvt->gtmatchfilep)) != '\n') {
        if (buffer[i] == EOF)
          return GT_MATCHER_STATUS_END;
        i++;
      }
      buffer[i+1] = '\0';
      if (buffer[0] == '#') {
        mpi->pvt->curpos++;
        i = 0;
      } else break;
    }
    if (sscanf(buffer," %ld %s %ld %*c %ld %s %ld %*d %*e %ld %*f\n",
               &storeinteger[0],
               seqid1,
               &storeinteger[1],
               &storeinteger[2],
               seqid2,
               &storeinteger[3],
               &storeinteger[4])
               != READNUMS + 2) {
      GT_MATCHER_OPEN_CANNOTPARSELINE("invalid format");
      had_err = -1;
    }
  }

  for (columncount = 0; columncount < (unsigned long) (READNUMS);
       columncount++) {
    if (storeinteger[columncount] < 0) {
         GT_MATCHER_OPEN_CANNOTPARSECOLUMN("non-negative integer expected");
          had_err = -1;
     }
   }

  if (!had_err) {
    *match = gt_match_open_new(seqid1,
                               seqid2,
                               storeinteger[1],
                               storeinteger[1] + storeinteger[0] - 1,
                               storeinteger[3],
                               storeinteger[3] + storeinteger[2] - 1,
                               storeinteger[4],
                               matchtype == 'D' ? GT_MATCH_DIRECT :
                                                  GT_MATCH_REVERSE);
    mpi->pvt->curpos++;
    return GT_MATCHER_STATUS_OK;
  }
  else {
    return GT_MATCHER_STATUS_ERROR;
  }
}

GtMatchIterator* gt_match_iterator_open_new(const char *matchfile, GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorOpen *mpo;
  mp = gt_match_iterator_create(gt_match_iterator_open_class());
  mpo = (GtMatchIteratorOpen*) mp;
  mpo->pvt = gt_calloc(1, sizeof (GtMatchIteratorOpenMembers));
  GtFileMode mode;
  if (gt_file_exists(matchfile)) {
    mode = gt_file_mode_determine(matchfile);
    if (mode == GT_FILE_MODE_UNCOMPRESSED) {
      mpo->pvt->matchfilep = fopen(matchfile, "r");
      mpo->pvt->gtmatchfilep = NULL;
      if (!mpo->pvt->matchfilep) {
        gt_error_set(err, "Could not open %s", matchfile);
        return NULL;
      }
    } else {
      mpo->pvt->gtmatchfilep = gt_file_open(mode, matchfile, "r", err);
      mpo->pvt->matchfilep = NULL;
      if (!mpo->pvt->gtmatchfilep)
        return NULL;
    }
    mpo->pvt->matchfile = matchfile;
    return mp;
  } else {
    gt_error_set(err, "No such file or directory %s", matchfile);
    return NULL;
  }
}

static void gt_match_iterator_open_free(GtMatchIterator *mp)
{
  GtMatchIteratorOpen *mpo;
  if (!mp) return;
  mpo = gt_match_iterator_open_cast(mp);
  if (mpo->pvt->matchfilep != NULL)
    fclose(mpo->pvt->matchfilep);
  if (mpo->pvt->gtmatchfilep != NULL)
    gt_file_delete(mpo->pvt->gtmatchfilep);
  gt_free(mpo->pvt);
}

const GtMatchIteratorClass* gt_match_iterator_open_class(void)
{
  static const GtMatchIteratorClass *mpc;
  gt_class_alloc_lock_enter();
  if (!mpc) {
    mpc = gt_match_iterator_class_new(sizeof (GtMatchIteratorOpen),
                               gt_match_iterator_open_next,
                               gt_match_iterator_open_free);
  }
  gt_class_alloc_lock_leave();
  return mpc;
}
