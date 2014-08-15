/*
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
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

#include <errno.h>
#include <stdio.h>
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/bittab_api.h"
#include "core/array_api.h"
#include "core/file.h"
#include "core/fileutils_api.h"
#include "core/undef_api.h"
#include "extended/match.h"
#include "extended/match_blast.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/match_iterator_rep.h"

#define READNUMS 5
#define READVALUES 10

#define GT_MATCHER_BLAST_CANNOTPARSECOLUMN(S)\
        gt_error_set(err, "file %s, line " GT_WU ", column " GT_WU ": %s", \
                     m->pvt->matchfile, m->pvt->curpos, columncount + 1, S)

#define GT_MATCHER_BLAST_CANNOTPARSELINE(S)\
        gt_error_set(err, "file %s, line " GT_WU ": %s", \
                     m->pvt->matchfile, m->pvt->curpos, S)

#define gt_match_iterator_blast_cast(M)\
        gt_match_iterator_cast(gt_match_iterator_blast_class(), M)

typedef struct {
  GtUword curpos;
  FILE *matchfilep;
  GtFile *gtmatchfilep;
  const char *matchfile;
  bool process;
} GtMatchIteratorBlastMembers;

struct GtMatchIteratorBlast {
  const GtMatchIterator parent_instance;
  GtMatchIteratorBlastMembers *pvt;
};

const GtMatchIteratorClass* gt_match_iterator_blast_class(void);

static GtMatchIteratorStatus gt_match_iterator_blast_next(GtMatchIterator *gm,
                                                          GtMatch **match,
                                                          GtError *err)
{
  gt_assert(gm);
  GtUword columncount = 0;
  GtWord storeinteger[READNUMS], tmp;
  double e_value;
  float bitscore, identity;
  bool reverse = false;
  char query_seq[BUFSIZ], db_seq[BUFSIZ], buffer[BUFSIZ];
  int had_err = 0, i = 0, readvalues = 0;
  GtMatchIteratorBlast *m = gt_match_iterator_blast_cast(gm);

  if (m->pvt->matchfilep) {
    if (!m->pvt->process) {
      while (true) {
        if (fgetc(m->pvt->matchfilep) == '#') {
          GT_UNUSED char *l = fgets(buffer, BUFSIZ, m->pvt->matchfilep);
          m->pvt->curpos++;
        } else break;
      }
    }
    if (!m->pvt->process)
      fseek(m->pvt->matchfilep, -1, SEEK_CUR);
    readvalues = fscanf(m->pvt->matchfilep,
                        "%s %s %f " GT_WD " %*d %*d " GT_WD " " GT_WD " " GT_WD
                        " " GT_WD " %lg %f\n", query_seq, db_seq, &identity,
                        &storeinteger[0],
                        &storeinteger[1], &storeinteger[2], &storeinteger[3],
                        &storeinteger[4], &e_value, &bitscore);
    if (readvalues == EOF)
      return GT_MATCHER_STATUS_END;
    if (readvalues != READVALUES)
    {
      GT_MATCHER_BLAST_CANNOTPARSELINE("invalid format");
      had_err = -1;
    }
  } else {
    while (true) {
      while ((buffer[i] = gt_file_xfgetc(m->pvt->gtmatchfilep)) != '\n') {
        if (buffer[i] == EOF)
          return GT_MATCHER_STATUS_END;
        i++;
      }
      buffer[i+1] = '\0';
      if (buffer[0] == '#') {
        m->pvt->curpos++;
        i = 0;
      } else break;
    }
    if ((readvalues = sscanf(buffer, "%s %s %f " GT_WD " %*d %*d " GT_WD " "
                             GT_WD " " GT_WD " " GT_WD " %lg " "%f\n",
                             query_seq, db_seq, &identity,
                             &storeinteger[0],
                             &storeinteger[1], &storeinteger[2],
                             &storeinteger[3], &storeinteger[4], &e_value,
                             &bitscore)) != READVALUES) {
      GT_MATCHER_BLAST_CANNOTPARSELINE("invalid format");
      had_err = -1;
    }
  }

  for (columncount = 0; !had_err && columncount < (GtUword) (READNUMS);
       columncount++) {
    if (storeinteger[columncount] < 0) {
         GT_MATCHER_BLAST_CANNOTPARSECOLUMN("non-negative integer expected");
          had_err = -1;
     }
   }

  if (!had_err) {
    if (storeinteger[1] > storeinteger[2]) {
      tmp = storeinteger[1];
      gt_assert(!reverse);
      reverse = true;
      storeinteger[1] = storeinteger[2];
      storeinteger[2] = tmp;
    }
    if (storeinteger[3] > storeinteger[4]) {
      tmp = storeinteger[3];
      gt_assert(!reverse);
      reverse = true;
      storeinteger[3] = storeinteger[4];
      storeinteger[4] = tmp;
    }
    *match = gt_match_blast_new(query_seq,
                                db_seq,
                                storeinteger[1],
                                storeinteger[2],
                                storeinteger[3],
                                storeinteger[4],
                                e_value,
                                bitscore,
                                storeinteger[0],
                                identity,
                                reverse ? GT_MATCH_REVERSE : GT_MATCH_DIRECT);
    m->pvt->curpos++;
    return GT_MATCHER_STATUS_OK;
  }
  else
  {
    return GT_MATCHER_STATUS_ERROR;
  }
}

GtMatchIterator* gt_match_iterator_blast_file_new(const char *matchfile,
                                                  GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorBlast *mpb;
  mp = gt_match_iterator_create(gt_match_iterator_blast_class());
  mpb = gt_match_iterator_blast_cast(mp);
  mpb->pvt = gt_calloc(1, sizeof (GtMatchIteratorBlastMembers));
  GtFileMode mode;

  if (gt_file_exists(matchfile)) {
    mode = gt_file_mode_determine(matchfile);
    if (mode == GT_FILE_MODE_UNCOMPRESSED) {
      mpb->pvt->matchfilep = fopen(matchfile, "r");
      mpb->pvt->gtmatchfilep = NULL;
      if (!mpb->pvt->matchfilep) {
        gt_error_set(err, "could not open %s", matchfile);
        return NULL;
      }
    } else {
      mpb->pvt->gtmatchfilep = gt_file_open(mode, matchfile, "r", err);
      mpb->pvt->matchfilep = NULL;
      if (!mpb->pvt->gtmatchfilep)
        return NULL;
    }
    mpb->pvt->matchfile = matchfile;
    mpb->pvt->process = false;
    return mp;
  } else {
    gt_error_set(err, "no such file or directory %s", matchfile);
    return NULL;
  }
}

GtMatchIterator* gt_match_iterator_blastalln_process_new(const char *query,
                                                         const char *db_name,
                                                         double evalue,
                                                         bool dust,
                                                         int word_size,
                                                         int gapopen,
                                                         int gapextend,
                                                         int penalty,
                                                         int reward,
                                                         double threshold,
                                                         int num_threads,
                                                         int xdrop_gap_final,
                                                         GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorBlast *mpb;
  char blast_call[BUFSIZ], *env;

  mp = gt_match_iterator_create(gt_match_iterator_blast_class());
  mpb = gt_match_iterator_blast_cast(mp);
  mpb->pvt = gt_calloc(1, sizeof (GtMatchIteratorBlastMembers));
  mpb->pvt->matchfile = query;
  mpb->pvt->process = true;

  env = getenv("GT_BLAST_PATH");
  if (env) {
    sprintf(blast_call, "%s/blastall -p blastn", env);
  } else {
    sprintf(blast_call, "blastall -p blastn");
  }
  if (evalue != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -e %.6e", blast_call, evalue);
  if (dust)
    sprintf(blast_call, "%s -F", blast_call);
  if (word_size != GT_UNDEF_INT)
    sprintf(blast_call, "%s -W %d", blast_call, word_size);
  if (gapopen != GT_UNDEF_INT)
    sprintf(blast_call, "%s -G %d", blast_call, gapopen);
  if (gapextend != GT_UNDEF_INT)
    sprintf(blast_call, "%s -E %d", blast_call, gapextend);
  if (penalty != GT_UNDEF_INT)
    sprintf(blast_call, "%s -q %d", blast_call, penalty);
  if (reward != GT_UNDEF_INT)
    sprintf(blast_call, "%s -r %d", blast_call, reward);
  if (threshold != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -f %.3f", blast_call, threshold);
  if (num_threads != GT_UNDEF_INT)
    sprintf(blast_call, "%s -a %d", blast_call, num_threads);
  if (xdrop_gap_final != GT_UNDEF_INT)
    sprintf(blast_call, "%s -Z %d", blast_call, xdrop_gap_final);
  sprintf(blast_call, "%s -i %s -d %s -m 8", blast_call, query,
           db_name);

  mpb->pvt->matchfilep = popen(blast_call, "r");
  if (!mpb->pvt->matchfilep) {
    gt_error_set(err, "Could not run BLAST process: %s", strerror(errno));
    return NULL;
  }
  mpb->pvt->gtmatchfilep = NULL;
  return mp;
}

GtMatchIterator* gt_match_iterator_blastallp_process_new(const char *query,
                                                         const char *db_name,
                                                         double evalue,
                                                         int word_size,
                                                         int gapopen,
                                                         int gapextend,
                                                         int xdrop_gap_final,
                                                         GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorBlast *mpb;
  char blast_call[BUFSIZ], *env;

  mp = gt_match_iterator_create(gt_match_iterator_blast_class());
  mpb = gt_match_iterator_blast_cast(mp);
  mpb->pvt = gt_calloc(1, sizeof (GtMatchIteratorBlastMembers));
  mpb->pvt->matchfile = query;
  mpb->pvt->process = true;

  env = getenv("GT_BLAST_PATH");
  if (env) {
    sprintf(blast_call, "%s/blastall -p blastp", env);
  } else {
    sprintf(blast_call, "blastall -p blastp");
  }
  if (evalue != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -e %.6e", blast_call, evalue);
  if (word_size != GT_UNDEF_INT)
    sprintf(blast_call, "%s -W %d", blast_call, word_size);
  if (gapopen != GT_UNDEF_INT)
    sprintf(blast_call, "%s -G %d", blast_call, gapopen);
  if (gapextend != GT_UNDEF_INT)
    sprintf(blast_call, "%s -E %d", blast_call, gapextend);
  if (xdrop_gap_final != GT_UNDEF_INT)
    sprintf(blast_call, "%s -Z %d", blast_call, xdrop_gap_final);
  sprintf(blast_call, "%s -i %s -d %s -m 8", blast_call, query,
           db_name);

  mpb->pvt->matchfilep = popen(blast_call, "r");
  if (!mpb->pvt->matchfilep) {
    gt_error_set(err, "Could not run BLAST process: %s", strerror(errno));
    return NULL;
  }
  mpb->pvt->gtmatchfilep = NULL;
  return mp;
}

GtMatchIterator* gt_match_iterator_blastn_process_new(const char *query,
                                                      const char *db_name,
                                                      double evalue,
                                                      bool dust,
                                                      int word_size,
                                                      int gapopen,
                                                      int gapextend,
                                                      int penalty,
                                                      int reward,
                                                      double perc_identity,
                                                      int num_threads,
                                                      double xdrop_gap_final,
                                                      const char *moreblast,
                                                      GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorBlast *mpb;
  char blast_call[BUFSIZ], *env;

  mp = gt_match_iterator_create(gt_match_iterator_blast_class());
  mpb = gt_match_iterator_blast_cast(mp);
  mpb->pvt = gt_calloc(1, sizeof (GtMatchIteratorBlastMembers));
  mpb->pvt->matchfile = query;
  mpb->pvt->process = true;

  env = getenv("GT_BLAST_PATH");
  if (env) {
    sprintf(blast_call, "%s/blastn", env);
  } else {
    sprintf(blast_call, "blastn");
  }
  if (evalue != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -evalue %.6e", blast_call, evalue);
  if (dust)
    sprintf(blast_call, "%s -dust yes", blast_call);
  if (word_size != GT_UNDEF_INT)
    sprintf(blast_call, "%s -word_size %d", blast_call, word_size);
  if (gapopen != GT_UNDEF_INT)
    sprintf(blast_call, "%s -gapopen %d", blast_call, gapopen);
  if (gapextend != GT_UNDEF_INT)
    sprintf(blast_call, "%s -gapextend %d", blast_call, gapextend);
  if (penalty != GT_UNDEF_INT)
    sprintf(blast_call, "%s -penalty %d", blast_call, penalty);
  if (reward != GT_UNDEF_INT)
    sprintf(blast_call, "%s -reward %d", blast_call, reward);
  if (perc_identity != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -perc_identity %.2f", blast_call,
            perc_identity);
  if (num_threads != GT_UNDEF_INT)
    sprintf(blast_call, "%s -num_threads %d", blast_call, num_threads);
  if (xdrop_gap_final != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -xdrop_gap_final %.2f", blast_call,
            xdrop_gap_final);
  if (moreblast)
    sprintf(blast_call, "%s %s", blast_call, moreblast);
  sprintf(blast_call, "%s -query %s -db %s -outfmt 6", blast_call,
          query, db_name);

  mpb->pvt->matchfilep = popen(blast_call, "r");
  if (!mpb->pvt->matchfilep) {
    gt_error_set(err, "Could not run BLAST process: %s", strerror(errno));
    return NULL;
  }
  mpb->pvt->gtmatchfilep = NULL;
  return mp;
}

GtMatchIterator* gt_match_iterator_blastp_process_new(const char *query,
                                                      const char *db_name,
                                                      double evalue,
                                                      int word_size,
                                                      int gapopen,
                                                      int gapextend,
                                                      int num_threads,
                                                      double xdrop_gap_final,
                                                      GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorBlast *mpb;
  char blast_call[BUFSIZ], *env;

  mp = gt_match_iterator_create(gt_match_iterator_blast_class());
  mpb = gt_match_iterator_blast_cast(mp);
  mpb->pvt = gt_calloc(1, sizeof (GtMatchIteratorBlastMembers));
  mpb->pvt->matchfile = query;
  mpb->pvt->process = true;

  env = getenv("GT_BLAST_PATH");
  if (env) {
    sprintf(blast_call, "%s/blastp", env);
  } else {
    sprintf(blast_call, "blastp");
  }
  if (evalue != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -evalue %.6e", blast_call, evalue);
  if (word_size != GT_UNDEF_INT)
    sprintf(blast_call, "%s -word_size %d", blast_call, word_size);
  if (gapopen != GT_UNDEF_INT)
    sprintf(blast_call, "%s -gapopen %d", blast_call, gapopen);
  if (gapextend != GT_UNDEF_INT)
    sprintf(blast_call, "%s -gapextend %d", blast_call, gapextend);
  if (num_threads != GT_UNDEF_INT)
    sprintf(blast_call, "%s -num_threads %d", blast_call, num_threads);
  if (xdrop_gap_final != GT_UNDEF_DOUBLE)
    sprintf(blast_call, "%s -xdrop_gap_final %.2f", blast_call,
             xdrop_gap_final);
  sprintf(blast_call, "%s -query %s -db %s -outfmt 6", blast_call,
           query, db_name);

  mpb->pvt->matchfilep = popen(blast_call, "r");
  if (!mpb->pvt->matchfilep) {
    gt_error_set(err, "Could not run BLAST process: %s", strerror(errno));
    return NULL;
  }
  mpb->pvt->gtmatchfilep = NULL;
  return mp;
}

static void gt_match_iterator_blast_free(GtMatchIterator *mp)
{
  GtMatchIteratorBlast *mpb;
  if (!mp) return;
  mpb = gt_match_iterator_blast_cast(mp);
  if (mpb->pvt->matchfilep) {
    if (mpb->pvt->process) {
      pclose(mpb->pvt->matchfilep);
    } else {
      fclose(mpb->pvt->matchfilep);
    }
    mpb->pvt->matchfilep = NULL;
  }
  if (mpb->pvt->gtmatchfilep != NULL)
    gt_file_delete(mpb->pvt->gtmatchfilep);
  gt_free(mpb->pvt);
}

const GtMatchIteratorClass* gt_match_iterator_blast_class(void)
{
  static const GtMatchIteratorClass *mpc;
  gt_class_alloc_lock_enter();
  if (!mpc) {
    mpc = gt_match_iterator_class_new(sizeof (GtMatchIteratorBlast),
                                      gt_match_iterator_blast_next,
                                      gt_match_iterator_blast_free);
  }
  gt_class_alloc_lock_leave();
  return mpc;
}
