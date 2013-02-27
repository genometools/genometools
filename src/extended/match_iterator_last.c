/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include "core/class_alloc_lock.h"
#include "core/cstr.h"
#include "core/cstr_array.h"
#include "core/encseq.h"
#include "core/fasta.h"
#include "core/fileutils.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/md5_fingerprint_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/match.h"
#include "extended/match_last_api.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_last.h"
#include "extended/match_iterator_rep.h"

const GtMatchIteratorClass* gt_match_iterator_last_class(void);

#define gt_match_iterator_last_cast(M)\
        gt_match_iterator_cast(gt_match_iterator_last_class(), M)

typedef struct {
  int match_score,
      mismatch_cost,
      gap_open_cost,
      gap_ext_cost,
      xdrop,
      ydrop,
      zdrop,
      mscoregapped,
      mscoregapless,
      k;
} GtLastOptions;

typedef struct {
  GtEncseq *es1, *es2;
  GtLastOptions op;
  char *tmpdir;
  GtStr *idxfilename,
        *queryfilename,
        *matchfilename,
        *linebuf;
  GtFile *matchfile;
  GtHashmap *desc_to_seqno;
} GtMatchIteratorLastMembers;

struct GtMatchIteratorLast {
  const GtMatchIterator parent_instance;
  GtMatchIteratorLastMembers *pvt;
};

static int get_index_parameterization(GtMatchIteratorLast *mil, GtStr *out,
                                GT_UNUSED GtError *err)
{
  int had_err = 0;
  gt_assert(mil);
  gt_error_check(err);
  gt_str_reset(out);
  gt_str_append_cstr(out, mil->pvt->tmpdir);
  gt_str_append_cstr(out, "/");
  gt_str_append_cstr(out, gt_encseq_indexname(mil->pvt->es1));
  return had_err;
}

static int get_run_parameterization(GtMatchIteratorLast *mil, GtStr *out,
                                    GtStr *lastbinary,
                                    const char *indexname, const char *qryname,
                                    GtError *err)
{
  int had_err = 0;
  const char *env;
  gt_assert(mil);
  gt_error_check(err);
  gt_str_reset(out);
  gt_str_reset(lastbinary);

  env = getenv("GT_LAST_PATH");
  if (env) {
    gt_str_append_cstr(out, env);
    gt_str_append_cstr(out, "/lastal");
    if (!gt_file_exists(gt_str_get(out))) {
      gt_error_set(err, "cannot find LAST executable at %s", gt_str_get(out));
      had_err = -1;
    }
  } else {
    gt_str_append_cstr(out, "lastal");
  }
  gt_str_append_str(lastbinary, out);

  if (!had_err) {
    if (mil->pvt->op.match_score != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -r ");
      gt_str_append_int(out, mil->pvt->op.match_score);
    }
    if (mil->pvt->op.mismatch_cost != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -q ");
      gt_str_append_int(out, mil->pvt->op.mismatch_cost);
    }
    if (mil->pvt->op.gap_open_cost != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -a ");
      gt_str_append_int(out, mil->pvt->op.gap_open_cost);
    }
    if (mil->pvt->op.gap_ext_cost != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -b ");
      gt_str_append_int(out, mil->pvt->op.gap_ext_cost);
    }
    if (mil->pvt->op.xdrop != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -x ");
      gt_str_append_int(out, mil->pvt->op.xdrop);
    }
    if (mil->pvt->op.ydrop != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -y ");
      gt_str_append_int(out, mil->pvt->op.ydrop);
    }
    if (mil->pvt->op.zdrop != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -z ");
      gt_str_append_int(out, mil->pvt->op.zdrop);
    }
    if (mil->pvt->op.k != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -k ");
      gt_str_append_int(out, mil->pvt->op.k);
    }
    if (mil->pvt->op.mscoregapped != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -e ");
      gt_str_append_int(out, mil->pvt->op.mscoregapped);
    }
    if (mil->pvt->op.mscoregapless != GT_UNDEF_INT) {
      gt_str_append_cstr(out, " -d ");
      gt_str_append_int(out, mil->pvt->op.mscoregapless);
    }
  }
  gt_str_append_cstr(out," -f 0 ");
  gt_str_append_cstr(out, indexname);
  gt_str_append_cstr(out," ");
  gt_str_append_cstr(out, qryname);
  return had_err;
}

static int last_prepare_fasta_seqs(const char *filename, GtEncseq *encseq,
                                   GtHashmap *desc_to_seqno,
                                   GT_UNUSED GtError *err)
{
  char *seq;
  const char *desc;
  int had_err = 0;
  GtFile *fp;
  GtStr *header = NULL;
  unsigned long desclen;
  gt_assert(filename && encseq);
  gt_error_check(err);

  if (!(fp = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, filename, "w", err)))
    had_err = -1;
  header = gt_str_new();
  if (!had_err) {
    unsigned long i, length, startpos, endpos;
    for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++) {
      char *desccpy;
      desc = gt_encseq_description(encseq, &desclen, i);
      desccpy = gt_calloc(desclen+1, sizeof (char));
      strncpy(desccpy, desc, desclen);
      startpos = gt_encseq_seqstartpos(encseq, i);
      length = gt_encseq_seqlength(encseq, i);
      gt_hashmap_add(desc_to_seqno, desccpy, (void*) i);
      endpos = startpos + length - 1;
      seq = gt_malloc(sizeof (char) * length);
      gt_str_reset(header);
      gt_str_append_cstr_nt(header, desc, desclen);
      gt_encseq_extract_decoded(encseq, seq, startpos, endpos);
      gt_fasta_show_entry(gt_str_get(header), seq, length, 80, fp);
      gt_free(seq);
    }
    gt_file_delete(fp);
  }
  gt_str_delete(header);

  return had_err;
}

static bool last_parse_comment_line(GtMatchIteratorLast *mil)
{
  char c;
  gt_assert(mil && mil->pvt->matchfile);
  c = gt_file_xfgetc(mil->pvt->matchfile);
  if (c == '#') {
    while (gt_file_xfgetc(mil->pvt->matchfile) != '\n');
  } else {
    gt_file_unget_char(mil->pvt->matchfile, c);
    return false;
  }
  return true;
}

static int last_parse_match(GtMatchIteratorLast *mil, GtMatch **match,
                            GtError *err)
{
  int had_err = 0;
  unsigned long score = GT_UNDEF_ULONG,
                start1 = GT_UNDEF_ULONG,
                start2 = GT_UNDEF_ULONG,
                mlength1 = GT_UNDEF_ULONG,
                mlength2 = GT_UNDEF_ULONG,
                slength1 = GT_UNDEF_ULONG,
                slength2 = GT_UNDEF_ULONG,
                seqno1 = GT_UNDEF_ULONG,
                seqno2 = GT_UNDEF_ULONG;
  char strand1 = GT_UNDEF_CHAR,
       strand2 = GT_UNDEF_CHAR,
       seqid1[BUFSIZ],
       seqid2[BUFSIZ];

  gt_assert(mil && mil->pvt->matchfile);
  gt_error_check(err);

  gt_str_reset(mil->pvt->linebuf);
  had_err = gt_str_read_next_line_generic(mil->pvt->linebuf,
                                          mil->pvt->matchfile);
  if (!had_err) {
    if (11 != sscanf(gt_str_get(mil->pvt->linebuf), "%lu %s %lu %lu %c %lu "
                                     "%s %lu %lu %c %lu",
                                     &score,
                                     seqid1, &start1, &mlength1,
                                     &strand1, &slength1,
                                     seqid2, &start2, &mlength2,
                                     &strand2, &slength2)) {
      gt_error_set(err, "error parsing LAST output: "
                        "could not parse line '%s'",
                   gt_str_get(mil->pvt->linebuf));
      had_err = -1;
    }
    if (!had_err) {
      seqno1 = (unsigned long) gt_hashmap_get(mil->pvt->desc_to_seqno, seqid1);
      gt_assert(seqno1 != GT_UNDEF_ULONG);
      seqno2 = (unsigned long) gt_hashmap_get(mil->pvt->desc_to_seqno, seqid2);
      gt_assert(seqno2 != GT_UNDEF_ULONG);
      *match = gt_match_last_new(seqid1, seqid2, score,
                                 seqno1, seqno2,
                                 start1, start2,
                                 start1 + mlength1 - 1, start2 + mlength2 - 1,
                                 strand1 == strand2 ? GT_MATCH_DIRECT :
                                                      GT_MATCH_REVERSE);
    }
  }
  return had_err;
}

static GtMatchIteratorStatus gt_match_iterator_last_next(GtMatchIterator *gmpi,
                                                         GtMatch **match,
                                                         GtError *err)
{
  int had_err = 0, rval;
  GtMatchIteratorLast *mil;
  GtStr *tmp = NULL,
        *lastbinary = NULL,
        *matchfilename = NULL;
  gt_assert(gmpi && match);
  gt_error_check(err);
  mil = gt_match_iterator_last_cast(gmpi);

  if (!mil->pvt->queryfilename) {
    GtStr *fn = gt_str_clone(mil->pvt->idxfilename);
    gt_str_append_cstr(fn, ".qry");
    if (!gt_file_exists(gt_str_get(fn))) {
      last_prepare_fasta_seqs(gt_str_get(fn), mil->pvt->es2,
                              mil->pvt->desc_to_seqno, err);
    }
    mil->pvt->queryfilename = gt_str_ref(fn);
    gt_str_delete(fn);
  }

  if (!mil->pvt->matchfile) {
    int fdout;
    pid_t pid;
    gt_assert(mil->pvt->idxfilename
                && gt_str_length(mil->pvt->idxfilename) > 0);

    tmp = gt_str_new();
    lastbinary = gt_str_new();
    matchfilename = gt_str_new();
    had_err = get_run_parameterization(mil, tmp, lastbinary,
                                       gt_str_get(mil->pvt->idxfilename),
                                       gt_str_get(mil->pvt->queryfilename),
                                       err);
    if (!had_err) {
      char *matchfilehash;
      gt_assert(gt_str_length(tmp) > 0);
      matchfilehash = gt_md5_fingerprint(gt_str_get(tmp), gt_str_length(tmp));
      gt_assert(matchfilehash && strlen(matchfilehash) > 0);
      gt_str_append_cstr(matchfilename, mil->pvt->tmpdir);
      gt_str_append_cstr(matchfilename, "/");
      gt_str_append_cstr(matchfilename, matchfilehash);
      gt_str_append_cstr(matchfilename, ".match");
    }

    if (!had_err && !gt_file_exists(gt_str_get(matchfilename))) {
      char **args = NULL;
      args = gt_cstr_split(gt_str_get(tmp), ' ');
      fdout = open(gt_str_get(matchfilename), O_WRONLY | O_CREAT, 0600);
      switch ((pid = fork())) {
        case -1:
               had_err = -1;
               break;
        case 0:
          dup2(fdout, STDOUT_FILENO);
          close(fdout);
          execvp(gt_str_get(lastbinary), args);
        default:
        {
          int status;
          while (-1 == waitpid(pid, &status, 0));
        }
      }
      if (had_err == -1)
        gt_error_set(err, "error forking the LAST process");
      gt_cstr_array_delete(args);
    }
    if (!had_err) {
      mil->pvt->matchfile = gt_file_open(GT_FILE_MODE_UNCOMPRESSED,
                                         gt_str_get(matchfilename), "r", err);
      if (!mil->pvt->matchfile) {
        had_err = -1;
      }
    }
    gt_str_delete(lastbinary);
    gt_str_delete(tmp);

  }
  gt_str_delete(matchfilename);

  if (!had_err) {
    /* result file should now be open */
    gt_assert(mil->pvt->matchfile);
    /* skip all comment lines */
    while (last_parse_comment_line(mil));
    /* parse next match */
    if ((rval = last_parse_match(mil, match, err)) == EOF) {
      *match = NULL;
      return GT_MATCHER_STATUS_END;
    }
  }
  /* return error if required */
  if (had_err) {
    *match = NULL;
    return GT_MATCHER_STATUS_ERROR;
  }
  return GT_MATCHER_STATUS_OK;
}

static int last_prepare_indices(GtMatchIteratorLast *mil,
                                GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtStr *param = gt_str_new(),
        *idxfilename = gt_str_new();
  gt_assert(mil);
  gt_error_check(err);

  had_err = get_index_parameterization(mil, param, err);
  if (!had_err) {
    GtStr *cmdline = NULL,
          *fn = NULL;
    char *seq, *hash;

    seq = gt_malloc(sizeof (char) * gt_encseq_total_length(mil->pvt->es1)+1);
    seq[gt_encseq_total_length(mil->pvt->es1)] = '\0';
    gt_encseq_extract_decoded(mil->pvt->es1,
                              seq, 0,
                              gt_encseq_total_length(mil->pvt->es1)-1);
    hash = gt_md5_fingerprint(seq, gt_encseq_total_length(mil->pvt->es1));
    gt_free(seq);

    /* maybe hash the sequence file! */
    gt_str_append_cstr(idxfilename, mil->pvt->tmpdir);
    gt_str_append_cstr(idxfilename, "/");
    gt_str_append_cstr(idxfilename, hash);
    mil->pvt->idxfilename = gt_str_ref(idxfilename);
    fn = gt_str_clone(mil->pvt->idxfilename);
    gt_str_append_cstr(fn, "0.suf");
    if (!gt_file_exists(gt_str_get(fn))) {
      const char *env;
      had_err = last_prepare_fasta_seqs(gt_str_get(mil->pvt->idxfilename),
                                        mil->pvt->es1, mil->pvt->desc_to_seqno,
                                        err);
      cmdline = gt_str_new();
      env = getenv("GT_LAST_PATH");
      if (env) {
        gt_str_append_cstr(cmdline, env);
        gt_str_append_cstr(cmdline, "/lastdb ");
        if (!gt_file_exists(gt_str_get(cmdline))) {
          gt_error_set(err, "cannot find LASTDB executable at %s",
                       gt_str_get(cmdline));
          had_err = -1;
        }
      } else {
        gt_str_append_cstr(cmdline, "lastdb ");
      }
      gt_str_append_str(cmdline, idxfilename);
      gt_str_append_cstr(cmdline, " ");
      gt_str_append_str(cmdline, idxfilename);
      had_err = system(gt_str_get(cmdline));
      if (had_err)
        gt_error_set(err, "error forking the LAST process");
      gt_str_delete(cmdline);
    }
    gt_str_delete(fn);
    gt_free(hash);
  }
  gt_str_delete(param);
  gt_str_delete(idxfilename);
  return had_err;
}

GtMatchIterator* gt_match_iterator_last_new(GtEncseq *es1, GtEncseq *es2,
                                            int match_score,
                                            int mismatch_cost,
                                            int gap_open_cost,
                                            int gap_ext_cost,
                                            int xdrop,
                                            int ydrop,
                                            int zdrop,
                                            int k,
                                            int mscoregapped,
                                            int mscoregapless,
                                            GtError *err)
{
  GtMatchIterator *mp;
  GtMatchIteratorLast *mil;
  int had_err = 0;
  gt_assert(es1 && es2);
  gt_error_check(err);

  mp = gt_match_iterator_create(gt_match_iterator_last_class());
  mil = (GtMatchIteratorLast*) mp;
  mil->pvt = gt_calloc(1, sizeof (GtMatchIteratorLastMembers));
  mil->pvt->es1 = gt_encseq_ref(es1);
  mil->pvt->es2 = gt_encseq_ref(es2);
  mil->pvt->op.match_score = match_score;
  mil->pvt->op.mismatch_cost = mismatch_cost;
  mil->pvt->op.gap_open_cost = gap_open_cost;
  mil->pvt->op.gap_ext_cost = gap_ext_cost;
  mil->pvt->op.xdrop = xdrop;
  mil->pvt->op.ydrop = ydrop;
  mil->pvt->op.zdrop = zdrop;
  mil->pvt->op.mscoregapped = mscoregapped;
  mil->pvt->op.mscoregapless = mscoregapless;
  mil->pvt->op.k = k;
  mil->pvt->tmpdir = getenv("TMPDIR");
  if (!mil->pvt->tmpdir)
    mil->pvt->tmpdir = "/tmp";
  mil->pvt->linebuf = gt_str_new();
  mil->pvt->desc_to_seqno = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);

  had_err = last_prepare_indices(mil, err);
  if (had_err) {
    gt_match_iterator_delete(mp);
    return NULL;
  }
  gt_assert(mil->pvt->idxfilename);
  return mp;
}

static void gt_match_iterator_last_free(GtMatchIterator *mp)
{
  GtMatchIteratorLast *mil;
  if (!mp) return;
  mil = gt_match_iterator_last_cast(mp);
  gt_file_delete(mil->pvt->matchfile);
  gt_str_delete(mil->pvt->matchfilename);
  gt_str_delete(mil->pvt->idxfilename);
  gt_str_delete(mil->pvt->queryfilename);
  gt_str_delete(mil->pvt->linebuf);
  gt_encseq_delete(mil->pvt->es1);
  gt_encseq_delete(mil->pvt->es2);
  gt_hashmap_delete(mil->pvt->desc_to_seqno);
  gt_free(mil->pvt);
}

const GtMatchIteratorClass* gt_match_iterator_last_class(void)
{
  static const GtMatchIteratorClass *mpc;
  gt_class_alloc_lock_enter();
  if (!mpc) {
    mpc = gt_match_iterator_class_new(sizeof (GtMatchIteratorLast),
                               gt_match_iterator_last_next,
                               gt_match_iterator_last_free);
  }
  gt_class_alloc_lock_leave();
  return mpc;
}
