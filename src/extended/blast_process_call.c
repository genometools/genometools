/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "core/ma.h"
#include "core/str_api.h"
#include "extended/blast_process_call.h"
#include "core/log_api.h"

struct GtBlastProcessCall
{
  GtStr *str;
  char *version_call;
  bool all,
       db,
       evalue,
       gapextend,
       gapopen,
       nucl,
       num_threads,
       outfmt,
       penalty,
       query,
       reward,
       word_size,
       xdrop_gap_final;
};

static
GtBlastProcessCall *gt_blast_process_call_new(const char *which)
{
  GtBlastProcessCall *call = gt_malloc(sizeof (*call));
  char *env;

  env = getenv("GT_BLAST_PATH");
  if (env != NULL) {
    call->str = gt_str_new_cstr(env);
    gt_str_append_cstr(call->str, "/");
    gt_str_append_cstr(call->str, which);
  }
  else
    call->str = gt_str_new_cstr(which);
  call->all = false;
  call->db = false;
  call->evalue = false;
  call->gapextend = false;
  call->gapopen = false;
  call->nucl = false;
  call->num_threads = false;
  call->outfmt = false;
  call->penalty = false;
  call->query = false;
  call->reward = false;
  call->word_size = false;
  call->xdrop_gap_final = false;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_all_nucl(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastall -p blastn");
  call->nucl = call->all = true;
  call->version_call = "blastall -";
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_nucl(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastn");
  call->version_call = "blastn -version";
  call->nucl = true;
  call->all = false;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_all_prot(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastall -p blastp");
  call->version_call = "blastall -";
  call->nucl = false;
  call->all = true;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_prot(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastp");
  call->version_call = "blastp -version";
  call->nucl = call->all = false;
  return call;
}

void gt_blast_process_call_set_query(GtBlastProcessCall *call,
                                     const char *query)
{
  gt_assert(!call->query);
  call->query = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -i ");
  else
    gt_str_append_cstr(call->str, " -query ");
  gt_str_append_cstr(call->str, query);
}

void gt_blast_process_call_set_db(GtBlastProcessCall *call, const char *db)
{
  gt_assert(!call->db);
  call->db = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -d ");
  else
    gt_str_append_cstr(call->str, " -db ");
  gt_str_append_cstr(call->str, db);
}

void gt_blast_process_call_set_evalue(GtBlastProcessCall *call, double evalue)
{
  gt_assert(!call->evalue);
  call->evalue = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -e ");
  else
    gt_str_append_cstr(call->str, " -evalue ");
  gt_str_append_sci_double(call->str, evalue, 6);
}

void gt_blast_process_call_set_wordsize(GtBlastProcessCall *call, int word_size)
{
  gt_assert(!call->word_size);
  call->word_size = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -W ");
  else
    gt_str_append_cstr(call->str, " -word_size ");
  gt_str_append_int(call->str, word_size);
}

void gt_blast_process_call_set_gapopen(GtBlastProcessCall *call, int gapopen)
{
  gt_assert(!call->gapopen);
  call->gapopen = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -G ");
  else
    gt_str_append_cstr(call->str, " -gapopen ");
  gt_str_append_int(call->str, gapopen);
}

void gt_blast_process_call_set_gapextend(GtBlastProcessCall *call,
                                         int gapextend)
{
  gt_assert(!call->gapextend);
  call->gapextend = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -E ");
  else
    gt_str_append_cstr(call->str, " -gapextend ");
  gt_str_append_int(call->str, gapextend);
}

void gt_blast_process_call_set_penalty(GtBlastProcessCall *call, int penalty)
{
  gt_assert(!call->penalty);
  call->penalty = true;
  gt_assert(call->nucl);
  if (call->all)
    gt_str_append_cstr(call->str, " -q ");
  else
    gt_str_append_cstr(call->str, " -penalty ");
  gt_str_append_int(call->str, penalty);
}

void gt_blast_process_call_set_reward(GtBlastProcessCall *call, int reward)
{
  gt_assert(!call->reward);
  call->reward = true;
  gt_assert(call->nucl);
  if (call->all)
    gt_str_append_cstr(call->str, " -r ");
  else
    gt_str_append_cstr(call->str, " -reward ");
  gt_str_append_int(call->str, reward);
}

void gt_blast_process_call_set_num_threads(GtBlastProcessCall *call,
                                           int num_threads)
{
  gt_assert(!call->num_threads);
  call->num_threads = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -a ");
  else
    gt_str_append_cstr(call->str, " -num_threads ");
  gt_str_append_int(call->str, num_threads);
}

void gt_blast_process_call_set_xdrop_gap_final(GtBlastProcessCall *call,
                                               double xdrop_gap_final)
{
  gt_assert(!call->xdrop_gap_final);
  call->xdrop_gap_final = true;
  gt_assert(call->nucl);
  if (call->all)
    gt_str_append_cstr(call->str, " -Z ");
  else
    gt_str_append_cstr(call->str, " -xdrop_gap_final ");
  gt_str_append_double(call->str, xdrop_gap_final, 2);
}

void gt_blast_process_call_set_outfmt(GtBlastProcessCall *call,
                                      int outfmt)
{
  gt_assert(!call->outfmt);
  call->outfmt = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -m ");
  else
    gt_str_append_cstr(call->str, " -outfmt ");
  gt_str_append_int(call->str, outfmt);
}

void gt_blast_process_call_set_outfmt_tabular(GtBlastProcessCall *call)
{
  gt_assert(!call->outfmt);
  call->outfmt = true;
  if (call->all)
    gt_str_append_cstr(call->str, " -m 8");
  else
    gt_str_append_cstr(call->str, " -outfmt 6");
}

void gt_blast_process_call_set_opt(GtBlastProcessCall *call,
                                   const char *opt)
{
  gt_str_append_cstr(call->str, " ");
  gt_str_append_cstr(call->str, opt);
}

FILE *gt_blast_process_call_run(GtBlastProcessCall *call, GtError *err)
{
  int had_err = 0;
#ifndef _WIN32
  int errcode;
#endif
  FILE *blastout = NULL,
       *installcheck = NULL;
  gt_assert(call->query && call->db);

  installcheck = popen(call->version_call, "r");
  if (installcheck == NULL) {
    gt_error_set(err, "Could not open pipe to run %s: %s",
                 call->version_call, strerror(errno));
    had_err = -1;
  }
  /* this assures that we get the output if debugging is set, and also we
     prevent BROKEN_PIPE error if pclose(3) is called before the version call
     exits */
  if (!had_err) {
    char *newline = NULL;
    char line[BUFSIZ + 1];
    line[BUFSIZ] = '\0';
    while (fgets(line, (int) BUFSIZ, installcheck) != NULL) {
      if ((newline = strrchr(line, '\n')) != NULL) {
        *newline = '\0';
        newline = NULL;
      }
      gt_log_log("%.*s", BUFSIZ, line);
    }
  }
  if (!had_err) {
#ifndef _WIN32
    errcode = pclose(installcheck);
    if ((call->all && WEXITSTATUS(errcode) != 1) ||
        errcode != 0) {
      if (errno == ECHILD)
        gt_error_set(err, "Error calling %s.", call->version_call);
      else if (WEXITSTATUS(errcode) == 127)
        gt_error_set(err, "shell returned 127, BLAST not installed?");
      else
        gt_error_set(err, "%s error, returned %d",
                     call->version_call, WEXITSTATUS(errcode));
      had_err = -1;
    }
#endif
  }
  if (!had_err) {
    blastout = popen(gt_str_get(call->str), "r");
    if (blastout == NULL) {
      gt_error_set(err, "Could not open pipe to run BLAST process: %s",
                   strerror(errno));
    }
  }
  return blastout;
}

const char *gt_blast_process_call_get_call(GtBlastProcessCall *call)
{
  return gt_str_get(call->str);
}

void gt_blast_process_call_delete(GtBlastProcessCall *call)
{
  if (call != NULL) {
    gt_str_delete(call->str);
    gt_free(call);
  }
}
