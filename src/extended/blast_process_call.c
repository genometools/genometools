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

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "core/ma.h"
#include "core/str_api.h"
#include "extended/blast_process_call.h"

struct GtBlastProcessCall
{
  GtStr *str;
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
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_all_nucl(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastall -p blastn");
  call->nucl = call->all = true;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_nucl(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastn");
  call->nucl = true;
  call->all = false;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_all_prot(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastall -p blastp");
  call->nucl = false;
  call->all = true;
  return call;
}

GtBlastProcessCall *gt_blast_process_call_new_prot(void)
{
  GtBlastProcessCall *call =
    gt_blast_process_call_new("blastall -p blastp");
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
  gt_assert(call->nucl);
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

void gt_blast_process_call_set_opt(GtBlastProcessCall *call,
                                   const char *opt)
{
  gt_str_append_cstr(call->str, opt);
}

FILE *gt_blast_process_call_run(GtBlastProcessCall *call, GtError *err)
{
  FILE *blastout;
  gt_assert(call->query && call->db);
  blastout = popen(gt_str_get(call->str), "r");
  if (blastout == NULL) {
    gt_error_set(err, "Could not run BLAST process: %s", strerror(errno));
    return NULL;
  }
  return blastout;
}

void gt_blast_process_call_delete(GtBlastProcessCall *call)
{
  if (call != NULL) {
    gt_str_delete(call->str);
    gt_free(call);
  }
}
