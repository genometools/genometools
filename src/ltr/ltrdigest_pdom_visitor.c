/*
  Copyright (c) 2013-2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2013      Center for Bioinformatics, University of Hamburg
  Copyright (c)      2014 Genome Research Ltd.

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

#include <signal.h>
#include <string.h>
#ifndef S_SPLINT_S
#include <ctype.h>
#include <sys/types.h>
#include <unistd.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif
#endif
#include "core/array_api.h"
#include "core/codon_api.h"
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/cstr_api.h"
#include "core/cstr_array.h"
#include "core/grep_api.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/symbol_api.h"
#include "core/thread_api.h"
#include "core/translator_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/globalchaining.h"
#include "extended/reverse_api.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_pdom_visitor.h"

#define GT_HMMER_BUF_LEN  122

struct GtLTRdigestPdomVisitor {
  const GtNodeVisitor parent_instance;
  GtPdomModelSet *model;
  GtRegionMapping *rmap;
  double eval_cutoff;
  GtFeatureNode *ltr_retrotrans;
  GtStr *fwd[3], *rev[3];
  unsigned int chain_max_gap_length;
  GtUword leftLTR_5, rightLTR_3;
  GtPdomCutoff cutoff;
  GtStr *cmdline, *tag;
  bool output_all_chains;
  char **args;
  const char *root_type;
};

typedef struct {
  GtStrand strand;
  unsigned int frame;
  GtStr *cur_model;
  GtHashmap *models;
} GtHMMERParseStatus;

typedef struct {
  GtArray *fwd_hits,
          *rev_hits;
  GtUword last_array_size_fwd,
                last_array_size_rev;
  double best_rev,
         best_fwd;
  char *modelname;
} GtHMMERModelHit;

typedef struct {
  GtUword hmmfrom, hmmto, alifrom, alito, frame;
  double evalue, score;
  GtStrand strand;
  GtArray *chains;
  bool reported;
  GtStr *alignment, *aastring;
} GtHMMERSingleHit;

#ifndef _WIN32
static void gt_hmmer_model_hit_delete(GtHMMERModelHit *mh);
#endif

#ifndef _WIN32
static GtHMMERParseStatus* gt_hmmer_parse_status_new(void)
{
  GtHMMERParseStatus *s;
  s = gt_calloc((size_t) 1, sizeof (GtHMMERParseStatus));
  s->cur_model = gt_str_new();
  s->models = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                             (GtFree) gt_hmmer_model_hit_delete);
  return s;
}
#endif

#ifndef _WIN32
static void gt_hmmer_parse_status_add_hit(GtHMMERParseStatus *s,
                                          GtHMMERSingleHit *hit)
{
  GtHMMERModelHit *mh;
  gt_assert(s);
  if (!(mh = gt_hashmap_get(s->models, gt_str_get(s->cur_model)))) {
    mh = gt_calloc((size_t) 1, sizeof (*mh));
    mh->fwd_hits = gt_array_new(sizeof (GtHMMERSingleHit*));
    mh->rev_hits = gt_array_new(sizeof (GtHMMERSingleHit*));
    mh->best_rev = mh->best_fwd = DBL_MAX;
    mh->modelname = gt_cstr_dup(gt_str_get(s->cur_model));
    gt_hashmap_add(s->models, mh->modelname, mh);
  }
  gt_assert(mh && mh->fwd_hits &&mh->rev_hits);
  if (hit->strand == GT_STRAND_FORWARD) {
    if (gt_double_compare(mh->best_fwd, hit->evalue) > 0)
      mh->best_fwd = hit->evalue;
    gt_array_add(mh->fwd_hits, hit);
  } else {
    if (gt_double_compare(mh->best_rev, hit->evalue) > 0)
      mh->best_rev = hit->evalue;
    gt_array_add(mh->rev_hits, hit);
  }
}
#endif

#ifndef _WIN32
static void gt_hmmer_parse_status_mark_frame_finished(GtHMMERParseStatus *s)
{
  GtHMMERModelHit *mh;
  gt_assert(s && s->models);
  mh = gt_hashmap_get(s->models, gt_str_get(s->cur_model));
  if (mh != NULL) {
    mh->last_array_size_fwd = gt_array_size(mh->fwd_hits);
    mh->last_array_size_rev = gt_array_size(mh->rev_hits);
  }
}
#endif

#ifndef _WIN32
static GtHMMERSingleHit* gt_hmmer_parse_status_get_hit(GtHMMERParseStatus *s,
                                                       GtUword i)
{
  GtHMMERModelHit *mh;
  gt_assert(s && s->models);
  mh = gt_hashmap_get(s->models, gt_str_get(s->cur_model));
  if (mh != NULL) {
    if (s->strand == GT_STRAND_FORWARD) {
      i += mh->last_array_size_fwd;
      gt_assert(i < gt_array_size(mh->fwd_hits));
      return *(GtHMMERSingleHit**) gt_array_get(mh->fwd_hits, i);
    } else {
      i += mh->last_array_size_rev;
      gt_assert(i < gt_array_size(mh->rev_hits));
      return *(GtHMMERSingleHit**) gt_array_get(mh->rev_hits, i);
    }
  } else return NULL;
}
#endif

GT_UNUSED static int pdom_printvals(void *key, void *val, GT_UNUSED void *data,
                                    GT_UNUSED GtError *err) {
  GtHMMERModelHit *mh = (GtHMMERModelHit*) val;
  GtUword i;
  printf(">>  %s\n", (char*) key);
  for (i = 0; i < gt_array_size(mh->fwd_hits); i++) {
    GtHMMERSingleHit *h = *(GtHMMERSingleHit**) gt_array_get(mh->fwd_hits, i);
    printf(">>f    "GT_WU"-"GT_WU" ("GT_WU"%c, eval %e)\n", h->alifrom,
           h->alito, h->frame, GT_STRAND_CHARS[h->strand], h->evalue);
  }
  for (i = 0; i < gt_array_size(mh->rev_hits); i++) {
    GtHMMERSingleHit *h = *(GtHMMERSingleHit**) gt_array_get(mh->rev_hits, i);
    printf(">>r    "GT_WU"-"GT_WU" ("GT_WU"%c, eval %e)\n", h->alifrom,
           h->alito, h->frame, GT_STRAND_CHARS[h->strand], h->evalue);
  }
  printf("bestf: %e, bestr: %e\n", mh->best_fwd, mh->best_rev);
  return 0;
}

GT_UNUSED static void gt_hmmer_parse_status_show(GtHMMERParseStatus *s)
{
  gt_assert(s);
  (void) gt_hashmap_foreach(s->models, pdom_printvals, NULL, NULL);
}

#ifndef _WIN32
static void gt_hmmer_model_hit_delete(GtHMMERModelHit *mh)
{
  GtUword i;
  if (!mh) return;
  for (i = 0; i < gt_array_size(mh->fwd_hits); i++) {
    GtHMMERSingleHit *h = *(GtHMMERSingleHit**)  gt_array_get(mh->fwd_hits, i);
    gt_str_delete(h->alignment);
    gt_str_delete(h->aastring);
    gt_array_delete(h->chains);
    gt_free(h);
  }
  gt_array_delete(mh->fwd_hits);
  for (i = 0; i < gt_array_size(mh->rev_hits); i++) {
    GtHMMERSingleHit *h = *(GtHMMERSingleHit**)  gt_array_get(mh->rev_hits, i);
    gt_str_delete(h->alignment);
    gt_str_delete(h->aastring);
    gt_array_delete(h->chains);
    gt_free(h);
  }
  gt_array_delete(mh->rev_hits);
  gt_free(mh);
}
#endif

#ifndef _WIN32
static void gt_hmmer_parse_status_delete(GtHMMERParseStatus *s)
{
  if (!s) return;
  gt_str_delete(s->cur_model);
  gt_hashmap_delete(s->models);
  gt_free(s);
}
#endif

const GtNodeVisitorClass* gt_ltrdigest_pdom_visitor_class(void);

#define gt_ltrdigest_pdom_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrdigest_pdom_visitor_class(), GV)

static inline int pdom_parser_get_next_line(char *buf, FILE *instream,
                                            GtError *err) {
  int had_err = 0;
  char *bufp = buf;
  if (fgets(buf, GT_HMMER_BUF_LEN, instream) == NULL) {
    if (feof(instream)) {
      memset(buf, (int) '\0', (size_t) GT_HMMER_BUF_LEN);
      return 0;
    } else if (ferror(instream)) {
      gt_error_set(err, "error reading from input stream");
      return -1;
    }
  }
  while (!had_err && bufp && (buf[0] == '#' || buf[0] == '\n'))
    bufp = fgets(buf, GT_HMMER_BUF_LEN, instream);
  (void) gt_cstr_rtrim(buf, '\n');
  (void) gt_cstr_rtrim(buf, ' ');
  return 0;
}

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_statistics(GT_UNUSED
                                                      GtLTRdigestPdomVisitor
                                                                            *lv,
                                                     char *buf, FILE *instream,
                                                     GtError *err)
{
  int had_err = 0;
  gt_assert(lv && instream);
  gt_error_check(err);
  while (!had_err && (buf[0] != '/' && buf[1] != '/'))
    had_err = pdom_parser_get_next_line(buf, instream, err);
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_scores(GT_UNUSED
                                                  GtLTRdigestPdomVisitor *lv,
                                                  char *buf, FILE *instream,
                                                  GtError *err)
{
  int had_err = 0;
  gt_assert(lv && instream);
  gt_error_check(err);
  had_err = pdom_parser_get_next_line(buf, instream, err);
  gt_assert(had_err || buf != NULL);
  if (!had_err && strncmp("Scores", buf, (size_t) 6) != 0) {
    gt_error_set(err, "expected 'Scores:' at beginning of new scores "
                      "section, '%s' read instead", buf);
    had_err = -1;
  }
  while (!had_err && strncmp("Domain annotation", buf, (size_t) 17))
    had_err = pdom_parser_get_next_line(buf, instream, err);
  return had_err;
}
#endif

#define gt_ltrdigest_pdom_visitor_isgap(c) \
        ((c) == ' ' || (c) == '.' || (c) == '_' \
                    || (c) == '-' || (c) == '~')

GT_UNUSED static void gt_ltrdigest_pdom_visitor_add_aaseq(const char *str,
                                                          GtStr *dest)
{
  GtUword i;
  gt_assert(str && dest);
  for (i = 0; str[i] != '\0'; i++) {
    if (!gt_ltrdigest_pdom_visitor_isgap(str[i])) {
      /* replace stop codons by 'X'es */
      if (str[i] == '*') {
        gt_str_append_char(dest, 'X');
      } else {
        gt_str_append_char(dest, toupper(str[i]));
      }
    }
  }
}

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_alignments(GT_UNUSED
                                                      GtLTRdigestPdomVisitor
                                                                            *lv,
                                                     GtHMMERParseStatus *status,
                                                     char *buf,
                                                     FILE *instream,
                                                     GtError *err)
{
  int had_err = 0, cur_domain = GT_UNDEF_INT, line = GT_UNDEF_INT;
  int mod_val = 4;
  GtHMMERSingleHit *hit = NULL;
  gt_assert(lv && instream && status);
  gt_error_check(err);
  had_err = pdom_parser_get_next_line(buf, instream, err);
  gt_assert(buf != NULL);
  while (!had_err && strncmp("Internal pipeline statistics",
                             buf, (size_t) 28) &&
                     strncmp(">>", buf, (size_t) 2)) {
    if ((buf[2] == '=' && buf[3] == '=')) {
      buf[17] = '\0';
      cur_domain = atoi(buf+12);
      gt_assert(cur_domain != GT_UNDEF_INT && cur_domain > 0);
      hit = gt_hmmer_parse_status_get_hit(status,
                                          (GtUword) cur_domain - 1);
      gt_assert(hit && !hit->alignment);
      hit->alignment = gt_str_new();
      hit->aastring = gt_str_new();
      mod_val = 4;
      line = 0;
    } else {
      bool junk_match = false;
      (void) gt_grep(&junk_match, "(CS|RF)$", buf, NULL);
      if (!junk_match) {
        gt_assert(hit && hit->alignment);
        gt_str_append_cstr(hit->alignment, buf);
        gt_str_append_char(hit->alignment, '\n');
        switch (line % mod_val) {
          case 1:
            gt_str_append_char(hit->alignment, '\n');
            break;
          case 2:
            {
              GT_UNUSED char *b = buf;
              b = strtok(buf, " ");
              gt_assert(strspn(b, "012+-") == (size_t) 2);
              b = strtok(NULL, " ");
              gt_assert(strlen(b) > 0);
              b = strtok(NULL, " ");
              gt_ltrdigest_pdom_visitor_add_aaseq(b, hit->aastring);
            }
            break;
        }
        line++;
      }
    }
    had_err = pdom_parser_get_next_line(buf, instream, err);
  }
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_domainhits(GtLTRdigestPdomVisitor
                                                                            *lv,
                                                     GtHMMERParseStatus *status,
                                                     char *buf,
                                                     FILE *instream,
                                                     GtError *err)
{
  int had_err = 0;
  GtUword i, nof_targets = 0, nof_hits = 0;
  gt_assert(lv && instream && status);
  gt_error_check(err);

  had_err = pdom_parser_get_next_line(buf, instream, err);
  gt_assert(buf != NULL);
  while (!had_err && strncmp("Internal", buf, (size_t) 8)) {
    GtUword no, hmmfrom, hmmto, alifrom, alito;
    double score, evalue;
    char threshold_ok = '-';
    if ((buf[0] == '>' && buf[1] == '>')) {
      char *b = buf;
      b = strtok(buf+3, " ");
      gt_str_reset(status->cur_model);
      gt_str_append_cstr(status->cur_model, b);
      had_err = pdom_parser_get_next_line(buf, instream, err);
      if (!had_err && strncmp("   [No individual", buf, (size_t) 17)) {
        for (i = 0UL; i < 2UL && !had_err; i++)
          had_err = pdom_parser_get_next_line(buf, instream, err);
      }
      nof_targets++;
      nof_hits = 0UL;
      gt_hmmer_parse_status_mark_frame_finished(status);
    }
    while (!had_err &&
             8 == sscanf(buf, ""GT_WU" %c %lf %*f %*f %lf "GT_WU" "GT_WU" %*s "
                         GT_WU" "GT_WU"", &no,  &threshold_ok, &score, &evalue,
                         &hmmfrom, &hmmto, &alifrom, &alito)) {
      GtHMMERSingleHit *shit = gt_calloc((size_t) 1, sizeof (*shit));
      shit->hmmfrom = hmmfrom;
      shit->hmmto = hmmto;
      shit->alifrom = alifrom;
      shit->alito = alito;
      shit->score = score;
      shit->evalue = evalue;
      shit->strand = status->strand;
      shit->frame = (GtUword) status->frame;
      shit->reported = (threshold_ok == '!');
      shit->chains = gt_array_new(sizeof (GtUword));
      gt_hmmer_parse_status_add_hit(status, shit);
      nof_hits++;
      had_err = pdom_parser_get_next_line(buf, instream, err);
    }
    if (!had_err) {
      if (nof_hits > 0)
        had_err = gt_ltrdigest_pdom_visitor_parse_alignments(lv, status, buf,
                                                             instream, err);
      else
        had_err = pdom_parser_get_next_line(buf, instream, err);
    }
  }
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_query(GtLTRdigestPdomVisitor *lv,
                                                 GtHMMERParseStatus *status,
                                                 bool *end,
                                                 FILE *instream, GtError *err)
{
  int had_err = 0;
  char buf[GT_HMMER_BUF_LEN];
  gt_assert(lv && instream && status);
  gt_error_check(err);

  had_err = pdom_parser_get_next_line(buf, instream, err);
  if (!had_err && strncmp("Query:", buf, (size_t) 6) != 0) {
    *end = true;
  }
  if (!had_err && !(*end)) {
    status->strand = gt_strand_get(buf[14]);
    buf[14] = '\0';
    status->frame = (unsigned) atoi(buf+13);
  }
  if (!had_err && !(*end)) {
    had_err = gt_ltrdigest_pdom_visitor_parse_scores(lv, buf, instream, err);
  }
  if (!had_err && !(*end)) {
    had_err = gt_ltrdigest_pdom_visitor_parse_domainhits(lv, status, buf,
                                                         instream, err);
  }
  if (!had_err && !(*end)) {
    had_err = gt_ltrdigest_pdom_visitor_parse_statistics(lv, buf, instream,
                                                         err);
  }
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_parse_output(GtLTRdigestPdomVisitor *lv,
                                                  GtHMMERParseStatus *status,
                                                  FILE *instream, GtError *err)
{
  int had_err = 0;
  bool end = false;
  gt_assert(lv && instream && status);
  gt_error_check(err);
  while (!had_err && !end) {
    had_err = gt_ltrdigest_pdom_visitor_parse_query(lv, status, &end,
                                                    instream, err);
  }
  /* gt_hmmer_parse_status_show(status); */
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_fragcmp(const void *frag1,
                                             const void *frag2)
{
  GtFragment *f1 = (GtFragment*) frag1;
  GtFragment *f2 = (GtFragment*) frag2;
  if (f1->startpos2 == f2->startpos2)
    return 0;
  else return (f1->startpos2 < f2->startpos2 ? -1 : 1);
}
#endif

#ifndef _WIN32
static void gt_ltrdigest_pdom_visitor_chainproc(GtChain *c, GtFragment *f,
                                             GT_UNUSED GtUword nof_frags,
                                             GT_UNUSED GtUword gap_length,
                                             void *data)
{
  GtUword i,
                *chainno = (GtUword*) data;
  gt_log_log("resulting chain has "GT_WD" GtFragments, score "GT_WD"",
             gt_chain_size(c),
             gt_chain_get_score(c));
  for (i = 0; i < gt_chain_size(c); i++) {
    GtFragment frag;
    frag = f[gt_chain_get_fragnum(c, i)];
    gt_log_log("("GT_WU" "GT_WU") ("GT_WU" "GT_WU")", frag.startpos1,
               frag.endpos1, frag.startpos2, frag.endpos2);
    gt_array_add(((GtHMMERSingleHit*) frag.data)->chains, *chainno);
  }
  (*chainno)++;
  gt_log_log("\n");
}
#endif

#ifndef _WIN32
static GtRange gt_ltrdigest_pdom_visitor_coords(GtLTRdigestPdomVisitor *lv,
                                              const GtHMMERSingleHit *singlehit)
{
  GtRange retrng;
  gt_assert(singlehit);
  switch (singlehit->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      retrng.start = lv->leftLTR_5 + (singlehit->alifrom - 1) * GT_CODON_LENGTH
                        + (GtUword) singlehit->frame;
      retrng.end   =  retrng.start + (singlehit->alito - singlehit->alifrom + 1)
                           * GT_CODON_LENGTH;
      break;
    case GT_STRAND_REVERSE:
      retrng.start = lv->rightLTR_3 - (singlehit->alito) * GT_CODON_LENGTH
                        - (GtUword) singlehit->frame;
      retrng.end   =  retrng.start + (singlehit->alito - singlehit->alifrom + 1)
                           * GT_CODON_LENGTH;
      break;
  }
  retrng.start++; retrng.end++;  /* GFF3 is 1-based */
  return retrng;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_attach_hit(GtLTRdigestPdomVisitor *lv,
                                                GtHMMERModelHit *modelhit,
                                                GtHMMERSingleHit *singlehit)
{
  GT_UNUSED GtUword i;
  GtGenomeNode *gf;
  int had_err = 0;
  GtRange rrng;
  gt_assert(lv && singlehit);

  rrng = gt_ltrdigest_pdom_visitor_coords(lv, singlehit);

  if (gt_array_size(singlehit->chains) > 0 || lv->output_all_chains) {
    char buf[32];
    gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                      lv->ltr_retrotrans),
                             gt_ft_protein_match,
                             rrng.start,
                             rrng.end,
                             singlehit->strand);
    gt_genome_node_add_user_data((GtGenomeNode*) gf, "pdom_alignment",
                                 gt_str_ref(singlehit->alignment),
                                 (GtFree) gt_str_delete);
    gt_genome_node_add_user_data((GtGenomeNode*) gf, "pdom_aaseq",
                                 gt_str_ref(singlehit->aastring),
                                 (GtFree) gt_str_delete);
    gt_feature_node_set_source((GtFeatureNode*) gf, lv->tag);
    gt_feature_node_set_score((GtFeatureNode*) gf, (float) singlehit->evalue);
    (void) snprintf(buf, (size_t) 32, "%d", (int) singlehit->frame);
    gt_feature_node_add_attribute((GtFeatureNode*) gf,
                                    "reading_frame", buf);
    if (modelhit->modelname != NULL) {
      gt_feature_node_add_attribute((GtFeatureNode*) gf, "name",
                                    modelhit->modelname);
    }
    if (gt_array_size(singlehit->chains) > 1UL && lv->output_all_chains) {
      GtStr *buffer;
      GtUword j;
      gt_assert(singlehit->chains != NULL);
      buffer = gt_str_new();
      for (j = 0UL; j < gt_array_size(singlehit->chains); j++) {
        gt_str_append_cstr(buffer, modelhit->modelname);
        gt_str_append_char(buffer, ':');
        gt_str_append_ulong(buffer,
                          *(GtUword*) gt_array_get(singlehit->chains, j));
        if (j != gt_array_size(singlehit->chains) - 1) {
          gt_str_append_char(buffer, ',');
        }
      }
      gt_feature_node_set_attribute((GtFeatureNode*) gf, "chains",
                                    gt_str_get(buffer));
      gt_str_delete(buffer);
    }
    gt_feature_node_add_child(lv->ltr_retrotrans, (GtFeatureNode*) gf);
  }
  gt_array_delete(singlehit->chains);
  singlehit->chains = NULL;
  return had_err;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_process_hit(GT_UNUSED void *key, void *val,
                                                 void *data,
                                                 GT_UNUSED GtError *err)
{
  GtHMMERModelHit *mh = (GtHMMERModelHit*) val;
  GtLTRdigestPdomVisitor *lv = (GtLTRdigestPdomVisitor*) data;
  const char *mdl = (const char*) key;
  GtArray *hits = NULL;
  GtUword nof_hits;
  GtFragment *frags;

  if (gt_double_compare(mh->best_fwd, mh->best_rev) <= 0)
    hits = mh->fwd_hits;
  else
    hits = mh->rev_hits;
  gt_assert(hits);
  nof_hits = gt_array_size(hits);
  if (nof_hits == 0) return 0;

  if (nof_hits > 1UL) {
    GtUword i, chainno;
    frags = gt_malloc((size_t) nof_hits * sizeof (GtFragment));
    for (i = 0; i < nof_hits; i++) {
      GtHMMERSingleHit *h = *(GtHMMERSingleHit**) gt_array_get(hits, i);
      gt_assert(h);
      frags[i].startpos1 = h->hmmfrom;
      frags[i].endpos1   = h->hmmto;
      frags[i].startpos2 = h->alifrom;
      frags[i].endpos2   = h->alito;
      frags[i].weight    = (GtWord) (h->alito - h->alifrom + 1) * h->score;
      frags[i].data      = h;
    }
    qsort(frags, (size_t) nof_hits, sizeof (GtFragment),
          gt_ltrdigest_pdom_visitor_fragcmp);
    gt_log_log("%s: chaining "GT_WU" frags", mdl, nof_hits);
    gt_globalchaining_max(frags, nof_hits,
                         (GtUword) lv->chain_max_gap_length,
                         gt_ltrdigest_pdom_visitor_chainproc, &chainno);
    gt_free(frags);
    for (i = 0; i < nof_hits; i++) {
      GtHMMERSingleHit *h = *(GtHMMERSingleHit**) gt_array_get(hits, i);
      (void) gt_ltrdigest_pdom_visitor_attach_hit(lv, mh, h);
    }
  } else {
    GtUword chainno = 0UL;
    GtHMMERSingleHit *h = *(GtHMMERSingleHit**) gt_array_get(hits, 0);
    gt_array_add(h->chains, chainno);
    (void) gt_ltrdigest_pdom_visitor_attach_hit(lv, mh, h);
  }

  return 0;
}
#endif

#ifndef _WIN32
static int gt_ltrdigest_pdom_visitor_process_hits(GtLTRdigestPdomVisitor *lv,
                                                  GtHMMERParseStatus *status,
                                                  GtError *err)
{
  int had_err = 0;
  gt_assert(lv && status && status->models);
  gt_error_check(err);

  had_err = gt_hashmap_foreach(status->models,
                               gt_ltrdigest_pdom_visitor_process_hit,
                               lv, err);

  return had_err;
}
#endif

static int gt_ltrdigest_pdom_visitor_choose_strand(GtLTRdigestPdomVisitor *lv)
{
  int had_err = 0;
  double log_eval_fwd = 0.0,
         log_eval_rev = 0.0;
  GtFeatureNodeIterator *fni;
  GtStrand strand;
  double score;
  bool seen_fwd = false,
       seen_rev = false;
  GtFeatureNode *curnode = NULL;
  GtUword i;
  GtArray *to_delete;

  fni = gt_feature_node_iterator_new(lv->ltr_retrotrans);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_protein_match) == 0) {
      strand = gt_feature_node_get_strand(curnode);
      score = (double) gt_feature_node_get_score(curnode);
      if (strand == GT_STRAND_FORWARD) {
        log_eval_fwd += log(score);
        seen_fwd = true;
      } else if (strand == GT_STRAND_REVERSE) {
        log_eval_rev += log(score);
        seen_rev = true;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (seen_rev && !seen_fwd)
    gt_feature_node_set_strand(lv->ltr_retrotrans, GT_STRAND_REVERSE);
  else if (!seen_rev && seen_fwd)
    gt_feature_node_set_strand(lv->ltr_retrotrans, GT_STRAND_FORWARD);
  else if (!seen_rev && !seen_fwd)
    return had_err;
  else {
    gt_assert(seen_rev && seen_fwd);
    if (gt_double_compare(log_eval_fwd, log_eval_rev) < 0)
      strand = GT_STRAND_FORWARD;
    else
      strand = GT_STRAND_REVERSE;
    gt_feature_node_set_strand(lv->ltr_retrotrans, strand);

    to_delete = gt_array_new(sizeof (GtFeatureNode*));
    fni = gt_feature_node_iterator_new(lv->ltr_retrotrans);
    while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
      if (strcmp(gt_feature_node_get_type(curnode),
                 gt_ft_protein_match) == 0) {
        if (strand != gt_feature_node_get_strand(curnode)) {
          gt_array_add(to_delete, curnode);
        }
      }
    }
    gt_feature_node_iterator_delete(fni);
    gt_assert(gt_array_size(to_delete) > 0);
    for (i = 0; i < gt_array_size(to_delete); i++) {
      gt_feature_node_remove_leaf(lv->ltr_retrotrans,
                                  *(GtFeatureNode**) gt_array_get(to_delete,
                                                                  i));
    }
    gt_array_delete(to_delete);
  }
  return had_err;
}

static int gt_ltrdigest_pdom_visitor_feature_node(GtNodeVisitor *nv,
                                                  GtFeatureNode *fn,
                                                  GtError *err)
{
  GtLTRdigestPdomVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL;
  int had_err = 0;
  GtRange rng;
  GtUword i;
  lv = gt_ltrdigest_pdom_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  /* traverse annotation subgraph and find LTR element */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (strcmp(gt_feature_node_get_type(curnode), lv->root_type) == 0) {
      lv->ltr_retrotrans = curnode;
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (!had_err && lv->ltr_retrotrans != NULL) {
    GtCodonIterator *ci;
    GtTranslator *tr;
    GtTranslatorStatus status;
    GtUword seqlen;
    char translated, *rev_seq;
#ifndef _WIN32
    FILE *instream;
    GtHMMERParseStatus *pstatus;
#endif
    unsigned int frame;
    GtStr *seq;

    seq = gt_str_new();
    rng = gt_genome_node_get_range((GtGenomeNode*) lv->ltr_retrotrans);
    lv->leftLTR_5 = rng.start - 1;
    lv->rightLTR_3 = rng.end - 1;
    seqlen = gt_range_length(&rng);

    had_err = gt_extract_feature_sequence(seq,
                                          (GtGenomeNode*) lv->ltr_retrotrans,
                                          lv->root_type,
                                          false, NULL, NULL, lv->rmap, err);

    if (!had_err) {
      for (i = 0UL; i < 3UL; i++) {
        gt_str_reset(lv->fwd[i]);
        gt_str_reset(lv->rev[i]);
      }

      /* create translations */
      ci = gt_codon_iterator_simple_new(gt_str_get(seq), seqlen, NULL);
      gt_assert(ci);
      tr = gt_translator_new(ci);
      status = gt_translator_next(tr, &translated, &frame, err);
      while (status == GT_TRANSLATOR_OK && translated) {
        gt_str_append_char(lv->fwd[frame], translated);
        status = gt_translator_next(tr, &translated, &frame, NULL);
      }
      if (status == GT_TRANSLATOR_ERROR) had_err = -1;
      if (!had_err) {
        rev_seq = gt_malloc((size_t) seqlen * sizeof (char));
        strncpy(rev_seq, gt_str_get(seq), (size_t) seqlen * sizeof (char));
        (void) gt_reverse_complement(rev_seq, seqlen, NULL);
        gt_codon_iterator_delete(ci);
        ci = gt_codon_iterator_simple_new(rev_seq, seqlen, NULL);
        gt_translator_set_codon_iterator(tr, ci);
        status = gt_translator_next(tr, &translated, &frame, err);
        while (status == GT_TRANSLATOR_OK && translated) {
          gt_str_append_char(lv->rev[frame], translated);
          status = gt_translator_next(tr, &translated, &frame, NULL);
        }
        if (status == GT_TRANSLATOR_ERROR) had_err = -1;
        gt_free(rev_seq);
      }
      gt_codon_iterator_delete(ci);
      gt_translator_delete(tr);
    }

    /* run HMMER and handle results */
    if (!had_err) {
#ifndef _WIN32
      int pid, pc[2], cp[2];
      GT_UNUSED int rval;

      (void) signal(SIGCHLD, SIG_IGN); /* XXX: for now, ignore child's
                                               exit status */
      rval = pipe(pc);
      gt_assert(rval == 0);
      rval = pipe(cp);
      gt_assert(rval == 0);

      switch ((pid = (int) fork())) {
        case -1:
          perror("Can't fork");
          exit(1);   /* XXX: error handling */
        case 0:    /* child */
          (void) close(1);    /* close current stdout. */
          rval = dup(cp[1]);  /* make stdout go to write end of pipe. */
          (void) close(0);    /* close current stdin. */
          rval = dup(pc[0]);  /* make stdin come from read end of pipe. */
          (void) close(pc[0]);
          (void) close(pc[1]);
          (void) close(cp[0]);
          (void) close(cp[1]);
          (void) execvp("hmmscan", lv->args); /* XXX: read path from env */
          perror("couldn't execute hmmscan!");
          exit(1);
        default:    /* parent */
          for (i = 0UL; i < 3UL; i++) {
            char buf[5];
            GT_UNUSED ssize_t written;
            (void) sprintf(buf, ">"GT_WU"%c\n", i, '+');
            written = write(pc[1], buf, 4 * sizeof (char));
            written = write(pc[1], gt_str_get(lv->fwd[i]),
                            (size_t) gt_str_length(lv->fwd[i]) * sizeof (char));
            written = write(pc[1], "\n", 1 * sizeof (char));
            (void) sprintf(buf, ">"GT_WU"%c\n", i, '-');
            written = write(pc[1], buf, 4 * sizeof (char));
            written = write(pc[1], gt_str_get(lv->rev[i]),
                            (size_t) gt_str_length(lv->rev[i]) * sizeof (char));
            written = write(pc[1], "\n", 1 * sizeof (char));
          }
          (void) close(pc[0]);
          (void) close(pc[1]);
          (void) close(cp[1]);
          instream = fdopen(cp[0], "r");
          pstatus = gt_hmmer_parse_status_new();
          had_err = gt_ltrdigest_pdom_visitor_parse_output(lv, pstatus,
                                                           instream, err);
          (void) fclose(instream);
          if (!had_err)
            had_err = gt_ltrdigest_pdom_visitor_process_hits(lv, pstatus, err);
          gt_hmmer_parse_status_delete(pstatus);
      }
#else
      /* XXX */
      gt_error_set(err, "HMMER call not implemented on Windows\n");
      had_err = -1;
#endif
    }
    gt_str_delete(seq);
  }
  if (!had_err)
    had_err = gt_ltrdigest_pdom_visitor_choose_strand(lv);
  return had_err;
}

void gt_ltrdigest_pdom_visitor_free(GtNodeVisitor *nv)
{
  GtLTRdigestPdomVisitor *lv;
  GtUword i;
  if (!nv) return;
  lv = gt_ltrdigest_pdom_visitor_cast(nv);
  for (i = 0UL; i < 3UL; i++) {
    gt_str_delete(lv->fwd[i]);
    gt_str_delete(lv->rev[i]);
  }
  gt_str_delete(lv->cmdline);
  gt_str_delete(lv->tag);
  gt_cstr_array_delete(lv->args);
}

const GtNodeVisitorClass* gt_ltrdigest_pdom_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRdigestPdomVisitor),
                                   gt_ltrdigest_pdom_visitor_free,
                                   NULL,
                                   gt_ltrdigest_pdom_visitor_feature_node,
                                   NULL,
                                   NULL,
                                   NULL);
  }
  return nvc;
}

void gt_ltrdigest_pdom_visitor_output_all_chains(GtLTRdigestPdomVisitor *lv)
{
  gt_assert(lv);
  lv->output_all_chains = true;
}

void gt_ltrdigest_pdom_visitor_set_root_type(GtLTRdigestPdomVisitor *lv,
                                             const char *type)
{
  gt_assert(lv && type);
  lv->root_type = gt_symbol(type);
}

void gt_ltrdigest_pdom_visitor_set_source_tag(GtLTRdigestPdomVisitor *lv,
                                              const char *tag)
{
  gt_assert(lv && tag && lv->tag);
  gt_str_reset(lv->tag);
  gt_str_append_cstr(lv->tag, tag);
}

GtNodeVisitor* gt_ltrdigest_pdom_visitor_new(GtPdomModelSet *model,
                                             double eval_cutoff,
                                             unsigned int chain_max_gap_length,
                                             GtPdomCutoff cutoff,
                                             GtRegionMapping *rmap,
                                             GtError *err)
{
  GtNodeVisitor *nv;
  GtLTRdigestPdomVisitor *lv;
  GtStr *cmd;
  int had_err = 0, i, rval;
  gt_assert(model && rmap);

  rval = system("hmmscan -h > /dev/null");
  if (rval == -1) {
    gt_error_set(err, "error executing system(hmmscan)");
    return NULL;
  }
#ifndef _WIN32
  if (WEXITSTATUS(rval) != 0) {
    gt_error_set(err, "cannot find the hmmscan executable in PATH");
    return NULL;
  }
#else
  /* XXX */
  gt_error_set(err, "hmmscan for Windows not implemented");
  return NULL;
#endif

  nv = gt_node_visitor_create(gt_ltrdigest_pdom_visitor_class());
  lv = gt_ltrdigest_pdom_visitor_cast(nv);
  lv->eval_cutoff = eval_cutoff;
  lv->cutoff = cutoff;
  lv->chain_max_gap_length = chain_max_gap_length;
  lv->rmap = rmap;
  lv->output_all_chains = false;
  lv->tag = gt_str_new_cstr("GenomeTools");
  lv->root_type = gt_symbol(gt_ft_LTR_retrotransposon);

  for (i = 0; i < 3; i++) {
    lv->fwd[i] = gt_str_new();
    lv->rev[i] = gt_str_new();
  }

  if (!had_err) {
    cmd = gt_str_new_cstr("hmmscan --cpu ");
    gt_str_append_uint(cmd, gt_jobs);
    gt_str_append_cstr(cmd, " ");
    switch (cutoff) {
      case GT_PHMM_CUTOFF_GA:
        gt_str_append_cstr(cmd, "--cut_ga");
        break;
      case GT_PHMM_CUTOFF_TC:
        gt_str_append_cstr(cmd, "--cut_tc");
        break;
      case GT_PHMM_CUTOFF_NONE:
        gt_str_append_cstr(cmd, "--domE ");
        gt_str_append_double(cmd, eval_cutoff, 50);
        break;
    }
    gt_str_append_cstr(cmd, " ");
    gt_str_append_cstr(cmd, gt_pdom_model_set_get_filename(model));
    gt_str_append_cstr(cmd, " -");
    lv->cmdline = cmd;
    lv->args = gt_cstr_split(gt_str_get(lv->cmdline), ' ');
    gt_log_log("HMMER cmdline: %s", gt_str_get(cmd));
  }
  return nv;
}
