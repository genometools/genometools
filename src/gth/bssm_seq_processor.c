/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#include "core/assert_api.h"
#include "core/bittab_api.h"
#include "core/compat.h"
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/file_api.h"
#include "core/fileutils_api.h"
#include "core/output_file.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "core/xposix.h"
#include "extended/reverse_api.h"
#include "gth/bssm_helper.h"
#include "gth/bssm_seq_processor.h"

#define GT_DIR  "GT_donor"
#define GC_DIR  "GC_donor"
#define AG_DIR  "AG_acceptor"

struct GthBSSMSeqProcessor {
  char *outdir;
  bool force,
       gcdonor;
  GtArray *exon_0,          /* squash */
          *exon_1,          /* squash */
          *exon_2,          /* squash */
          *intron_0,        /* squash */
          *intron_1,        /* squash */
          *intron_2,        /* squash */
          *intron_all,      /* squash */
          *i0_true_don_gt,  /* true donor and acceptor sites */
          *i0_true_don_gc,  /* true donor and acceptor sites */
          *i0_true_acc,     /* true donor and acceptor sites */
          *i1_true_don_gt,  /* true donor and acceptor sites */
          *i1_true_don_gc,  /* true donor and acceptor sites */
          *i1_true_acc,     /* true donor and acceptor sites */
          *i2_true_don_gt,  /* true donor and acceptor sites */
          *i2_true_don_gc,  /* true donor and acceptor sites */
          *i2_true_acc,     /* true donor and acceptor sites */
          *e0_false_don_gt, /* false donor and acceptor sites */
          *e0_false_don_gc, /* false donor and acceptor sites */
          *e0_false_acc,    /* false donor and acceptor sites */
          *e1_false_don_gt, /* false donor and acceptor sites */
          *e1_false_don_gc, /* false donor and acceptor sites */
          *e1_false_acc,    /* false donor and acceptor sites */
          *e2_false_don_gt, /* false donor and acceptor sites */
          *e2_false_don_gc, /* false donor and acceptor sites */
          *e2_false_acc,    /* false donor and acceptor sites */
          *i_false_don_gt,  /* false donor and acceptor sites */
          *i_false_don_gc,  /* false donor and acceptor sites */
          *i_false_acc;     /* false donor and acceptor sites */
  GtFile *gt_t1_fp,
         *gt_t2_fp,
         *gt_t0_fp,
         *gt_f1_fp,
         *gt_f2_fp,
         *gt_f0_fp,
         *gt_fi_fp,
         *gc_t1_fp,
         *gc_t2_fp,
         *gc_t0_fp,
         *gc_f1_fp,
         *gc_f2_fp,
         *gc_f0_fp,
         *gc_fi_fp,
         *ag_t1_fp,
         *ag_t2_fp,
         *ag_t0_fp,
         *ag_f1_fp,
         *ag_f2_fp,
         *ag_f0_fp,
         *ag_fi_fp;
};

typedef struct {
  GtFile *exon_0_fp,       /* squash */
         *exon_1_fp,       /* squash */
         *exon_2_fp,       /* squash */
         *intron_0_fp,     /* squash */
         *intron_1_fp,     /* squash */
         *intron_2_fp,     /* squash */
         *intron_all_fp,   /* squash */
         *i0_true_don_fp,  /* true donor and acceptor sites */
         *i0_true_acc_fp,  /* true donor and acceptor sites */
         *i1_true_don_fp,  /* true donor and acceptor sites */
         *i1_true_acc_fp,  /* true donor and acceptor sites */
         *i2_true_don_fp,  /* true donor and acceptor sites */
         *i2_true_acc_fp,  /* true donor and acceptor sites */
         *e0_false_don_fp, /* false donor and acceptor sites */
         *e0_false_acc_fp, /* false donor and acceptor sites */
         *e1_false_don_fp, /* false donor and acceptor sites */
         *e1_false_acc_fp, /* false donor and acceptor sites */
         *e2_false_don_fp, /* false donor and acceptor sites */
         *e2_false_acc_fp, /* false donor and acceptor sites */
         *i_false_don_fp,  /* false donor and acceptor sites */
         *i_false_acc_fp;  /* false donor and acceptor sites */
} IntermediateFiles;

void intermediate_files_delete(IntermediateFiles *ifp)
{
  if (!ifp) return;

  gt_file_delete(ifp->i_false_acc_fp);
  gt_file_delete(ifp->i_false_don_fp);
  gt_file_delete(ifp->e2_false_acc_fp);
  gt_file_delete(ifp->e2_false_don_fp);
  gt_file_delete(ifp->e1_false_acc_fp);
  gt_file_delete(ifp->e1_false_don_fp);
  gt_file_delete(ifp->e0_false_acc_fp);
  gt_file_delete(ifp->e0_false_don_fp);

  gt_file_delete(ifp->i2_true_acc_fp);
  gt_file_delete(ifp->i2_true_don_fp);
  gt_file_delete(ifp->i1_true_acc_fp);
  gt_file_delete(ifp->i1_true_don_fp);
  gt_file_delete(ifp->i0_true_acc_fp);
  gt_file_delete(ifp->i0_true_don_fp);

  gt_file_delete(ifp->intron_all_fp);
  gt_file_delete(ifp->intron_2_fp);
  gt_file_delete(ifp->intron_1_fp);
  gt_file_delete(ifp->intron_0_fp);
  gt_file_delete(ifp->exon_2_fp);
  gt_file_delete(ifp->exon_1_fp);
  gt_file_delete(ifp->exon_0_fp);

  gt_free(ifp);
}

typedef enum {
  SQUASH,
  TRUE_DON,
  TRUE_ACC,
  FALSE_DON,
  FALSE_ACC
} SeqFileType;

static GtFile* open_seq_file(const char *outdir, const char *filesuffix,
                             bool force, SeqFileType type, GtError *err)
{
  GtFile *file;
  GtStr *path;
  gt_error_check(err);
  gt_assert(filesuffix);
  path = gt_str_new_cstr(outdir);
  gt_str_append_char(path, GT_PATH_SEPARATOR);
  gt_str_append_cstr(path, "file");
  gt_str_append_cstr(path, filesuffix);
  switch (type) {
    case SQUASH:
      break;
    case TRUE_DON:
      gt_str_append_cstr(path, ".GT_AT.truedons");
      break;
    case TRUE_ACC:
      gt_str_append_cstr(path, ".GT_AT.trueaccs");
      break;
     case FALSE_DON:
      gt_str_append_cstr(path, ".GT_AT.falsedons");
      break;
    case FALSE_ACC:
      gt_str_append_cstr(path, ".GT_AT.falseaccs");
      break;
    default:
      gt_assert(0);
  }
  file = gt_output_file_xopen_forcecheck(gt_str_get(path), "w", force,
                                         err);
  gt_str_delete(path);
  return file;
}

IntermediateFiles* intermediate_files_new(GthBSSMSeqProcessor *bsp,
                                          GtError *err)
{
  IntermediateFiles *ifp;

  gt_error_check(err);
  gt_assert(bsp);

  ifp = gt_calloc(1, sizeof *ifp);

  /* squash */
  if (!(ifp->exon_0_fp = open_seq_file(bsp->outdir, "0", bsp->force, SQUASH,
                                       err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->exon_1_fp = open_seq_file(bsp->outdir, "1", bsp->force, SQUASH,
                                       err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->exon_2_fp = open_seq_file(bsp->outdir, "2", bsp->force, SQUASH,
                                       err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->intron_0_fp = open_seq_file(bsp->outdir, "I0", bsp->force, SQUASH,
                                         err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->intron_1_fp = open_seq_file(bsp->outdir, "I1", bsp->force, SQUASH,
                                         err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->intron_2_fp = open_seq_file(bsp->outdir, "I2", bsp->force, SQUASH,
                                         err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->intron_all_fp = open_seq_file(bsp->outdir, "I", bsp->force, SQUASH,
                                           err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }

  /* true donor and acceptor sites */
  if (!(ifp->i0_true_don_fp = open_seq_file(bsp->outdir, "I0", bsp->force,
                                            TRUE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i0_true_acc_fp = open_seq_file(bsp->outdir, "I0", bsp->force,
                                            TRUE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i1_true_don_fp = open_seq_file(bsp->outdir, "I1", bsp->force,
                                            TRUE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i1_true_acc_fp = open_seq_file(bsp->outdir, "I1", bsp->force,
                                            TRUE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i2_true_don_fp = open_seq_file(bsp->outdir, "I2", bsp->force,
                                            TRUE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i2_true_acc_fp = open_seq_file(bsp->outdir, "I2", bsp->force,
                                            TRUE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }

  /* false donor and acceptor sites */
  if (!(ifp->e0_false_don_fp = open_seq_file(bsp->outdir, "E0", bsp->force,
                                             FALSE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->e0_false_acc_fp = open_seq_file(bsp->outdir, "E0", bsp->force,
                                             FALSE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->e1_false_don_fp = open_seq_file(bsp->outdir, "E1", bsp->force,
                                             FALSE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->e1_false_acc_fp = open_seq_file(bsp->outdir, "E1", bsp->force,
                                             FALSE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->e2_false_don_fp = open_seq_file(bsp->outdir, "E2", bsp->force,
                                             FALSE_DON, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->e2_false_acc_fp = open_seq_file(bsp->outdir, "E2", bsp->force,
                                             FALSE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i_false_don_fp = open_seq_file(bsp->outdir, "I", bsp->force,
                                            FALSE_DON,
                                            err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }
  if (!(ifp->i_false_acc_fp = open_seq_file(bsp->outdir, "I", bsp->force,
                                            FALSE_ACC, err))) {
    intermediate_files_delete(ifp);
    return NULL;
  }

  return ifp;
}

typedef struct {
  GtStr *seqid,
        *seq,
        *desc;
  GtRange range;
  bool reverse;
  unsigned int phase;
} BSSMSeq;

static GtStr* construct_seq_description(const GtRange *range, bool reverse,
                                        unsigned int phase, GtStr *seqid)
{
  GtStr *desc;
  gt_assert(range && phase <= 2);
  desc = gt_str_new();
  gt_str_append_ulong(desc, reverse ? range->end : range->start);
  gt_str_append_char(desc, ' ');
  gt_str_append_ulong(desc, reverse ? range->start : range->end);
  gt_str_append_char(desc, ' ');
  gt_str_append_uint(desc, phase);
  gt_str_append_char(desc, ' ');
  gt_str_append_str(desc, seqid);
  gt_str_append_char(desc, reverse ? '-' : '+');
  return desc;
}

static BSSMSeq* bssm_seq_new(GtStr *seqid, const GtRange *range, bool reverse,
                             unsigned int phase, const GtStr *seq)
{
  BSSMSeq *s;
  gt_assert(seqid && range && phase <= 2 && seqid);
  s = gt_malloc(sizeof *s);
  s->seqid = gt_str_ref(seqid);
  s->range = *range;
  s->reverse = reverse;
  s->phase = phase;
  s->seq = gt_str_clone(seq);
  s->desc = construct_seq_description(range, reverse, phase, seqid);
  return s;
}

static void bssm_seq_delete(BSSMSeq *s)
{
  if (!s) return;
  gt_str_delete(s->desc);
  gt_str_delete(s->seq);
  gt_str_delete(s->seqid);
  gt_free(s);
}

static GtStr* bssm_seq_seq(BSSMSeq *s)
{
  gt_assert(s && s->seq);
  return s->seq;
}

static void bssm_seq_append_desc(BSSMSeq *s_1, BSSMSeq *s_2)
{
  gt_assert(s_1 && s_2);
  gt_str_append_cstr(s_1->desc, " AND *** ");
  gt_str_append_str(s_1->desc, s_2->desc);
}

static int bssm_seq_compare(BSSMSeq **seq_1, BSSMSeq **seq_2)
{
  return gt_str_cmp((*seq_1)->seq, (*seq_2)->seq);
}

static void bssm_seq_write(BSSMSeq *s, GtFile *outfp)
{
  gt_assert(s && outfp);
  gt_fasta_show_entry(gt_str_get(s->desc), gt_str_get(s->seq),
                      gt_str_length(s->seq), 0, outfp);
}

static GtFile* open_res_file(const char *outdir, const char *resdir,
                             const char *filename, GtFileMode file_mode,
                             bool force, GtError *err)
{
  GtFile *file;
  GtStr *path;
  gt_error_check(err);
  gt_assert(outdir && resdir && filename);
  path = gt_str_new_cstr(outdir);
  gt_str_append_char(path, GT_PATH_SEPARATOR);
  gt_str_append_cstr(path, resdir);
  gt_str_append_char(path, GT_PATH_SEPARATOR);
  gt_str_append_cstr(path, filename);
  gt_str_append_cstr(path, gt_file_mode_suffix(file_mode));
  file = gt_output_file_xopen_forcecheck(gt_str_get(path), "w", force,
                                         err);
  gt_str_delete(path);
  return file;
}

GthBSSMSeqProcessor* gth_bssm_seq_processor_new(const char *outdir,
                                                GtFileMode fm, bool force,
                                                bool gcdonor, GtError *err)
{
  GthBSSMSeqProcessor *bsp;
  GtStr *dir;
  gt_error_check(err);
  gt_assert(outdir);

  bsp = gt_calloc(1, sizeof *bsp);
  bsp->outdir = gt_cstr_dup(outdir);
  bsp->force = force;
  bsp->gcdonor = gcdonor;

  bsp->exon_0 = gt_array_new(sizeof (BSSMSeq*));
  bsp->exon_1 = gt_array_new(sizeof (BSSMSeq*));
  bsp->exon_2 = gt_array_new(sizeof (BSSMSeq*));
  bsp->intron_0 = gt_array_new(sizeof (BSSMSeq*));
  bsp->intron_1 = gt_array_new(sizeof (BSSMSeq*));
  bsp->intron_2 = gt_array_new(sizeof (BSSMSeq*));
  bsp->intron_all = gt_array_new(sizeof (BSSMSeq*));

  /* creat outdir, if necessary */
  dir = gt_str_new_cstr(bsp->outdir);
  if (!gt_file_exists(gt_str_get(dir)))
    gt_xmkdir(gt_str_get(dir));

  /* create GT directory, if necessary */
  gt_str_append_char(dir, GT_PATH_SEPARATOR);
  gt_str_append_cstr(dir, GT_DIR);
  if (!gt_file_exists(gt_str_get(dir)))
    gt_xmkdir(gt_str_get(dir));

  /* create GC directory, if necessary */
  if (gcdonor) {
    gt_str_set(dir, bsp->outdir);
    gt_str_append_char(dir, GT_PATH_SEPARATOR);
    gt_str_append_cstr(dir, GC_DIR);
    if (!gt_file_exists(gt_str_get(dir)))
      gt_xmkdir(gt_str_get(dir));
  }

  /* create AG directory, if necessary */
  gt_str_set(dir, bsp->outdir);
  gt_str_append_char(dir, GT_PATH_SEPARATOR);
  gt_str_append_cstr(dir, AG_DIR);
  if (!gt_file_exists(gt_str_get(dir)))
    gt_xmkdir(gt_str_get(dir));
  gt_str_delete(dir);

  /* open result files */
  if (!(bsp->gt_t1_fp = open_res_file(outdir, GT_DIR, "T1", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_t2_fp = open_res_file(outdir, GT_DIR, "T2", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_t0_fp = open_res_file(outdir, GT_DIR, "T0", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_f1_fp = open_res_file(outdir, GT_DIR, "F1", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_f2_fp = open_res_file(outdir, GT_DIR, "F2", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_f0_fp = open_res_file(outdir, GT_DIR, "F0", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->gt_fi_fp = open_res_file(outdir, GT_DIR, "Fi", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }

  if (bsp->gcdonor) {
    if (!(bsp->gc_t1_fp = open_res_file(outdir, GC_DIR, "T1", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_t2_fp = open_res_file(outdir, GC_DIR, "T2", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_t0_fp = open_res_file(outdir, GC_DIR, "T0", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_f1_fp = open_res_file(outdir, GC_DIR, "F1", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_f2_fp = open_res_file(outdir, GC_DIR, "F2", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_f0_fp = open_res_file(outdir, GC_DIR, "F0", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
    if (!(bsp->gc_fi_fp = open_res_file(outdir, GC_DIR, "Fi", fm, force,
                                        err))) {
      gth_bssm_seq_processor_delete(bsp);
      return NULL;
    }
  }

  if (!(bsp->ag_t1_fp = open_res_file(outdir, AG_DIR, "T1", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_t2_fp = open_res_file(outdir, AG_DIR, "T2", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_t0_fp = open_res_file(outdir, AG_DIR, "T0", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_f1_fp = open_res_file(outdir, AG_DIR, "F1", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_f2_fp = open_res_file(outdir, AG_DIR, "F2", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_f0_fp = open_res_file(outdir, AG_DIR, "F0", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }
  if (!(bsp->ag_fi_fp = open_res_file(outdir, AG_DIR, "Fi", fm, force, err))) {
    gth_bssm_seq_processor_delete(bsp);
    return NULL;
  }

  return bsp;
}

static void bssm_seqs_delete(GtArray *seqs)
{
  GtUword i;
  for (i = 0; i < gt_array_size(seqs); i++)
    bssm_seq_delete(*(BSSMSeq**) gt_array_get(seqs, i));
  gt_array_delete(seqs);
}

void gth_bssm_seq_processor_delete(GthBSSMSeqProcessor *bsp)
{
  if (!bsp) return;

  gt_file_delete(bsp->ag_fi_fp);
  gt_file_delete(bsp->ag_f0_fp);
  gt_file_delete(bsp->ag_f2_fp);
  gt_file_delete(bsp->ag_f1_fp);
  gt_file_delete(bsp->ag_t0_fp);
  gt_file_delete(bsp->ag_t2_fp);
  gt_file_delete(bsp->ag_t1_fp);

  gt_file_delete(bsp->gc_fi_fp);
  gt_file_delete(bsp->gc_f0_fp);
  gt_file_delete(bsp->gc_f2_fp);
  gt_file_delete(bsp->gc_f1_fp);
  gt_file_delete(bsp->gc_t0_fp);
  gt_file_delete(bsp->gc_t2_fp);
  gt_file_delete(bsp->gc_t1_fp);

  gt_file_delete(bsp->gt_fi_fp);
  gt_file_delete(bsp->gt_f0_fp);
  gt_file_delete(bsp->gt_f2_fp);
  gt_file_delete(bsp->gt_f1_fp);
  gt_file_delete(bsp->gt_t0_fp);
  gt_file_delete(bsp->gt_t2_fp);
  gt_file_delete(bsp->gt_t1_fp);

  bssm_seqs_delete(bsp->i_false_acc);
  bssm_seqs_delete(bsp->i_false_don_gc);
  bssm_seqs_delete(bsp->i_false_don_gt);
  bssm_seqs_delete(bsp->e2_false_acc);
  bssm_seqs_delete(bsp->e2_false_don_gc);
  bssm_seqs_delete(bsp->e2_false_don_gt);
  bssm_seqs_delete(bsp->e1_false_acc);
  bssm_seqs_delete(bsp->e1_false_don_gc);
  bssm_seqs_delete(bsp->e1_false_don_gt);
  bssm_seqs_delete(bsp->e0_false_acc);
  bssm_seqs_delete(bsp->e0_false_don_gc);
  bssm_seqs_delete(bsp->e0_false_don_gt);

  bssm_seqs_delete(bsp->i2_true_acc);
  bssm_seqs_delete(bsp->i2_true_don_gc);
  bssm_seqs_delete(bsp->i2_true_don_gt);
  bssm_seqs_delete(bsp->i1_true_acc);
  bssm_seqs_delete(bsp->i1_true_don_gc);
  bssm_seqs_delete(bsp->i1_true_don_gt);
  bssm_seqs_delete(bsp->i0_true_acc);
  bssm_seqs_delete(bsp->i0_true_don_gc);
  bssm_seqs_delete(bsp->i0_true_don_gt);

  bssm_seqs_delete(bsp->intron_all);
  bssm_seqs_delete(bsp->intron_2);
  bssm_seqs_delete(bsp->intron_1);
  bssm_seqs_delete(bsp->intron_0);
  bssm_seqs_delete(bsp->exon_2);
  bssm_seqs_delete(bsp->exon_1);
  bssm_seqs_delete(bsp->exon_0);

  gt_free(bsp->outdir);
  gt_free(bsp);
}

void gth_bssm_seq_processor_proc_exon(GthBSSMSeqProcessor *bsp,
                                      unsigned int phase, GtStr *seqid,
                                      const GtRange *range, bool reverse,
                                      const GtStr *seq)
{
  GtArray *bssm_seqs = NULL;
  BSSMSeq *bssm_seq;
  gt_assert(bsp);
  gt_assert(phase <= 2);
  gt_assert(seq && gt_str_length(seq) == gt_range_length(range));
  switch (phase) {
    case 0: bssm_seqs = bsp->exon_0; break;
    case 1: bssm_seqs = bsp->exon_1; break;
    case 2: bssm_seqs = bsp->exon_2; break;
    default: gt_assert(0);
  }
  bssm_seq = bssm_seq_new(seqid, range, reverse, phase, seq);
  gt_array_add(bssm_seqs, bssm_seq);
}

void gth_bssm_seq_processor_proc_intron(GthBSSMSeqProcessor *bsp,
                                        unsigned int phase, GtStr *seqid,
                                        const GtRange *range, bool reverse,
                                        const GtStr *seq)
{
  GtArray *bssm_seqs = NULL;
  BSSMSeq *bssm_seq;
  gt_assert(bsp);
  gt_assert(phase <= 2);
  gt_assert(seq && gt_str_length(seq) == gt_range_length(range));
  switch (phase) {
    case 0: bssm_seqs = bsp->intron_0; break;
    case 1: bssm_seqs = bsp->intron_1; break;
    case 2: bssm_seqs = bsp->intron_2; break;
    default: gt_assert(0);
  }
  bssm_seq = bssm_seq_new(seqid, range, reverse, phase, seq);
  gt_array_add(bssm_seqs, bssm_seq);
  bssm_seq = bssm_seq_new(seqid, range, reverse, phase, seq);
  gt_array_add(bsp->intron_all, bssm_seq);
}

void bssm_seqs_squash(GtArray *seqs)
{
  BSSMSeq *cur_seq, *pre_seq;
  GtArray *tmp_seqs;
  GtUword i;
  gt_assert(seqs);
  if (!gt_array_size(seqs))
    return;
  /* sort sequences */
  qsort(gt_array_get_space(seqs), gt_array_size(seqs), gt_array_elem_size(seqs),
        (GtCompare) bssm_seq_compare);
  tmp_seqs = gt_array_new(sizeof (BSSMSeq*));
  pre_seq = *(BSSMSeq**) gt_array_get_first(seqs);
  for (i = 1; i < gt_array_size(seqs); i++) {
    cur_seq = *(BSSMSeq**) gt_array_get(seqs, i);
    if (!gt_str_cmp(bssm_seq_seq(pre_seq), bssm_seq_seq(cur_seq))) {
      /* sequences are equal */
      bssm_seq_append_desc(pre_seq, cur_seq);
      bssm_seq_delete(cur_seq);
    }
    else {
      gt_array_add(tmp_seqs, pre_seq);
      pre_seq = cur_seq;
    }
  }
  gt_array_add(tmp_seqs, pre_seq);
  gt_array_reset(seqs);
  gt_array_add_array(seqs, tmp_seqs);
  gt_array_delete(tmp_seqs);
}

void gth_bssm_seq_processor_squash(GthBSSMSeqProcessor *bsp)
{
  gt_assert(bsp);
  bssm_seqs_squash(bsp->exon_0);
  bssm_seqs_squash(bsp->exon_1);
  bssm_seqs_squash(bsp->exon_2);
  bssm_seqs_squash(bsp->intron_0);
  bssm_seqs_squash(bsp->intron_1);
  bssm_seqs_squash(bsp->intron_2);
  bssm_seqs_squash(bsp->intron_all);
}

static int get_true_seq(GtArray *true_sites, BSSMSeq *intron,
                        const char *sequence,
                        GT_UNUSED GtUword sequence_length,
                        GtStr *seq, const GtRange *site_range, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(true_sites && intron);
  /* 1-based coordinates */
  gt_assert(intron->range.start && intron->range.end);
  gt_assert(intron->range.end <= sequence_length);
  gt_str_reset(seq);
  gt_str_append_cstr_nt(seq, sequence, gt_range_length(site_range));
  if (intron->reverse)
    had_err = gt_reverse_complement(gt_str_get(seq), gt_str_length(seq), err);
  if (!had_err && !gth_seq_contains_wildcard(seq)) {
    BSSMSeq *true_seq = bssm_seq_new(intron->seqid, &intron->range,
                                     intron->reverse, intron->phase, seq);
    gt_array_add(true_sites, true_seq);
  }
  return had_err;
}

static int find_true_sites(GtArray *true_don_gt, GtArray *true_don_gc,
                           GtArray *true_acc, GtArray *introns,
                           GtRegionMapping *region_mapping, bool gcdonor,
                           GtError *err)
{
  GtUword i, len;
  bool don_underflow, acc_underflow;
  const char *cseq;
  GtRange don_range, acc_range;
  BSSMSeq *intron;
  GtStr *seq;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(true_don_gt && true_acc && introns && region_mapping);
  seq = gt_str_new();
  for (i = 0; !had_err && i < gt_array_size(introns); i++) {
    intron = *(BSSMSeq**) gt_array_get(introns, i);
    len = gt_str_length(intron->seq);
    if (len >= 4) {
      cseq = gt_str_get(intron->seq);
      if (  (cseq[0]     == 'G' || cseq[0]     == 'g') &&
           ((cseq[1]     == 'T' || cseq[1]     == 't') ||
            (gcdonor &&
            (cseq[1]     == 'C' || cseq[1]     == 'c'))) &&
            (cseq[len-2] == 'A' || cseq[len-2] == 'a') &&
            (cseq[len-1] == 'G' || cseq[len-1] == 'g')) {
        GtUword sequence_length;
        char *sequence = NULL;
        /* correct splice site found -> get flanking sequences */
        if (!intron->reverse) {
          if (intron->range.start >= 50) {
            don_range.start = intron->range.start - 50;
            don_underflow = false;
          }
          else
            don_underflow = true;
          don_range.end   = intron->range.start + 51;
          if (intron->range.end >= 51) {
            acc_range.start = intron->range.end - 51;
            acc_underflow = false;
          }
          else
            acc_underflow = true;
          acc_range.end = intron->range.end + 50;
        }
        else {
          if (intron->range.end >= 51) {
            don_range.start = intron->range.end - 51;
            don_underflow = false;
          }
          else
            don_underflow = true;
          don_range.end = intron->range.end + 50;
          if (intron->range.start >= 50) {
            acc_range.start = intron->range.start - 50;
            acc_underflow = false;
          }
          else
            acc_underflow = true;
          acc_range.end = intron->range.start + 51;
        }
        /* donor sequence */
        if (!had_err && !don_underflow) {
          had_err = gt_region_mapping_get_sequence_length(region_mapping,
                                                          &sequence_length,
                                                          intron->seqid, err);
          if (!had_err && don_range.end <= sequence_length) {
            had_err = gt_region_mapping_get_sequence(region_mapping, &sequence,
                                                     intron->seqid,
                                                     don_range.start,
                                                     don_range.end, err);
            if (!had_err) {
              if (cseq[1] == 'T' || cseq[1] == 't') {
                had_err = get_true_seq(true_don_gt, intron, sequence,
                                       sequence_length, seq, &don_range, err);
              }
              else {
                gt_assert(gcdonor && (cseq[1] == 'C' || cseq[1] == 'c'));
                had_err = get_true_seq(true_don_gc, intron, sequence,
                                       sequence_length, seq, &don_range, err);
              }
              gt_free(sequence);
            }
          }
        }
        /* acceptor sequence */
        if (!had_err && !acc_underflow) {
          had_err = gt_region_mapping_get_sequence_length(region_mapping,
                                                          &sequence_length,
                                                          intron->seqid, err);
          if (!had_err && acc_range.end <= sequence_length) {
            had_err = gt_region_mapping_get_sequence(region_mapping, &sequence,
                                                     intron->seqid,
                                                     acc_range.start,
                                                     acc_range.end, err);
            if (!had_err) {
              had_err = get_true_seq(true_acc, intron, sequence,
                                     sequence_length, seq, &acc_range, err);
              gt_free(sequence);
            }
          }
        }
      }
    }
  }
  gt_str_delete(seq);
  return had_err;
}

int gth_bssm_seq_processor_find_true_sites(GthBSSMSeqProcessor *bsp,
                                           GtRegionMapping *region_mapping,
                                           GtError *err)
{
  int had_err;

  gt_error_check(err);
  gt_assert(bsp && region_mapping);
  gt_assert(!bsp->i0_true_don_gt);
  gt_assert(!bsp->i0_true_don_gc);
  gt_assert(!bsp->i0_true_acc);
  gt_assert(!bsp->i1_true_don_gt);
  gt_assert(!bsp->i1_true_don_gc);
  gt_assert(!bsp->i1_true_acc);
  gt_assert(!bsp->i2_true_don_gt);
  gt_assert(!bsp->i2_true_don_gc);
  gt_assert(!bsp->i2_true_acc);

  bsp->i0_true_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->i0_true_acc = gt_array_new(sizeof (BSSMSeq*));
  bsp->i1_true_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->i1_true_acc = gt_array_new(sizeof (BSSMSeq*));
  bsp->i2_true_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->i2_true_acc = gt_array_new(sizeof (BSSMSeq*));

  if (bsp->gcdonor) {
    bsp->i0_true_don_gc = gt_array_new(sizeof (BSSMSeq*));
    bsp->i1_true_don_gc = gt_array_new(sizeof (BSSMSeq*));
    bsp->i2_true_don_gc = gt_array_new(sizeof (BSSMSeq*));
  }

  had_err = find_true_sites(bsp->i0_true_don_gt, bsp->i0_true_don_gc,
                            bsp->i0_true_acc, bsp->intron_0, region_mapping,
                            bsp->gcdonor, err);
  if (!had_err) {
    had_err = find_true_sites(bsp->i1_true_don_gt, bsp->i1_true_don_gc,
                              bsp->i1_true_acc, bsp->intron_1, region_mapping,
                              bsp->gcdonor, err);
  }
  if (!had_err) {
    had_err = find_true_sites(bsp->i2_true_don_gt, bsp->i2_true_don_gc,
                              bsp->i2_true_acc, bsp->intron_2, region_mapping,
                              bsp->gcdonor, err);
  }

  return had_err;
}

static int get_false_don_seq(GtArray *false_don_0_gt, GtArray *false_don_0_gc,
                             GtArray *false_don_1_gt, GtArray *false_don_1_gc,
                             GtArray *false_don_2_gt, GtArray *false_don_2_gc,
                             BSSMSeq *intron, const char *sequence,
                             GT_UNUSED GtUword sequence_length,
                             GtStr *seq, const GtRange *don_range,
                             bool proc_exons, GT_UNUSED bool gcdonor,
                             GtUword j, GtError *err)
{
  unsigned int phase = 0;
  BSSMSeq *false_seq;
  int had_err = 0;
  gt_error_check(err);

  /* 1-based coordinates */
  gt_assert(intron->range.start && intron->range.end);
  gt_assert(intron->range.end <= sequence_length);
  gt_str_reset(seq);
  gt_str_append_cstr_nt(seq, sequence, gt_range_length(don_range));
  if (intron->reverse)
    had_err = gt_reverse_complement(gt_str_get(seq), gt_str_length(seq), err);
  if (!had_err && !gth_seq_contains_wildcard(seq)) {
    const char *iseq = gt_str_get(seq);
    gt_assert(gt_str_length(seq) == 102);
    gt_assert(iseq[50] == 'G' || iseq[50] == 'g');
    gt_assert(iseq[51] == 'T' || iseq[51] == 't' ||
              (gcdonor && (iseq[51] == 'C' || iseq[51] == 'c')));
    false_seq = bssm_seq_new(intron->seqid, &intron->range, intron->reverse,
                             intron->phase, seq);
    if (proc_exons)
      phase = (intron->phase + j) % 3;
    if (iseq[51] == 'T' || iseq[51] == 't') {
      switch (phase) {
        case 0: gt_array_add(false_don_0_gt, false_seq); break;
        case 1: gt_array_add(false_don_1_gt, false_seq); break;
        case 2: gt_array_add(false_don_2_gt, false_seq); break;
        default: gt_assert(0);
      }
    }
    else {
      gt_assert(gcdonor && (iseq[51] == 'C' || iseq[51] == 'c'));
      switch (phase) {
        case 0: gt_array_add(false_don_0_gc, false_seq); break;
        case 1: gt_array_add(false_don_1_gc, false_seq); break;
        case 2: gt_array_add(false_don_2_gc, false_seq); break;
        default: gt_assert(0);
      }
    }
  }
  return had_err;
}

static int get_false_acc_seq(GtArray *false_acc_0, GtArray *false_acc_1,
                             GtArray *false_acc_2, BSSMSeq *intron,
                             const char *sequence,
                             GT_UNUSED GtUword sequence_length,
                             GtStr *seq, const GtRange *acc_range,
                             bool proc_exons, GtUword j, GtError *err)
{
  unsigned int phase = 0;
  BSSMSeq *false_seq;
  int had_err = 0;
  gt_error_check(err);

  /* 1-based coordinates */
  gt_assert(intron->range.start && intron->range.end);
  gt_assert(intron->range.end <= sequence_length);
  gt_str_reset(seq);
  gt_str_append_cstr_nt(seq, sequence, gt_range_length(acc_range));
  if (intron->reverse)
    had_err = gt_reverse_complement(gt_str_get(seq), gt_str_length(seq), err);
  if (!had_err && !gth_seq_contains_wildcard(seq)) {
#ifndef NDEBUG
    const char *iseq = gt_str_get(seq);
#endif
    gt_assert(gt_str_length(seq) == 102);
    gt_assert(iseq[50] == 'A' || iseq[50] == 'a');
    gt_assert(iseq[51] == 'G' || iseq[51] == 'g');
    false_seq = bssm_seq_new(intron->seqid, &intron->range, intron->reverse,
                             intron->phase, seq);
    if (proc_exons)
      phase = (intron->phase + j) % 3;
    switch (phase) {
      case 0: gt_array_add(false_acc_0, false_seq); break;
      case 1: gt_array_add(false_acc_1, false_seq); break;
      case 2: gt_array_add(false_acc_2, false_seq); break;
      default: gt_assert(0);
    }
  }
  return had_err;
}

static int find_false_sites(GtArray *false_don_0_gt, GtArray *false_don_0_gc,
                            GtArray *false_acc_0,
                            GtArray *false_don_1_gt, GtArray *false_don_1_gc,
                            GtArray *false_acc_1,
                            GtArray *false_don_2_gt, GtArray *false_don_2_gc,
                            GtArray *false_acc_2,
                            GtArray *seqs, bool proc_exons,
                            GtRegionMapping *region_mapping, bool gcdonor,
                            GtError *err)
{
  GtUword i, j, len;
  bool don_underflow, acc_underflow;
  const char *cseq;
  GtRange don_range = {0}, acc_range = {0};
  BSSMSeq *intron;
  GtStr *seq;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(false_don_0_gt && false_acc_0 && seqs && region_mapping && err);
  seq = gt_str_new();
  for (i = 0; !had_err && i < gt_array_size(seqs); i++) {
    intron = *(BSSMSeq**) gt_array_get(seqs, i);
    len = gt_str_length(intron->seq);
    if (len >= 2) {
      cseq = gt_str_get(intron->seq);
      for (j = 0; !had_err && j < len - 1; j++) {
        GtUword sequence_length;
        char *sequence = NULL;
        if ((proc_exons || j) && /* skip true donor site */
             (cseq[j]   == 'G' || cseq[j]   == 'g') &&
            ((cseq[j+1] == 'T' || cseq[j+1] == 't') ||
             (gcdonor &&
             (cseq[j+1] == 'C' || cseq[j+1] == 'c')))) {
          if (!intron->reverse) {
            if (intron->range.start + j >= 50) {
              don_range.start = intron->range.start + j - 50;
              don_underflow = false;
            }
            else
              don_underflow = true;
            don_range.end = intron->range.start + j + 51;
          }
          else {
            if (intron->range.end >= j + 51) {
              don_range.start = intron->range.end - j - 51;
              don_underflow = false;
            }
            else
              don_underflow = true;
            don_range.end = intron->range.end - j + 50;
          }
          /* donor sequence */
          if (!had_err && !don_underflow) {
            had_err = gt_region_mapping_get_sequence_length(region_mapping,
                                                            &sequence_length,
                                                            intron->seqid, err);
            if (!had_err && don_range.end < sequence_length) {
              had_err = gt_region_mapping_get_sequence(region_mapping,
                                                       &sequence, intron->seqid,
                                                       don_range.start,
                                                       don_range.end, err);
              if (!had_err) {
               had_err = get_false_don_seq(false_don_0_gt, false_don_0_gc,
                                           false_don_1_gt, false_don_1_gc,
                                           false_don_2_gt, false_don_2_gc,
                                           intron, sequence, sequence_length,
                                           seq, &don_range, proc_exons, gcdonor,
                                           j, err);
              }
            }
          }
        }
        else  if ((proc_exons || j < len - 2) && /* skip true acceptor sites */
                  (cseq[j]   == 'A' || cseq[j]   == 'a') &&
                  (cseq[j+1] == 'G' || cseq[j+1] == 'g')) {
          if (!intron->reverse) {
            if (intron->range.start + j >= 50) {
              acc_range.start = intron->range.start + j - 50;
              acc_underflow = false;
            }
            else
              acc_underflow = true;
            acc_range.end = intron->range.start + j + 51;
          }
          else {
            if (intron->range.end >= j + 51) {
              acc_range.start = intron->range.end - j - 51;
              acc_underflow = false;
            }
            else
              acc_underflow = true;
            acc_range.end = intron->range.end - j + 50;
          }
          /* acceptor sequence */
          if (!had_err && !acc_underflow) {
            had_err = gt_region_mapping_get_sequence_length(region_mapping,
                                                            &sequence_length,
                                                            intron->seqid, err);
            if (!had_err && acc_range.end < sequence_length) {
              had_err = gt_region_mapping_get_sequence(region_mapping,
                                                       &sequence, intron->seqid,
                                                       acc_range.start,
                                                       acc_range.end, err);
              if (!had_err) {
                had_err = get_false_acc_seq(false_acc_0, false_acc_1,
                                            false_acc_2, intron, sequence,
                                            sequence_length, seq, &acc_range,
                                            proc_exons, j, err);
              }
            }
          }
        }
        gt_free(sequence);
      }
    }
  }
  gt_str_delete(seq);
  return had_err;
}

int gth_bssm_seq_processor_find_false_sites(GthBSSMSeqProcessor *bsp,
                                            GtRegionMapping *region_mapping,
                                            GtError *err)
{
  int had_err = 0;

  gt_error_check(err);
  gt_assert(bsp && region_mapping);
  gt_assert(!bsp->e0_false_don_gt);
  gt_assert(!bsp->e0_false_don_gc);
  gt_assert(!bsp->e0_false_acc);
  gt_assert(!bsp->e1_false_don_gt);
  gt_assert(!bsp->e1_false_don_gc);
  gt_assert(!bsp->e1_false_acc);
  gt_assert(!bsp->e2_false_don_gt);
  gt_assert(!bsp->e2_false_don_gc);
  gt_assert(!bsp->e2_false_acc);
  gt_assert(!bsp->i_false_don_gt);
  gt_assert(!bsp->i_false_don_gc);
  gt_assert(!bsp->i_false_acc);

  bsp->e0_false_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->e0_false_acc = gt_array_new(sizeof (BSSMSeq*));
  bsp->e1_false_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->e1_false_acc = gt_array_new(sizeof (BSSMSeq*));
  bsp->e2_false_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->e2_false_acc = gt_array_new(sizeof (BSSMSeq*));
  bsp->i_false_don_gt = gt_array_new(sizeof (BSSMSeq*));
  bsp->i_false_acc = gt_array_new(sizeof (BSSMSeq*));

  if (bsp->gcdonor) {
    bsp->e0_false_don_gc = gt_array_new(sizeof (BSSMSeq*));
    bsp->e1_false_don_gc = gt_array_new(sizeof (BSSMSeq*));
    bsp->e2_false_don_gc = gt_array_new(sizeof (BSSMSeq*));
    bsp->i_false_don_gc = gt_array_new(sizeof (BSSMSeq*));
  }

  had_err = find_false_sites(bsp->i_false_don_gt, bsp->i_false_don_gc,
                             bsp->i_false_acc, NULL, NULL, NULL, NULL, NULL,
                             NULL, bsp->intron_all, false, region_mapping,
                             bsp->gcdonor, err);

  if (!had_err) {
    had_err = find_false_sites(bsp->e0_false_don_gt, bsp->e0_false_don_gc,
                               bsp->e0_false_acc,
                               bsp->e1_false_don_gt, bsp->e1_false_don_gc,
                               bsp->e1_false_acc,
                               bsp->e2_false_don_gt, bsp->e2_false_don_gc,
                               bsp->e2_false_acc, bsp->exon_0, true,
                               region_mapping, bsp->gcdonor, err);
  }

  if (!had_err) {
    had_err = find_false_sites(bsp->e0_false_don_gt, bsp->e0_false_don_gc,
                               bsp->e0_false_acc,
                               bsp->e1_false_don_gt, bsp->e1_false_don_gc,
                               bsp->e1_false_acc,
                               bsp->e2_false_don_gt, bsp->e2_false_don_gc,
                               bsp->e2_false_acc, bsp->exon_1, true,
                               region_mapping, bsp->gcdonor, err);
  }

  if (!had_err) {
    had_err = find_false_sites(bsp->e0_false_don_gt, bsp->e0_false_don_gc,
                               bsp->e0_false_acc,
                               bsp->e1_false_don_gt, bsp->e1_false_don_gc,
                               bsp->e1_false_acc,
                               bsp->e2_false_don_gt, bsp->e2_false_don_gc,
                               bsp->e2_false_acc, bsp->exon_2, true,
                               region_mapping, bsp->gcdonor, err);
  }

  return had_err;
}

static void bssm_seqs_write(GtArray *seqs, GtFile *outfp)
{
  GtUword i;
  gt_assert(seqs && outfp);
  for (i = 0; i < gt_array_size(seqs); i++)
    bssm_seq_write(*(BSSMSeq**) gt_array_get(seqs, i), outfp);
}

static void sample_bssm_seqs(GtArray *seqs, GtUword target_size)
{
  GtUword i, randnum, original_size, sample_size = 0;
  GtBittab *samples;
  GtArray *tmp_seqs;

  gt_assert(seqs);
  gt_assert(target_size <= gt_array_size(seqs));

  original_size = gt_array_size(seqs);
  if (!original_size)
    return; /* nothing to do */
  tmp_seqs = gt_array_new(sizeof (BSSMSeq*));
  samples = gt_bittab_new(original_size);

  /* compute random samples */
  while (sample_size < target_size) {
    randnum = gt_rand_max(original_size - 1);
    if (!gt_bittab_bit_is_set(samples, randnum)) {
      gt_bittab_set_bit(samples, randnum);
      sample_size++;
    }
  }

  gt_assert(gt_bittab_count_set_bits(samples) == target_size);

  /* process samples */
  for (i = 0; i < gt_bittab_size(samples); i++) {
    if (gt_bittab_bit_is_set(samples, i)) {
      BSSMSeq *seq = (*(BSSMSeq**) gt_array_get(seqs, i));
      gt_array_add(tmp_seqs, seq);
    }
    else
      bssm_seq_delete(*(BSSMSeq**) gt_array_get(seqs, i));
  }
  gt_array_reset(seqs);
  gt_array_add_array(seqs, tmp_seqs);

  gt_bittab_delete(samples);
  gt_array_delete(tmp_seqs);
}

int gth_bssm_seq_processor_write_intermediate(GthBSSMSeqProcessor *bsp,
                                              GtError *err)
{
  IntermediateFiles *ifp;

  gt_error_check(err);
  gt_assert(bsp);

  if (!(ifp = intermediate_files_new(bsp, err)))
    return -1;

  bssm_seqs_write(bsp->exon_0, ifp->exon_0_fp);
  bssm_seqs_write(bsp->exon_1, ifp->exon_1_fp);
  bssm_seqs_write(bsp->exon_2, ifp->exon_2_fp);
  bssm_seqs_write(bsp->intron_0, ifp->intron_0_fp);
  bssm_seqs_write(bsp->intron_1, ifp->intron_1_fp);
  bssm_seqs_write(bsp->intron_2, ifp->intron_2_fp);
  bssm_seqs_write(bsp->intron_all, ifp->intron_all_fp);

  bssm_seqs_write(bsp->i0_true_don_gt, ifp->i0_true_don_fp);
  bssm_seqs_write(bsp->i0_true_acc, ifp->i0_true_acc_fp);
  bssm_seqs_write(bsp->i1_true_don_gt, ifp->i1_true_don_fp);
  bssm_seqs_write(bsp->i1_true_acc, ifp->i1_true_acc_fp);
  bssm_seqs_write(bsp->i2_true_don_gt, ifp->i2_true_don_fp);
  bssm_seqs_write(bsp->i2_true_acc, ifp->i2_true_acc_fp);

  bssm_seqs_write(bsp->e0_false_don_gt, ifp->e0_false_don_fp);
  bssm_seqs_write(bsp->e0_false_acc, ifp->e0_false_acc_fp);
  bssm_seqs_write(bsp->e1_false_don_gt, ifp->e1_false_don_fp);
  bssm_seqs_write(bsp->e1_false_acc, ifp->e1_false_acc_fp);
  bssm_seqs_write(bsp->e2_false_don_gt, ifp->e2_false_don_fp);
  bssm_seqs_write(bsp->e2_false_acc, ifp->e2_false_acc_fp);
  bssm_seqs_write(bsp->i_false_don_gt, ifp->i_false_don_fp);
  bssm_seqs_write(bsp->i_false_acc, ifp->i_false_acc_fp);

  intermediate_files_delete(ifp);

  return 0;
}

#if 0
  "T1",
  "T2",
  "T0",
  "F1",
  "F2",
  "F0",
  "Fi"
#endif

static void show_sample_sizes(GthBSSMSeqProcessor *bsp, bool verbose,
                              GtFile *logfp)
{
  GtUword len0, len1, len2;

  len0 = gt_array_size(bsp->i0_true_don_gt);
  len1 = gt_array_size(bsp->i1_true_don_gt);
  len2 = gt_array_size(bsp->i2_true_don_gt);

  if (verbose) {
    printf("%s/T1: "GT_WU" seqs\n", GT_DIR, len0);
    printf("%s/T2: "GT_WU" seqs\n", GT_DIR, len1);
    printf("%s/T0: "GT_WU" seqs\n", GT_DIR, len2);
    printf("%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n", GT_DIR, len0,
           gt_array_size(bsp->e0_false_don_gt));
    printf("%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n", GT_DIR, len1,
           gt_array_size(bsp->e1_false_don_gt));
    printf("%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n", GT_DIR, len2,
           gt_array_size(bsp->e2_false_don_gt));
    printf("%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n", GT_DIR,
           MAX3(len0, len1, len2),
           gt_array_size(bsp->i_false_don_gt));
  }

  gt_file_xprintf(logfp, "%s/T1: "GT_WU" seqs\n", GT_DIR, len0);
  gt_file_xprintf(logfp, "%s/T2: "GT_WU" seqs\n", GT_DIR, len1);
  gt_file_xprintf(logfp, "%s/T0: "GT_WU" seqs\n", GT_DIR, len2);
  gt_file_xprintf(logfp, "%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  GT_DIR, len0, gt_array_size(bsp->e0_false_don_gt));
  gt_file_xprintf(logfp, "%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  GT_DIR, len1, gt_array_size(bsp->e1_false_don_gt));
  gt_file_xprintf(logfp, "%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  GT_DIR, len2, gt_array_size(bsp->e2_false_don_gt));
  gt_file_xprintf(logfp, "%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  GT_DIR, MAX3(len0, len1, len2),
                  gt_array_size(bsp->i_false_don_gt));

  if (bsp->gcdonor) {
    len0 = gt_array_size(bsp->i0_true_don_gc);
    len1 = gt_array_size(bsp->i1_true_don_gc);
    len2 = gt_array_size(bsp->i2_true_don_gc);

    if (verbose) {
      printf("%s/T1: "GT_WU" seqs\n", GC_DIR, len0);
      printf("%s/T2: "GT_WU" seqs\n", GC_DIR, len1);
      printf("%s/T0: "GT_WU" seqs\n", GC_DIR, len2);
      printf("%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n", GC_DIR, len0,
             gt_array_size(bsp->e0_false_don_gc));
      printf("%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n", GC_DIR, len1,
             gt_array_size(bsp->e1_false_don_gc));
      printf("%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n", GC_DIR, len2,
             gt_array_size(bsp->e2_false_don_gc));
      printf("%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n", GC_DIR,
             MAX3(len0, len1, len2),
             gt_array_size(bsp->i_false_don_gc));
    }

    gt_file_xprintf(logfp, "%s/T1: "GT_WU" seqs\n", GC_DIR, len0);
    gt_file_xprintf(logfp, "%s/T2: "GT_WU" seqs\n", GC_DIR, len1);
    gt_file_xprintf(logfp, "%s/T0: "GT_WU" seqs\n", GC_DIR, len2);
    gt_file_xprintf(logfp, "%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n",
                    GC_DIR, len0, gt_array_size(bsp->e0_false_don_gc));
    gt_file_xprintf(logfp, "%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n",
                    GC_DIR, len1, gt_array_size(bsp->e1_false_don_gc));
    gt_file_xprintf(logfp, "%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n",
                    GC_DIR, len2, gt_array_size(bsp->e2_false_don_gc));
    gt_file_xprintf(logfp, "%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n",
                    GC_DIR, MAX3(len0, len1, len2),
                    gt_array_size(bsp->i_false_don_gc));
  }

  len0 = gt_array_size(bsp->i0_true_acc);
  len1 = gt_array_size(bsp->i1_true_acc);
  len2 = gt_array_size(bsp->i2_true_acc);

  if (verbose) {
    printf("%s/T1: "GT_WU" seqs\n", AG_DIR, len0);
    printf("%s/T2: "GT_WU" seqs\n", AG_DIR, len1);
    printf("%s/T0: "GT_WU" seqs\n", AG_DIR, len2);
    printf("%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n", AG_DIR, len0,
           gt_array_size(bsp->e0_false_acc));
    printf("%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n", AG_DIR, len1,
           gt_array_size(bsp->e1_false_acc));
    printf("%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n", AG_DIR, len2,
           gt_array_size(bsp->e2_false_acc));
    printf("%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n", AG_DIR,
           MAX3(len0, len1, len2),
           gt_array_size(bsp->i_false_acc));
  }

  gt_file_xprintf(logfp, "%s/T1: "GT_WU" seqs\n", AG_DIR, len0);
  gt_file_xprintf(logfp, "%s/T2: "GT_WU" seqs\n", AG_DIR, len1);
  gt_file_xprintf(logfp, "%s/T0: "GT_WU" seqs\n", AG_DIR, len2);
  gt_file_xprintf(logfp, "%s/F1: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  AG_DIR, len0, gt_array_size(bsp->e0_false_acc));
  gt_file_xprintf(logfp, "%s/F2: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  AG_DIR, len1, gt_array_size(bsp->e1_false_acc));
  gt_file_xprintf(logfp, "%s/F0: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  AG_DIR, len2, gt_array_size(bsp->e2_false_acc));
  gt_file_xprintf(logfp, "%s/Fi: "GT_WU" seqs (sampled out of "GT_WU")\n",
                  AG_DIR, MAX3(len0, len1, len2),
                  gt_array_size(bsp->i_false_acc));

}

void gth_bssm_seq_processor_sample(GthBSSMSeqProcessor *bsp, bool verbose,
                                   GtFile *logfp)
{
  GtUword len0, len1, len2;
  gt_assert(bsp);

  show_sample_sizes(bsp, verbose, logfp);

  sample_bssm_seqs(bsp->e0_false_don_gt, gt_array_size(bsp->i0_true_don_gt));
  sample_bssm_seqs(bsp->e0_false_acc, gt_array_size(bsp->i0_true_acc));
  sample_bssm_seqs(bsp->e1_false_don_gt, gt_array_size(bsp->i1_true_don_gt));
  sample_bssm_seqs(bsp->e1_false_acc, gt_array_size(bsp->i1_true_acc));
  sample_bssm_seqs(bsp->e2_false_don_gt, gt_array_size(bsp->i2_true_don_gt));
  sample_bssm_seqs(bsp->e2_false_acc, gt_array_size(bsp->i2_true_acc));

  if (bsp->gcdonor) {
    sample_bssm_seqs(bsp->e0_false_don_gc, gt_array_size(bsp->i0_true_don_gc));
    sample_bssm_seqs(bsp->e1_false_don_gc, gt_array_size(bsp->i1_true_don_gc));
    sample_bssm_seqs(bsp->e2_false_don_gc, gt_array_size(bsp->i2_true_don_gc));
  }

  len0 = gt_array_size(bsp->i0_true_don_gt);
  len1 = gt_array_size(bsp->i1_true_don_gt);
  len2 = gt_array_size(bsp->i2_true_don_gt);
  sample_bssm_seqs(bsp->i_false_don_gt, MAX3(len0, len1, len2));

  if (bsp->gcdonor) {
    len0 = gt_array_size(bsp->i0_true_don_gc);
    len1 = gt_array_size(bsp->i1_true_don_gc);
    len2 = gt_array_size(bsp->i2_true_don_gc);
    sample_bssm_seqs(bsp->i_false_don_gc, MAX3(len0, len1, len2));
  }

  len0 = gt_array_size(bsp->i0_true_acc);
  len1 = gt_array_size(bsp->i1_true_acc);
  len2 = gt_array_size(bsp->i2_true_acc);
  sample_bssm_seqs(bsp->i_false_acc, MAX3(len0, len1, len2));
}

void gth_bssm_seq_processor_write(GthBSSMSeqProcessor *bsp)
{
  gt_assert(bsp);

  /* for the output files we use Volker Brendel's phase notation */
  bssm_seqs_write(bsp->i0_true_don_gt, bsp->gt_t1_fp);
  bssm_seqs_write(bsp->i1_true_don_gt, bsp->gt_t2_fp);
  bssm_seqs_write(bsp->i2_true_don_gt, bsp->gt_t0_fp);
  bssm_seqs_write(bsp->e0_false_don_gt, bsp->gt_f1_fp);
  bssm_seqs_write(bsp->e1_false_don_gt, bsp->gt_f2_fp);
  bssm_seqs_write(bsp->e2_false_don_gt, bsp->gt_f0_fp);
  bssm_seqs_write(bsp->i_false_don_gt, bsp->gt_fi_fp);

  if (bsp->gcdonor) {
    bssm_seqs_write(bsp->i0_true_don_gc, bsp->gc_t1_fp);
    bssm_seqs_write(bsp->i1_true_don_gc, bsp->gc_t2_fp);
    bssm_seqs_write(bsp->i2_true_don_gc, bsp->gc_t0_fp);
    bssm_seqs_write(bsp->e0_false_don_gc, bsp->gc_f1_fp);
    bssm_seqs_write(bsp->e1_false_don_gc, bsp->gc_f2_fp);
    bssm_seqs_write(bsp->e2_false_don_gc, bsp->gc_f0_fp);
    bssm_seqs_write(bsp->i_false_don_gc, bsp->gc_fi_fp);
  }

  bssm_seqs_write(bsp->i0_true_acc, bsp->ag_t1_fp);
  bssm_seqs_write(bsp->i1_true_acc, bsp->ag_t2_fp);
  bssm_seqs_write(bsp->i2_true_acc, bsp->ag_t0_fp);
  bssm_seqs_write(bsp->e0_false_acc, bsp->ag_f1_fp);
  bssm_seqs_write(bsp->e1_false_acc, bsp->ag_f2_fp);
  bssm_seqs_write(bsp->e2_false_acc, bsp->ag_f0_fp);
  bssm_seqs_write(bsp->i_false_acc, bsp->ag_fi_fp);
}
