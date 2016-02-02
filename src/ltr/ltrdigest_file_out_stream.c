/*
  Copyright (c) 2008-2015 Sascha Steinbiss <sascha@steinbiss.name>
  Copyright (c) 2008-2013 Center for Bioinformatics, University of Hamburg

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

#ifndef S_SPLINT_S
#include <ctype.h>
#include <unistd.h>
#endif
#include <stdio.h>
#include <string.h>
#include "core/class_alloc_lock.h"
#include "core/compat.h"
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/str.h"
#include "core/symbol.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_stream_api.h"
#include "extended/reverse_api.h"
#include "ltr/ltrdigest_file_out_stream.h"

#define GT_FSWIDTH         60UL
#define GT_MAXFILENAMELEN   256
#define GT_MAXFASTAHEADER   256

typedef struct GtLTRElement {
  GtUword leftLTR_3,
                leftLTR_5,
                rightLTR_3,
                rightLTR_5;
  GtFeatureNode *mainnode,
                *leftLTR,
                *rightLTR,
                *leftTSD,
                *rightTSD,
                *ppt,
                *pbs;
  char *seqid;
  GtArray *pdomorder;
  GtHashmap *pdoms;
} GtLTRElement;

typedef struct GtLTRVisitor {
  const GtNodeVisitor parent_instance;
  GtLTRElement *element;
} GtLTRVisitor;

struct GtLTRdigestFileOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtRegionMapping *rmap;
  const char *fileprefix;
  GtFile *tabout_file,
         *pbsout_file,
         *pptout_file,
         *ltr5out_file,
         *ltr3out_file,
         *elemout_file;
  GtLTRVisitor *lv;
  int tests_to_run;
  unsigned int seqnamelen;
  GtLTRElement element;
  bool write_pdom_alignments,
       write_pdom_aaseqs;
};

#define gt_ltrdigest_file_out_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrdigest_file_out_stream_class(), GS)

static inline GtUword gt_ltrelement_length(GtLTRElement *e)
{
  gt_assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->rightLTR_3 - e->leftLTR_5 + 1;
}

static inline GtUword gt_ltrelement_leftltrlen(GtLTRElement *e)
{
  gt_assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->leftLTR_3-e->leftLTR_5 + 1;
}

static inline GtUword gt_ltrelement_rightltrlen(GtLTRElement *e)
{
  gt_assert(e && (e->rightLTR_3 >= e->rightLTR_5));
  return e->rightLTR_3 - e->rightLTR_5 + 1;
}

static int gt_ltrelement_format_description(GtLTRElement *e,
                                            unsigned int seqnamelen,
                                            char *buf, size_t buflen)
{
  int ret;
  char tmpstr[BUFSIZ];
  gt_assert(buf && e);

  (void) snprintf(tmpstr, MIN(BUFSIZ, (size_t) seqnamelen+1), "%s", e->seqid);
  tmpstr[seqnamelen+1] = '\0';
  gt_cstr_rep(tmpstr, ' ', '_');
  ret = snprintf(buf, buflen, "%s_"GT_WU"_"GT_WU"", tmpstr, e->leftLTR_5+1,
                 e->rightLTR_3+1);
  return ret;
}

const GtNodeVisitorClass* gt_ltr_visitor_class(void);

#define gt_ltr_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltr_visitor_class(), GV)

static int gt_ltr_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                       GT_UNUSED GtError *err)
{
  GtLTRVisitor *lv;
  GtRange node_range;
  GtArray *pdomarr = NULL;
  const char *pfamname;
  const char *fnt;
  lv = gt_ltr_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  fnt = gt_feature_node_get_type(fn);

  if (strcmp(fnt, gt_ft_LTR_retrotransposon) == 0)
  {
    lv->element->mainnode = fn;
  } else if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0)
  {
    if (lv->element->leftLTR == NULL)
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) fn);
      lv->element->leftLTR = fn;
      /* compensate for 1-based node coords */
      lv->element->leftLTR_5 = node_range.start - 1;
      lv->element->leftLTR_3 = node_range.end - 1;
    }
    else
    {
      node_range = gt_genome_node_get_range((GtGenomeNode*) fn);
      lv->element->rightLTR = fn;
      /* compensate for 1-based node coords */
      lv->element->rightLTR_5 = node_range.start - 1;
      lv->element->rightLTR_3 = node_range.end - 1;
    }
  } else if (strcmp(fnt, gt_ft_target_site_duplication) == 0)
  {
    if (lv->element->leftTSD == NULL)
    {
      lv->element->leftTSD = fn;
    }
    else
    {
      lv->element->rightTSD = fn;
    }
  } else if (strcmp(fnt, gt_ft_RR_tract) == 0)
  {
    if (lv->element->ppt == NULL)
    {
      lv->element->ppt = fn;
    }
  } else if (strcmp(fnt, gt_ft_primer_binding_site) == 0)
  {
    if (lv->element->pbs == NULL)
    {
      lv->element->pbs = fn;
    }
  } else if (strcmp(fnt, gt_ft_protein_match) == 0)
  {
    char buf[BUFSIZ];
    if (!lv->element->pdoms)
    {
      lv->element->pdoms = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                          (GtFree) gt_array_delete);
    }
    pfamname = gt_feature_node_get_attribute(fn, "name");
    (void) snprintf(buf, BUFSIZ-1, "%s", pfamname);
    gt_cstr_rep(buf, '/', '_');
    if (!(pdomarr = (GtArray*) gt_hashmap_get(lv->element->pdoms, buf)))
    {
      char *pfamcpy = gt_cstr_dup(buf);
      pdomarr = gt_array_new(sizeof (GtFeatureNode*));
      gt_hashmap_add(lv->element->pdoms, pfamcpy, pdomarr);
      if (lv->element->pdomorder != NULL)
        gt_array_add(lv->element->pdomorder, pfamcpy);
    }
    gt_array_add(pdomarr, fn);
  }
  return 0;
}

const GtNodeVisitorClass* gt_ltr_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltr_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

static GtNodeVisitor* gt_ltr_visitor_new(GtLTRElement *element)
{
  GtNodeVisitor *nv;
  GtLTRVisitor *lv;
  gt_assert(element);
  nv = gt_node_visitor_create(gt_ltr_visitor_class());
  lv = gt_ltr_visitor_cast(nv);
  lv->element = element;
  gt_assert(lv);
  return nv;
}

static int write_pdom(GtLTRdigestFileOutStream *ls, GtArray *pdoms,
                      const char *pdomname, GT_UNUSED GtRegionMapping *rmap,
                      char *desc, GtError *err)
{
  int had_err = 0;
  GtFile *seqfile = NULL,
         *alifile = NULL,
         *aafile = NULL;
  GtUword i = 0,
                seq_length = 0;
  GtStr *pdom_seq,
        *pdom_aaseq;
  gt_error_check(err);
  char buffer[GT_MAXFILENAMELEN];

  pdom_seq = gt_str_new();
  pdom_aaseq = gt_str_new();

  (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1), "%s_pdom_%s.fas",
                  ls->fileprefix, pdomname);
  seqfile = gt_file_xopen(buffer, "a+");

  /* get protein alignment output file */
  if (ls->write_pdom_alignments)
  {
    (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1), "%s_pdom_%s.ali",
                    ls->fileprefix, pdomname);
    alifile = gt_file_xopen(buffer, "a+");
  }

  /* get amino acid sequence output file */
  if (ls->write_pdom_aaseqs)
  {
    (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_pdom_%s_aa.fas",
                    ls->fileprefix, pdomname);
    aafile = gt_file_xopen(buffer, "a+");
  }

  if (gt_array_size(pdoms) > 1UL)
  {
    for (i=1UL; i<gt_array_size(pdoms); i++)
    {
      gt_assert(gt_genome_node_cmp(*(GtGenomeNode**)gt_array_get(pdoms, i),
                                *(GtGenomeNode**)gt_array_get(pdoms, i-1))
                >= 0);
    }
    if (gt_feature_node_get_strand(*(GtFeatureNode**) gt_array_get(pdoms, 0UL))
        == GT_STRAND_REVERSE)
    {
      gt_array_reverse(pdoms);
    }
  }

  /* output protein domain data */
  for (i=0;i<gt_array_size(pdoms);i++)
  {
    GtRange pdom_rng;
    GtStr *ali,
          *aaseq;
    GtFeatureNode *fn;

    fn = *(GtFeatureNode**) gt_array_get(pdoms, i);

    ali = gt_genome_node_get_user_data((GtGenomeNode*) fn, "pdom_alignment");
    aaseq = gt_genome_node_get_user_data((GtGenomeNode*) fn, "pdom_aaseq");
    pdom_rng = gt_genome_node_get_range((GtGenomeNode*) fn);

    if (gt_extract_feature_sequence(pdom_seq, (GtGenomeNode*) fn,
                                    gt_symbol(gt_ft_protein_match), false,
                                    NULL, NULL, rmap, err))
    {
      had_err = -1;
      break;
    }

    if (ls->write_pdom_alignments && ali)
    {
      char buf[BUFSIZ];
      gt_assert(alifile);
      /* write away alignment */
      (void) snprintf(buf, BUFSIZ-1, "Protein domain alignment in translated "
                                     "sequence for candidate\n'%s':\n\n",
                                     desc);
      gt_file_xwrite(alifile, buf, (size_t) strlen(buf) * sizeof (char));
      gt_file_xwrite(alifile, gt_str_get(ali),
                        (size_t) gt_str_length(ali) * sizeof (char));
      gt_file_xwrite(alifile, "---\n\n", 5 * sizeof (char));
    }
    if (ls->write_pdom_aaseqs && aaseq)
    {
      /* append amino acid sequence */
      gt_str_append_str(pdom_aaseq, aaseq);
    }
    gt_genome_node_release_user_data((GtGenomeNode*) fn, "pdom_alignment");
    gt_genome_node_release_user_data((GtGenomeNode*) fn, "pdom_aaseq");
    seq_length += gt_range_length(&pdom_rng);
  }

  if (!had_err)
  {
    gt_assert(seqfile);
    gt_fasta_show_entry(desc,
                        gt_str_get(pdom_seq),
                        seq_length,
                        GT_FSWIDTH,
                        seqfile);
    if (ls->write_pdom_aaseqs)
    {
      gt_assert(aafile);
      gt_fasta_show_entry(desc,
                          gt_str_get(pdom_aaseq),
                          gt_str_length(pdom_aaseq),
                          GT_FSWIDTH,
                          aafile);
    }
  }
  gt_str_delete(pdom_seq);
  gt_str_delete(pdom_aaseq);
  gt_file_delete(aafile);
  gt_file_delete(alifile);
  gt_file_delete(seqfile);
  return had_err;
}

int gt_ltrfileout_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtLTRdigestFileOutStream *ls;
  GtFeatureNode *fn;
  GtRange lltr_rng = {GT_UNDEF_UWORD, GT_UNDEF_UWORD},
          rltr_rng = {GT_UNDEF_UWORD, GT_UNDEF_UWORD},
          ppt_rng = {GT_UNDEF_UWORD, GT_UNDEF_UWORD},
          pbs_rng = {GT_UNDEF_UWORD, GT_UNDEF_UWORD};
  int had_err;
  GtUword i=0;

  gt_error_check(err);
  ls = gt_ltrdigest_file_out_stream_cast(ns);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (GtLTRElement));

  /* get annotations from parser */
  had_err = gt_node_stream_next(ls->in_stream, gn, err);
  if (!had_err && *gn)
  {
    GtFeatureNodeIterator* gni;
    GtFeatureNode *mygn;

    /* only process feature nodes */
    if (!(fn = gt_feature_node_try_cast(*gn)))
      return 0;

    ls->element.pdomorder = gt_array_new(sizeof (const char*));

    /* fill LTRElement structure from GFF3 subgraph */
    gni = gt_feature_node_iterator_new(fn);
    for (mygn = fn; mygn != NULL; mygn = gt_feature_node_iterator_next(gni))
      (void) gt_genome_node_accept((GtGenomeNode*) mygn,
                                   (GtNodeVisitor*) ls->lv,
                                   err);
    gt_feature_node_iterator_delete(gni);
  }

  if (!had_err && ls->element.mainnode != NULL)
  {
    char desc[GT_MAXFASTAHEADER];
    GtFeatureNode *ltr3, *ltr5;
    GtStr *sdesc, *sreg, *seq;

    /* find sequence in GtEncseq */
    sreg = gt_genome_node_get_seqid((GtGenomeNode*) ls->element.mainnode);

    sdesc = gt_str_new();
    had_err = gt_region_mapping_get_description(ls->rmap, sdesc, sreg, err);

    if (!had_err) {
      GtRange rng;
      ls->element.seqid = gt_calloc((size_t) ls->seqnamelen+1, sizeof (char));
      (void) snprintf(ls->element.seqid,
                      MIN((size_t) gt_str_length(sdesc),
                          (size_t) ls->seqnamelen)+1,
                      "%s", gt_str_get(sdesc));
      gt_cstr_rep(ls->element.seqid, ' ', '_');
      if (gt_str_length(sdesc) > (GtUword) ls->seqnamelen)
        ls->element.seqid[ls->seqnamelen] = '\0';

      (void) gt_ltrelement_format_description(&ls->element,
                                              ls->seqnamelen,
                                              desc,
                                              (size_t) (GT_MAXFASTAHEADER-1));
      gt_str_delete(sdesc);

      /* output basic retrotransposon data */
      lltr_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.leftLTR);
      rltr_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.rightLTR);
      rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.mainnode);
      gt_file_xprintf(ls->tabout_file,
                      GT_WU"\t"GT_WU"\t"GT_WU"\t%s\t"GT_WU"\t"GT_WU"\t"GT_WU"\t"
                      GT_WU"\t"GT_WU"\t"GT_WU"\t",
                      rng.start, rng.end, gt_ltrelement_length(&ls->element),
                      ls->element.seqid, lltr_rng.start, lltr_rng.end,
                      gt_ltrelement_leftltrlen(&ls->element), rltr_rng.start,
                      rltr_rng.end, gt_ltrelement_rightltrlen(&ls->element));
    }
    seq = gt_str_new();

    /* output TSDs */
    if (!had_err && ls->element.leftTSD != NULL)
    {
      GtRange tsd_rng;
      tsd_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.leftTSD);
      had_err = gt_extract_feature_sequence(seq,
                                       (GtGenomeNode*) ls->element.leftTSD,
                                       gt_symbol(gt_ft_target_site_duplication),
                                       false,
                                       NULL, NULL, ls->rmap, err);
      if (!had_err) {
        gt_file_xprintf(ls->tabout_file,
                         ""GT_WU"\t"GT_WU"\t%s\t",
                         tsd_rng.start,
                         tsd_rng.end,
                         gt_str_get(seq));
      }
    gt_str_reset(seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t");

    if (!had_err && ls->element.rightTSD != NULL)
    {
      GtRange tsd_rng;

      tsd_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.rightTSD);
      had_err = gt_extract_feature_sequence(seq,
                                       (GtGenomeNode*) ls->element.rightTSD,
                                       gt_symbol(gt_ft_target_site_duplication),
                                       false,
                                       NULL, NULL, ls->rmap, err);
      if (!had_err) {
        gt_file_xprintf(ls->tabout_file,
                           ""GT_WU"\t"GT_WU"\t%s\t",
                           tsd_rng.start,
                           tsd_rng.end,
                           gt_str_get(seq));
      }
      gt_str_reset(seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t");

    /* output PPT */
    if (!had_err && ls->element.ppt != NULL)
    {
      GtStrand ppt_strand = gt_feature_node_get_strand(ls->element.ppt);

      ppt_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.ppt);
      had_err = gt_extract_feature_sequence(seq,
                                            (GtGenomeNode*) ls->element.ppt,
                                            gt_symbol(gt_ft_RR_tract), false,
                                            NULL, NULL, ls->rmap, err);
      if (!had_err) {
        gt_fasta_show_entry(desc, gt_str_get(seq), gt_range_length(&ppt_rng),
                            GT_FSWIDTH, ls->pptout_file);
        gt_file_xprintf(ls->tabout_file,
                           ""GT_WU"\t"GT_WU"\t%s\t%c\t%d\t",
                           ppt_rng.start,
                           ppt_rng.end,
                           gt_str_get(seq),
                           GT_STRAND_CHARS[ppt_strand],
                           (ppt_strand == GT_STRAND_FORWARD ?
                               abs((int) (rltr_rng.start - ppt_rng.end)) :
                               abs((int) (lltr_rng.end - ppt_rng.start))));
      }
      gt_str_reset(seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t\t\t");

    /* output PBS */
    if (!had_err && ls->element.pbs != NULL)
    {
      GtStrand pbs_strand;

      pbs_strand = gt_feature_node_get_strand(ls->element.pbs);
      pbs_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.pbs);
      had_err = gt_extract_feature_sequence(seq,
                                           (GtGenomeNode*) ls->element.pbs,
                                           gt_symbol(gt_ft_primer_binding_site),
                                           false, NULL, NULL, ls->rmap, err);
      if (!had_err) {
        gt_fasta_show_entry(desc, gt_str_get(seq), gt_range_length(&pbs_rng),
                            GT_FSWIDTH, ls->pbsout_file);
        gt_file_xprintf(ls->tabout_file,
                         ""GT_WU"\t"GT_WU"\t%c\t%s\t%s\t%s\t%s\t%s\t",
                         pbs_rng.start,
                         pbs_rng.end,
                         GT_STRAND_CHARS[pbs_strand],
                         gt_feature_node_get_attribute(ls->element.pbs, "trna"),
                         gt_str_get(seq),
                         gt_feature_node_get_attribute(ls->element.pbs,
                                                       "pbsoffset"),
                         gt_feature_node_get_attribute(ls->element.pbs,
                                                       "trnaoffset"),
                         gt_feature_node_get_attribute(ls->element.pbs,
                                                       "edist"));
      }
      gt_str_reset(seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t\t\t\t\t\t");

    /* output protein domains */
    if (!had_err && ls->element.pdoms != NULL)
    {
      GtStr *pdomorderstr = gt_str_new();
      for (i=0; !had_err && i<gt_array_size(ls->element.pdomorder); i++)
      {
        const char* key = *(const char**) gt_array_get(ls->element.pdomorder,
                                                       i);
        GtArray *entry = (GtArray*) gt_hashmap_get(ls->element.pdoms, key);
        had_err = write_pdom(ls, entry, key, ls->rmap, desc, err);
      }

      if (GT_STRAND_REVERSE == gt_feature_node_get_strand(ls->element.mainnode))
        gt_array_reverse(ls->element.pdomorder);

      for (i=0 ;!had_err && i<gt_array_size(ls->element.pdomorder); i++)
      {
        const char* name = *(const char**) gt_array_get(ls->element.pdomorder,
                                                        i);
        gt_str_append_cstr(pdomorderstr, name);
        if (i != gt_array_size(ls->element.pdomorder)-1)
          gt_str_append_cstr(pdomorderstr, "/");
      }
      gt_file_xprintf(ls->tabout_file, "%s", gt_str_get(pdomorderstr));
      gt_str_delete(pdomorderstr);
    }

    /* output LTRs (we just expect them to exist) */
    switch (gt_feature_node_get_strand(ls->element.mainnode))
    {
      case GT_STRAND_REVERSE:
        ltr5 = ls->element.rightLTR;
        ltr3 = ls->element.leftLTR;
        break;
      case GT_STRAND_FORWARD:
      default:
        ltr5 = ls->element.leftLTR;
        ltr3 = ls->element.rightLTR;
        break;
    }

    if (!had_err) {
      had_err = gt_extract_feature_sequence(seq, (GtGenomeNode*) ltr5,
                                          gt_symbol(gt_ft_long_terminal_repeat),
                                          false,
                                          NULL, NULL, ls->rmap, err);
    }
    if (!had_err) {
      gt_fasta_show_entry(desc, gt_str_get(seq), gt_str_length(seq),
                          GT_FSWIDTH, ls->ltr5out_file);
      gt_str_reset(seq);
    }
    if (!had_err) {
      had_err = gt_extract_feature_sequence(seq, (GtGenomeNode*) ltr3,
                                          gt_symbol(gt_ft_long_terminal_repeat),
                                          false,
                                          NULL, NULL, ls->rmap, err);
    }
    if (!had_err) {
      gt_fasta_show_entry(desc, gt_str_get(seq), gt_str_length(seq),
                          GT_FSWIDTH, ls->ltr3out_file);
      gt_str_reset(seq);
    }

    /* output complete oriented element */
    if (!had_err) {
      had_err = gt_extract_feature_sequence(seq,
                                           (GtGenomeNode*) ls->element.mainnode,
                                           gt_symbol(gt_ft_LTR_retrotransposon),
                                           false,
                                           NULL, NULL, ls->rmap, err);
    }
    if (!had_err) {
      gt_fasta_show_entry(desc,gt_str_get(seq), gt_str_length(seq),
                          GT_FSWIDTH, ls->elemout_file);
      gt_str_reset(seq);
    }
    gt_file_xprintf(ls->tabout_file, "\n");
    gt_str_delete(seq);
  }
  gt_hashmap_delete(ls->element.pdoms);
  gt_array_delete(ls->element.pdomorder);
  gt_free(ls->element.seqid);
  return had_err;
}

void gt_ltrfileout_stream_free(GtNodeStream *ns)
{
  GtLTRdigestFileOutStream *ls = gt_ltrdigest_file_out_stream_cast(ns);
  if (ls->tabout_file != NULL)
    gt_file_delete(ls->tabout_file);
  if (ls->pbsout_file != NULL)
    gt_file_delete(ls->pbsout_file);
  if (ls->pptout_file != NULL)
    gt_file_delete(ls->pptout_file);
  if (ls->ltr5out_file != NULL)
    gt_file_delete(ls->ltr5out_file);
  if (ls->ltr3out_file != NULL)
    gt_file_delete(ls->ltr3out_file);
  if (ls->elemout_file != NULL)
    gt_file_delete(ls->elemout_file);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_node_stream_delete(ls->in_stream);
}

const GtNodeStreamClass* gt_ltrdigest_file_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRdigestFileOutStream),
                                   gt_ltrfileout_stream_free,
                                   gt_ltrfileout_stream_next );
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

void gt_ltrdigest_file_out_stream_enable_pdom_alignment_output(GtNodeStream *ns)
{
  GtLTRdigestFileOutStream *ls;
  gt_assert(ns);
  ls = gt_ltrdigest_file_out_stream_cast(ns);
  ls->write_pdom_alignments = true;
}

void gt_ltrdigest_file_out_stream_enable_aa_sequence_output(GtNodeStream *ns)
{
  GtLTRdigestFileOutStream *ls;
  gt_assert(ns);
  ls = gt_ltrdigest_file_out_stream_cast(ns);
  ls->write_pdom_aaseqs = true;
}

int gt_ltrdigest_file_out_stream_write_metadata(GtLTRdigestFileOutStream *ls,
                                         int tests_to_run,
                                         const char *trnafilename,
                                         const char *gfffilename,
                                         GtRange ppt_len,
                                         GtRange ubox_len,
                                         unsigned int ppt_radius,
                                         GtRange alilen,
                                         unsigned int max_edist,
                                         GtRange offsetlen,
                                         GtRange trnaoffsetlen,
                                         unsigned int pbs_radius,
                                         GtStrArray *hmm_files,
                                         unsigned int chain_max_gap_length,
                                         double evalue_cutoff,
                                         GtError *err)
{
  int buflen = 1024;
  GtFile *metadata_file;
  char *buffer,
       fn[GT_MAXFILENAMELEN];

  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_conditions.csv", ls->fileprefix);
  metadata_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  if (!metadata_file)
    return -1;

  buffer = gt_calloc((size_t) (buflen+1), sizeof (char));
  /* get working directory */
  while (getcwd(buffer, (size_t) buflen) == NULL) {
    buflen += 1024;
    buffer = gt_realloc(buffer, (buflen+1) * sizeof (char));
  }
  gt_assert(buffer && strlen(buffer) > 0);

  /* append working dir to relative paths if necessary */
  if (gfffilename == NULL) {
    gt_file_xprintf(metadata_file,
                       "GFF3 input used\t<stdin>\n");
  } else {
    if (gfffilename[0] != GT_PATH_SEPARATOR)
      gt_file_xprintf(metadata_file,
                         "GFF3 input used\t%s/%s\n", buffer, gfffilename);
    else
      gt_file_xprintf(metadata_file,
                         "GFF3 input used\t%s\n", gfffilename);
  }

  if (tests_to_run & GT_LTRDIGEST_RUN_PPT)
  {
    gt_file_xprintf(metadata_file,
                       "PPT length\t"GT_WU"-"GT_WU"nt\t8-30nt\n",
                       ppt_len.start,
                       ppt_len.end);
    gt_file_xprintf(metadata_file,
                       "U-box length\t"GT_WU"-"GT_WU"nt\t3-30nt\n",
                       ubox_len.start,
                       ubox_len.end);
    gt_file_xprintf(metadata_file,
                       "PPT search radius\t%u\t30\n", ppt_radius);
  }

  if (tests_to_run & GT_LTRDIGEST_RUN_PBS)
  {
    if (trnafilename[0] != GT_PATH_SEPARATOR)
      gt_file_xprintf(metadata_file,
                         "tRNA library for PBS detection\t%s/%s\n",
                         buffer, trnafilename);
    else
      gt_file_xprintf(metadata_file,
                         "tRNA library for PBS detection\t%s\n",
                         trnafilename);
    gt_file_xprintf(metadata_file,
                       "allowed PBS/tRNA alignment length"
                       " range\t"GT_WU"-"GT_WU"nt\t11-30nt\n",
                       alilen.start,
                       alilen.end);
    gt_file_xprintf(metadata_file,
                       "PBS/tRNA maximum unit edit distance\t%u\t1\n",
                       max_edist);
    gt_file_xprintf(metadata_file,
                       "allowed PBS offset from 5' LTR range"
                       "\t"GT_WU"-"GT_WU"nt\t0-5nt\n",
                       offsetlen.start,
                       offsetlen.end);
    gt_file_xprintf(metadata_file,
                       "allowed PBS offset from 3' tRNA end"
                       " range\t"GT_WU"-"GT_WU"nt\t0-5nt\n",
                       trnaoffsetlen.start,
                       trnaoffsetlen.end);
    gt_file_xprintf(metadata_file,
                       "PBS search radius\t%d\t30\n", pbs_radius);
  }

  if (tests_to_run & GT_LTRDIGEST_RUN_PDOM)
  {
    GtUword i;
    gt_file_xprintf(metadata_file,
                       "Protein domain models\t"GT_WU" (",
                       gt_str_array_size(hmm_files));
    for (i=0;i<gt_str_array_size(hmm_files);i++)
    {
      gt_file_xprintf(metadata_file, "%s", gt_str_array_get(hmm_files, i));
      if (i != gt_str_array_size(hmm_files)-1)
        gt_file_xprintf(metadata_file, ", ");
    }
    gt_file_xprintf(metadata_file, ")\n");
    gt_file_xprintf(metadata_file,
                       "pHMM e-value cutoff \t%g\t%g\n",
                       evalue_cutoff, 0.000001);
    gt_file_xprintf(metadata_file,
                       "maximal allowed gap length between fragments to chain"
                       " \t%u\t%u\n",
                       chain_max_gap_length, 50);
  }

  gt_file_xprintf(metadata_file, "\n");
  if (metadata_file != NULL)
    gt_file_delete(metadata_file);
  gt_free(buffer);
  return 0;
}

GtNodeStream* gt_ltrdigest_file_out_stream_new(GtNodeStream *in_stream,
                                        int tests_to_run,
                                        GtRegionMapping *rmap,
                                        char *file_prefix,
                                        unsigned int seqnamelen,
                                        GtError* err)
{
  GtNodeStream *ns;
  GtLTRdigestFileOutStream *ls;
  char fn[GT_MAXFILENAMELEN];
  int had_err = 0;
  gt_error_check(err);

  gt_assert(file_prefix && in_stream && rmap);

  ns = gt_node_stream_create(gt_ltrdigest_file_out_stream_class(), false);
  ls = gt_ltrdigest_file_out_stream_cast(ns);

  /* ref GFF input stream and sequences*/
  ls->in_stream = gt_node_stream_ref(in_stream);
  ls->rmap = rmap;
  ls->tests_to_run = tests_to_run;
  ls->seqnamelen = seqnamelen;
  ls->write_pdom_alignments = false;
  ls->write_pdom_aaseqs = false;
  ls->tabout_file = ls->pptout_file = ls->pbsout_file = ls->ltr3out_file =
    ls->ltr5out_file = ls->elemout_file = NULL;

  /* open outfiles */
  ls->fileprefix = file_prefix;
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_tabout.csv", file_prefix);
  ls->tabout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  if (!ls->tabout_file)
    had_err = -1;
  if (!had_err && (tests_to_run & GT_LTRDIGEST_RUN_PPT))
  {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_ppt.fas", file_prefix);
    ls->pptout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
    if (!ls->pptout_file)
      had_err = -1;
  }
  if (!had_err && (tests_to_run & GT_LTRDIGEST_RUN_PBS))
  {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_pbs.fas", file_prefix);
    ls->pbsout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
    if (!ls->pbsout_file)
      had_err = -1;
  }
  if (!had_err) {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                     "%s_5ltr.fas", file_prefix);
    ls->ltr5out_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
    if (!ls->ltr5out_file)
      had_err = -1;
  }
  if (!had_err) {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_3ltr.fas", file_prefix);
    ls->ltr3out_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
    if (!ls->ltr3out_file)
      had_err = -1;
  }
  if (!had_err) {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_complete.fas", file_prefix);
    ls->elemout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
    if (!ls->elemout_file)
      had_err = -1;
  }

  if (!had_err) {
    /* print tabular outfile headline */
    gt_file_xprintf(ls->tabout_file,
                "element start\telement end\telement length\tsequence\t"
                "lLTR start\tlLTR end\tlLTR length\t"
                "rLTR start\trLTR end\trLTR length\t"
                "lTSD start\tlTSD end\tlTSD motif\t"
                "rTSD start\trTSD end\trTSD motif\t"
                "PPT start\tPPT end\tPPT motif\tPPT strand\tPPT offset");
    gt_file_xprintf(ls->tabout_file,
                "\tPBS start\tPBS end\tPBS strand\ttRNA\ttRNA motif\tPBS "
                "offset\t"
                "tRNA offset\tPBS/tRNA edist");
  #ifdef HAVE_HMMER
    gt_file_xprintf(ls->tabout_file, "\tProtein domain hits");
  #endif
    gt_file_xprintf(ls->tabout_file, "\n");

    /* create visitor */
    ls->lv = (GtLTRVisitor*) gt_ltr_visitor_new(&ls->element);
    gt_assert(ls->lv);
  }

  if (had_err) {
    gt_node_stream_delete(ns);
    ns = NULL;
  }
  return ns;
}
