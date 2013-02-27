/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/codon_api.h"
#include "core/complement.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma_api.h"
#include "core/strand_api.h"
#include "core/symbol.h"
#include "core/trans_table.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/gff3_defines.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"
#include "extended/reverse_api.h"
#include "extended/snp_annotator_visitor.h"

struct GtSNPAnnotatorVisitor {
  const GtNodeVisitor parent_instance;
  GtFeatureNode *gene;
  GtRegionMapping *rmap;
  GtHashmap *rnaseqs;
  GtTransTable *tt;
  bool own_tt;
  const char *mRNA_type,
             *CDS_type,
             *SNV_type,
             *SNP_type;
};

#define GT_SNP_MISSENSE_EFFECT         "non_conservative_missense_variant"
#define GT_SNP_NONSENSE_EFFECT         "stop_gained"
#define GT_SNP_STOP_LOST_EFFECT        "stop_lost"
#define GT_SNP_SYNONYMOUS_AMINO_EFFECT "synonymous_variant"
#define GT_SNP_SYNONYMOUS_STOP_EFFECT  "stop_retained_variant"

#define snp_annotator_visitor_cast(GV)\
        gt_node_visitor_cast(gt_snp_annotator_visitor_class(), GV)

static void snp_annotator_visitor_free(GtNodeVisitor *nv)
{
  GtSNPAnnotatorVisitor *sav;
  if (!nv) return;
  sav = snp_annotator_visitor_cast(nv);
  gt_genome_node_delete((GtGenomeNode*) sav->gene);
  if (sav->own_tt)
    gt_trans_table_delete(sav->tt);
  gt_region_mapping_delete(sav->rmap);
  gt_hashmap_delete(sav->rnaseqs);
}

static int snp_annotator_classify_snp(GtSNPAnnotatorVisitor *sav,
                                      GtFeatureNode *mRNA,
                                      GtFeatureNode *snp,
                                      unsigned long variant_pos,
                                      unsigned long variant_idx,
                                      char variant_char,
#ifndef NDEBUG
                                      GT_UNUSED char reference_char,
#endif
                                      GT_UNUSED GtError *err)
{
  int had_err = 0;
  char *mrnaseq;
  const char *variant_effect = NULL;
  gt_assert(mRNA && snp && sav);
  gt_log_log("processing variant char %c for SNP %s\n",
               variant_char, gt_feature_node_get_attribute(snp, "Dbxref"));
  mrnaseq = gt_hashmap_get(sav->rnaseqs, mRNA);
  gt_assert(mrnaseq);
  if (mrnaseq) {
    char codon[3],
         variant_codon[3];
    GtStr *effect_string;
    char oldamino,
         newamino;
    GT_UNUSED unsigned long mrnalen;
    unsigned long startpos = variant_pos / GT_CODON_LENGTH,
                  variantoffset = variant_pos % GT_CODON_LENGTH;
    mrnalen = strlen(mrnaseq);
    gt_assert(variant_pos < mrnalen);
    variant_codon[0] = codon[0] = mrnaseq[3*startpos];
    variant_codon[1] = codon[1] = mrnaseq[3*startpos+1];
    variant_codon[2] = codon[2] = mrnaseq[3*startpos+2];
    variant_codon[variantoffset] = variant_char;
#ifndef NDEBUG
    gt_assert(toupper(codon[variantoffset]) == toupper(reference_char));
#endif
    if (gt_trans_table_is_stop_codon(sav->tt, codon[0], codon[1], codon[2])) {
      if (gt_trans_table_is_stop_codon(sav->tt, variant_codon[0],
                                       variant_codon[1], variant_codon[2])) {
        variant_effect = gt_symbol(GT_SNP_SYNONYMOUS_STOP_EFFECT);
      } else {
        variant_effect = gt_symbol(GT_SNP_STOP_LOST_EFFECT);
      }
    } else {
      if (gt_trans_table_is_stop_codon(sav->tt, variant_codon[0],
                                       variant_codon[1], variant_codon[2])) {
        variant_effect = gt_symbol(GT_SNP_NONSENSE_EFFECT);
      } else {
        had_err = gt_trans_table_translate_codon(sav->tt, codon[0], codon[1],
                                                 codon[2], &oldamino, err);
        if (!had_err) {
          had_err = gt_trans_table_translate_codon(sav->tt, variant_codon[0],
                                                   variant_codon[1],
                                                   variant_codon[2],
                                                   &newamino, err);
        }
        if (!had_err) {
          if (newamino == oldamino) {
            variant_effect = gt_symbol(GT_SNP_SYNONYMOUS_AMINO_EFFECT);
          } else {
            variant_effect = gt_symbol(GT_SNP_MISSENSE_EFFECT);
          }
        }
      }
    }
    if (!had_err) {
      const char *var_attrib;
      gt_assert(variant_effect != NULL);
      if ((var_attrib = gt_feature_node_get_attribute(snp,
                                                      GT_GVF_VARIANT_EFFECT))) {
        effect_string = gt_str_new_cstr(var_attrib);
        gt_str_append_cstr(effect_string, ",");
        gt_str_append_cstr(effect_string, variant_effect);
      } else {
        effect_string = gt_str_new_cstr(variant_effect);
      }
      gt_str_append_cstr(effect_string, " ");
      gt_str_append_ulong(effect_string, variant_idx);
      gt_str_append_cstr(effect_string, " ");
      gt_str_append_cstr(effect_string, gt_feature_node_get_type(mRNA));
      gt_str_append_cstr(effect_string, " ");
      gt_str_append_cstr(effect_string,
                         gt_feature_node_get_attribute(mRNA, GT_GFF_ID));
      gt_feature_node_set_attribute(snp, GT_GVF_VARIANT_EFFECT,
                                    gt_str_get(effect_string));
      gt_str_reset(effect_string);
      gt_str_delete(effect_string);
    }
  }

  return had_err;
}

static int snp_annotator_visitor_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn,
                                              GtError *err)
{
  int had_err = 0;
  GtSNPAnnotatorVisitor *sav;
  GtFeatureNodeIterator *fni,
                        *mrnafni;
  GtFeatureNode *curnode,
                *curnode2;
  GtRange snp_rng;
  gt_error_check(err);
  sav = snp_annotator_visitor_cast(nv);

  /* ignore non-nodes */
  if (!fn) return 0;

  /* only process SNPs */
  if (!(gt_feature_node_get_type(fn) == sav->SNV_type ||
        gt_feature_node_get_type(fn) == sav->SNP_type)) {
    return 0;
  }

  fni = gt_feature_node_iterator_new_direct(sav->gene);
  snp_rng = gt_genome_node_get_range((GtGenomeNode*) fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_get_type(curnode) == sav->mRNA_type) {
      GtStrand mrna_strand = gt_feature_node_get_strand(curnode);
#ifndef NDEBUG
      const char *refstr;
#endif
      unsigned long mrnasnppos = 0;
      mrnafni = gt_feature_node_iterator_new(curnode);
      while (!had_err && (curnode2 = gt_feature_node_iterator_next(mrnafni))) {
        if (gt_feature_node_get_type(curnode2) == sav->CDS_type) {
          GtRange cds_rng = gt_genome_node_get_range((GtGenomeNode*) curnode2);
          if (gt_range_overlap(&snp_rng, &cds_rng)) {
            char *mRNA,
                 origchar;
            char *variantchars, *variantptr = NULL;
            GT_UNUSED char *refchars, *refptr = NULL;
            mRNA = (char*) gt_hashmap_get(sav->rnaseqs, curnode);
            gt_assert(mRNA);
            gt_assert(snp_rng.start >= cds_rng.start);
            mrnasnppos += (snp_rng.start - cds_rng.start);
            if (mrna_strand == GT_STRAND_REVERSE)
              mrnasnppos = strlen(mRNA) - mrnasnppos - 1;
            gt_assert(mrnasnppos < strlen(mRNA));
            origchar = mRNA[mrnasnppos];
#ifndef NDEBUG
            refstr = refptr = gt_cstr_dup(gt_feature_node_get_attribute(fn,
                                                         GT_GVF_REFERENCE_SEQ));
            if (!had_err && refstr) {
              if (gt_feature_node_get_strand(curnode) == GT_STRAND_REVERSE) {
                int rval = gt_complement(&origchar, origchar, err);
                gt_assert(rval == 0);
              }
              gt_assert(toupper(origchar) == toupper(refstr[0]));
            }
#endif
            variantchars = variantptr = gt_cstr_dup(
                         gt_feature_node_get_attribute(fn, GT_GVF_VARIANT_SEQ));
            if (!had_err && variantchars) {
              unsigned long i = 0;

              while (!had_err &&
                              (*variantchars != ';' && *variantchars != '\0')) {
                if (*variantchars != ',' && *variantchars != origchar) {
                  char variantchar = *variantchars;
#ifndef NDEBUG
                  char refchar = refstr[0];
                  if (!had_err && mrna_strand == GT_STRAND_REVERSE)
                    had_err = gt_complement(&refchar, refchar, err);
#endif
                  if (!had_err && mrna_strand == GT_STRAND_REVERSE)
                    had_err = gt_complement(&variantchar, variantchar, err);
                  if (!had_err) {
                    had_err = snp_annotator_classify_snp(sav, curnode, fn,
                                                         mrnasnppos,
                                                         i++,
                                                         variantchar,
#ifndef NDEBUG
                                                         refchar,
#endif
                                                         err);
                  }
                } else if (*variantchars == origchar) {
                  i++;
                }
                variantchars++;
              }
              gt_free(variantptr);
              gt_free(refptr);
            }
          } else {
            mrnasnppos += gt_range_length(&cds_rng);
          }
        }
      }
      gt_feature_node_iterator_delete(mrnafni);
    }
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

const GtNodeVisitorClass* gt_snp_annotator_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSNPAnnotatorVisitor),
                                    snp_annotator_visitor_free,
                                    NULL,
                                    snp_annotator_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

static int gt_snp_annotator_visitor_prepare_gene(GtSNPAnnotatorVisitor *sav,
                                                 GtError *err)
{
  GtFeatureNodeIterator *fni,
                        *mrnafni;
  GtFeatureNode *curnode,
                *last_mRNA = NULL;
  GtStr *mrnaseq,
        *seqid;
  int had_err = 0;

  mrnaseq = gt_str_new();
  seqid = gt_genome_node_get_seqid((GtGenomeNode*) sav->gene);
  fni = gt_feature_node_iterator_new(sav->gene);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_get_type(curnode) == sav->mRNA_type) {
      GtFeatureNode *curnode2;
      if (last_mRNA) {
        char *mrna_charseq = gt_calloc(gt_str_length(mrnaseq)+1, sizeof (char));
        (void) strncpy(mrna_charseq, gt_str_get(mrnaseq),
                       gt_str_length(mrnaseq));
        if (gt_feature_node_get_strand(sav->gene) == GT_STRAND_REVERSE) {
          had_err = gt_reverse_complement(mrna_charseq, gt_str_length(mrnaseq),
                                          err);
        }
        if (!had_err) {
          gt_hashmap_add(sav->rnaseqs, last_mRNA, mrna_charseq);
          last_mRNA = curnode;
          gt_str_reset(mrnaseq);
        }
      } else last_mRNA = curnode;
      if (!had_err) {
        mrnafni = gt_feature_node_iterator_new(curnode);
        while (!had_err && (curnode2 =
                                      gt_feature_node_iterator_next(mrnafni))) {
          if (gt_feature_node_get_type(curnode2) == sav->CDS_type) {
            char *tmp;
            GtRange rng = gt_genome_node_get_range((GtGenomeNode*) curnode2);
            had_err = gt_region_mapping_get_sequence(sav->rmap, &tmp, seqid,
                                                     rng.start, rng.end, err);
            if (!had_err) {
              gt_str_append_cstr_nt(mrnaseq, tmp, gt_range_length(&rng));
              gt_free(tmp);
            }
          }
        }
        gt_feature_node_iterator_delete(mrnafni);
      }
    }
  }
  if (!had_err && last_mRNA) {
    char *mrna_charseq = gt_calloc(gt_str_length(mrnaseq)+1, sizeof (char));
    (void) strncpy(mrna_charseq, gt_str_get(mrnaseq), gt_str_length(mrnaseq));
    if (gt_feature_node_get_strand(sav->gene) == GT_STRAND_REVERSE) {
      had_err = gt_reverse_complement(mrna_charseq, gt_str_length(mrnaseq),
                                      err);
    }
    if (!had_err) {
      gt_hashmap_add(sav->rnaseqs, last_mRNA, mrna_charseq);
    }
  }
  gt_feature_node_iterator_delete(fni);
  gt_str_delete(mrnaseq);
  return had_err;
}

GtNodeVisitor* gt_snp_annotator_visitor_new(GtFeatureNode *gene,
                                            GtTransTable *trans_table,
                                            GtRegionMapping *rmap,
                                            GtError *err)
{
  GtNodeVisitor *nv;
  GtSNPAnnotatorVisitor *sav;
  gt_assert(gene && gt_feature_node_get_type(gene) == gt_symbol(gt_ft_gene));
  nv = gt_node_visitor_create(gt_snp_annotator_visitor_class());
  sav = snp_annotator_visitor_cast(nv);
  sav->gene = (GtFeatureNode*) gt_genome_node_ref((GtGenomeNode*) gene);
  sav->rmap = gt_region_mapping_ref(rmap);
  sav->mRNA_type = gt_symbol(gt_ft_mRNA);
  sav->CDS_type = gt_symbol(gt_ft_CDS);
  sav->SNV_type = gt_symbol(gt_ft_SNV);
  sav->SNP_type = gt_symbol(gt_ft_SNP);
  sav->rnaseqs = gt_hashmap_new(GT_HASH_DIRECT, NULL, gt_free_func);
  if (trans_table) {
    sav->tt = trans_table;
    sav->own_tt = false;
  } else {
    sav->tt = gt_trans_table_new_standard(err);
    sav->own_tt = true;
  }
  if (!sav->tt || gt_snp_annotator_visitor_prepare_gene(sav, err) != 0) {
    gt_node_visitor_delete(nv);
    return NULL;
  }
  return nv;
}
