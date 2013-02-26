/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/fasta.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "ltr/ltrfileout_stream.h"
#include "ltr/ltr_visitor.h"

#define GT_FSWIDTH         60UL
#define GT_MAXFILENAMELEN 256
#define GT_MAXFASTAHEADER 256

struct GtLTRFileOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtEncseq *encseq;
  const char *fileprefix;
  GtFile *metadata_file,
            *tabout_file,
            *pbsout_file,
            *pptout_file,
            *ltr5out_file,
            *ltr3out_file,
            *elemout_file;
  GtHashmap *pdomout_files,
            *pdomali_files,
            *pdomaa_files;
  GtLTRVisitor *lv;
  int tests_to_run;
  unsigned int seqnamelen;
  GtLTRElement element;
  bool write_pdom_alignments,
       write_pdom_aaseqs;
};

#define gt_ltr_fileout_stream_cast(GS)\
        gt_node_stream_cast(gt_ltr_fileout_stream_class(), GS)

static int write_pdom(GtLTRFileOutStream *ls, GtArray *pdoms,
                      const char *pdomname, GtEncseq *seq,
                      unsigned long startpos, char *desc, GtError *err)
{
  int had_err = 0;
  GtFile *seqfile = NULL,
            *alifile = NULL,
            *aafile = NULL;
  unsigned long i = 0,
                seq_length = 0;
  GtStr *pdom_seq,
        *pdom_aaseq;
  gt_error_check(err);

  pdom_seq = gt_str_new();
  pdom_aaseq = gt_str_new();

  /* get protein domain output file */
  seqfile = (GtFile*) gt_hashmap_get(ls->pdomout_files, pdomname);
  if (seqfile == NULL)
  {
    /* no file opened for this domain yet, do it */
    char buffer[GT_MAXFILENAMELEN];
    (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1), "%s_pdom_%s.fas",
                    ls->fileprefix, pdomname);
    seqfile = gt_file_xopen(buffer, "w+");
    gt_hashmap_add(ls->pdomout_files, gt_cstr_dup(pdomname), seqfile);
  }

  /* get protein alignment output file */
  if (ls->write_pdom_alignments)
  {
    alifile = (GtFile*) gt_hashmap_get(ls->pdomali_files, pdomname);
    if (alifile == NULL)
    {
      /* no file opened for this domain yet, do it */
      char buffer[GT_MAXFILENAMELEN];
      (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1), "%s_pdom_%s.ali",
                      ls->fileprefix, pdomname);
      alifile = gt_file_xopen(buffer, "w+");
      gt_hashmap_add(ls->pdomali_files, gt_cstr_dup(pdomname), alifile);
    }
  }

  /* get amino acid sequence output file */
  if (ls->write_pdom_aaseqs)
  {
    aafile = (GtFile*) gt_hashmap_get(ls->pdomaa_files, pdomname);
    if (aafile == NULL)
    {
      /* no file opened for this domain yet, do it */
      char buffer[GT_MAXFILENAMELEN];
      (void) snprintf(buffer, (size_t) (GT_MAXFILENAMELEN-1),
                      "%s_pdom_%s_aa.fas",
                      ls->fileprefix, pdomname);
      aafile = gt_file_xopen(buffer, "w+");
      gt_hashmap_add(ls->pdomaa_files, gt_cstr_dup(pdomname), aafile);
    }
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
    char *tmpstr = NULL;

    fn = *(GtFeatureNode**) gt_array_get(pdoms, i);

    ali = gt_genome_node_get_user_data((GtGenomeNode*) fn, "pdom_alignment");
    aaseq = gt_genome_node_get_user_data((GtGenomeNode*) fn, "pdom_aaseq");
    pdom_rng = gt_genome_node_get_range((GtGenomeNode*) fn);
    tmpstr =  gt_ltrelement_get_sequence(pdom_rng.start-1,
                                         pdom_rng.end-1,
                                         gt_feature_node_get_strand(fn),
                                         seq, startpos, err);
    if (!tmpstr)
    {
      had_err = -1;
      break;
    }
    gt_str_append_cstr(pdom_seq, tmpstr);
    if (ls->write_pdom_alignments && ali)
    {
      char buf[BUFSIZ];
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
    gt_free(tmpstr);
  }
  if (!had_err)
  {
    gt_fasta_show_entry(desc,
                        gt_str_get(pdom_seq),
                        seq_length,
                        GT_FSWIDTH,
                        seqfile);
    if (ls->write_pdom_aaseqs)
    {
      gt_fasta_show_entry(desc,
                          gt_str_get(pdom_aaseq),
                          gt_str_length(pdom_aaseq),
                          GT_FSWIDTH,
                          aafile);
    }
  }
  gt_str_delete(pdom_seq);
  gt_str_delete(pdom_aaseq);
  return had_err;
}

int gt_ltrfileout_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtLTRFileOutStream *ls;
  GtFeatureNode *fn;
  GtRange lltr_rng, rltr_rng, rng, ppt_rng, pbs_rng;
  char *ppt_seq, *pbs_seq;
  int had_err;
  unsigned long i=0;

  gt_error_check(err);
  ls = gt_ltr_fileout_stream_cast(ns);

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

  if (ls->element.mainnode != NULL)
  {
    unsigned long seqnr, seqid_len, seqstartpos;
    GT_UNUSED unsigned long seqlength;
    char *outseq,
         desc[GT_MAXFASTAHEADER];
    GtRange ltr3_rng, ltr5_rng;
    GT_UNUSED GtRange elemrng;

    /* find sequence in GtEncseq */
    const char *sreg = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*)
                                                        ls->element.mainnode));
    (void) sscanf(sreg,"seq%lu", &seqnr);

    elemrng = gt_genome_node_get_range((GtGenomeNode*) ls->element.mainnode);
    seqstartpos = gt_encseq_seqstartpos(ls->encseq, seqnr);
    seqlength = gt_encseq_seqlength(ls->encseq, seqnr);

    ls->element.seqid = gt_calloc((size_t) ls->seqnamelen+1, sizeof (char));
    (void) snprintf(ls->element.seqid,
                    MIN((size_t) seqid_len, (size_t) ls->seqnamelen)+1,
                    "%s",
                    gt_encseq_description(ls->encseq, &seqid_len, seqnr));
    gt_cstr_rep(ls->element.seqid, ' ', '_');
    if (seqid_len > (unsigned long) ls->seqnamelen)
      ls->element.seqid[ls->seqnamelen] = '\0';

    (void) gt_ltrelement_format_description(&ls->element,
                                            ls->seqnamelen,
                                            desc,
                                            (size_t) (GT_MAXFASTAHEADER-1));

    /* output basic retrotransposon data */
    lltr_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.leftLTR);
    rltr_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.rightLTR);
    rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.mainnode);
    gt_file_xprintf(ls->tabout_file,
                       "%lu\t%lu\t%lu\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
                       rng.start,
                       rng.end,
                       gt_ltrelement_length(&ls->element),
                       ls->element.seqid,
                       lltr_rng.start,
                       lltr_rng.end,
                       gt_ltrelement_leftltrlen(&ls->element),
                       rltr_rng.start,
                       rltr_rng.end,
                       gt_ltrelement_rightltrlen(&ls->element));

    /* output TSDs */
    if (ls->element.leftTSD != NULL)
    {
      GtRange tsd_rng;
      GtStrand tsd_strand;
      char *tsd_seq;

      tsd_strand = gt_feature_node_get_strand(ls->element.leftTSD);
      tsd_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.leftTSD);
      tsd_seq = gt_ltrelement_get_sequence(tsd_rng.start-1, tsd_rng.end-1,
                                           tsd_strand, ls->encseq, seqstartpos,
                                           err);
      gt_file_xprintf(ls->tabout_file,
                         "%lu\t%lu\t%s\t",
                         tsd_rng.start,
                         tsd_rng.end,
                         tsd_seq);
      gt_free((char*) tsd_seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t");
    if (ls->element.rightTSD != NULL)
    {
      GtRange tsd_rng;
      GtStrand tsd_strand;
      char *tsd_seq;

      tsd_strand = gt_feature_node_get_strand(ls->element.rightTSD);
      tsd_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.rightTSD);
      tsd_seq = gt_ltrelement_get_sequence(tsd_rng.start-1, tsd_rng.end-1,
                                           tsd_strand, ls->encseq, seqstartpos,
                                           err);
      gt_file_xprintf(ls->tabout_file,
                         "%lu\t%lu\t%s\t",
                         tsd_rng.start,
                         tsd_rng.end,
                         tsd_seq);
      gt_free((char*) tsd_seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t");

    /* output PPT */
    if (ls->element.ppt != NULL)
    {
      GtStrand ppt_strand;
      ppt_strand = gt_feature_node_get_strand(ls->element.ppt);
      ppt_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.ppt);
      ppt_seq = gt_ltrelement_get_sequence(ppt_rng.start-1, ppt_rng.end-1,
                                           ppt_strand, ls->encseq, seqstartpos,
                                           err);
      gt_fasta_show_entry(desc, ppt_seq, gt_range_length(&ppt_rng),
                          GT_FSWIDTH, ls->pptout_file);
      gt_file_xprintf(ls->tabout_file,
                         "%lu\t%lu\t%s\t%c\t%d\t",
                         ppt_rng.start,
                         ppt_rng.end,
                         ppt_seq,
                         GT_STRAND_CHARS[ppt_strand],
                         (ppt_strand == GT_STRAND_FORWARD ?
                             abs((int) (rltr_rng.start - ppt_rng.end)) :
                             abs((int) (lltr_rng.end - ppt_rng.start))));
      gt_free((char*) ppt_seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t\t\t");

    /* output PBS */
    if (ls->element.pbs != NULL)
    {
      GtStrand pbs_strand;

      pbs_strand = gt_feature_node_get_strand(ls->element.pbs);
      pbs_rng = gt_genome_node_get_range((GtGenomeNode*) ls->element.pbs);
      pbs_seq = gt_ltrelement_get_sequence(pbs_rng.start-1, pbs_rng.end-1,
                                           pbs_strand, ls->encseq, seqstartpos,
                                           err);
      gt_fasta_show_entry(desc, pbs_seq, gt_range_length(&pbs_rng),
                          GT_FSWIDTH, ls->pbsout_file);
      gt_file_xprintf(ls->tabout_file,
                         "%lu\t%lu\t%c\t%s\t%s\t%s\t%s\t%s\t",
                         pbs_rng.start,
                         pbs_rng.end,
                         GT_STRAND_CHARS[pbs_strand],
                      gt_feature_node_get_attribute(ls->element.pbs, "trna"),
                      pbs_seq,
                      gt_feature_node_get_attribute(ls->element.pbs,
                                                    "pbsoffset"),
                      gt_feature_node_get_attribute(ls->element.pbs,
                                                    "trnaoffset"),
                      gt_feature_node_get_attribute(ls->element.pbs, "edist"));
      gt_free((char*) pbs_seq);
    } else gt_file_xprintf(ls->tabout_file, "\t\t\t\t\t\t\t\t");

    /* output protein domains */
    if (ls->element.pdoms != NULL)
    {
      GtStr *pdomorderstr = gt_str_new();
      for (i=0;i<gt_array_size(ls->element.pdomorder);i++)
      {
        const char* key = *(const char**) gt_array_get(ls->element.pdomorder,
                                                       i);
        GtArray *entry = (GtArray*) gt_hashmap_get(ls->element.pdoms, key);
        (void) write_pdom(ls, entry, key, ls->encseq, seqstartpos, desc, err);
      }

      if (GT_STRAND_REVERSE == gt_feature_node_get_strand(ls->element.mainnode))
        gt_array_reverse(ls->element.pdomorder);

      for (i=0;i<gt_array_size(ls->element.pdomorder);i++)
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
        ltr5_rng = rltr_rng;
        ltr3_rng = lltr_rng;
        break;
      case GT_STRAND_FORWARD:
      default:
        ltr5_rng = lltr_rng;
        ltr3_rng = rltr_rng;
        break;
    }
    outseq = gt_ltrelement_get_sequence(ltr5_rng.start-1,
                              ltr5_rng.end-1,
                              gt_feature_node_get_strand(ls->element.mainnode),
                              ls->encseq, seqstartpos, err);
    gt_fasta_show_entry(desc,
                        outseq,
                        gt_range_length(&ltr5_rng),
                        GT_FSWIDTH,
                        ls->ltr5out_file);
    gt_free(outseq);
    outseq = gt_ltrelement_get_sequence(rng.start-1,
                              rng.end-1,
                              gt_feature_node_get_strand(ls->element.mainnode),
                              ls->encseq, seqstartpos, err);
    gt_fasta_show_entry(desc,
                        outseq,
                        gt_range_length(&rng),
                        GT_FSWIDTH,
                        ls->elemout_file);
    gt_free(outseq);

    /* output complete oriented element */
    outseq = gt_ltrelement_get_sequence(ltr3_rng.start-1,
                              ltr3_rng.end-1,
                              gt_feature_node_get_strand(ls->element.mainnode),
                              ls->encseq, seqstartpos, err);
    gt_fasta_show_entry(desc,
                        outseq,
                        gt_range_length(&ltr3_rng),
                        GT_FSWIDTH,
                        ls->ltr3out_file);
    gt_free(outseq);
    gt_file_xprintf(ls->tabout_file, "\n");
  }
  gt_hashmap_delete(ls->element.pdoms);
  gt_array_delete(ls->element.pdomorder);
  gt_free(ls->element.seqid);
  return had_err;
}

static void write_metadata(GtFile *metadata_file,
                           int tests_to_run,
                           GtPPTOptions *ppt_opts,
                           GtPBSOptions *pbs_opts,
#ifdef HAVE_HMMER
                           GtPdomOptions *pdom_opts,
#endif
                           const char *trnafilename,
                           const char *seqfilename,
                           const char *gfffilename)
{
  int buflen = 1024;
  char *buffer;

  buffer = gt_calloc((size_t) (buflen+1), sizeof (char));
  /* get working directory */
  while (getcwd(buffer, (size_t) buflen) == NULL) {
    buflen += 1024;
    buffer = gt_realloc(buffer, (buflen+1) * sizeof (char));
  }
  gt_assert(buffer && strlen(buffer) > 0);

  /* append working dir to relative paths if necessary */
  if (seqfilename[0] != '/')
    gt_file_xprintf(metadata_file,
                       "Sequence file used\t%s/%s\n", buffer, seqfilename);
  else
    gt_file_xprintf(metadata_file,
                       "Sequence file used\t%s\n", seqfilename);
  if (gfffilename[0] != '/')
    gt_file_xprintf(metadata_file,
                       "GFF3 input used\t%s/%s\n", buffer, gfffilename);
  else
    gt_file_xprintf(metadata_file,
                       "GFF3 input used\t%s\n", gfffilename);

  if (tests_to_run & GT_LTRDIGEST_RUN_PPT)
  {
    gt_file_xprintf(metadata_file,
                       "PPT length\t%lu-%lunt\t8-30nt\n",
                       ppt_opts->ppt_len.start,
                       ppt_opts->ppt_len.end);
    gt_file_xprintf(metadata_file,
                       "U-box length\t%lu-%lunt\t3-30nt\n",
                       ppt_opts->ubox_len.start,
                       ppt_opts->ubox_len.end);
    gt_file_xprintf(metadata_file,
                       "PPT search radius\t%u\t30\n", ppt_opts->radius);
  }

  if (tests_to_run & GT_LTRDIGEST_RUN_PBS)
  {
    if (trnafilename[0] != '/')
      gt_file_xprintf(metadata_file,
                         "tRNA library for PBS detection\t%s/%s\n",
                         buffer, trnafilename);
    else
      gt_file_xprintf(metadata_file,
                         "tRNA library for PBS detection\t%s\n",
                         trnafilename);
    gt_file_xprintf(metadata_file,
                       "allowed PBS/tRNA alignment length"
                       " range\t%lu-%lunt\t11-30nt\n",
                       pbs_opts->alilen.start,
                       pbs_opts->alilen.end);
    gt_file_xprintf(metadata_file,
                       "PBS/tRNA maximum unit edit distance\t%u\t1\n",
                       pbs_opts->max_edist);
    gt_file_xprintf(metadata_file,
                       "allowed PBS offset from 5' LTR range"
                       "\t%lu-%lunt\t0-5nt\n",
                       pbs_opts->offsetlen.start,
                       pbs_opts->offsetlen.end);
    gt_file_xprintf(metadata_file,
                       "allowed PBS offset from 3' tRNA end"
                       " range\t%lu-%lunt\t0-5nt\n",
                       pbs_opts->trnaoffsetlen.start,
                       pbs_opts->trnaoffsetlen.end);
    gt_file_xprintf(metadata_file,
                       "PBS search radius\t%d\t30\n", pbs_opts->radius);
  }

#ifdef HAVE_HMMER
  if (tests_to_run & GT_LTRDIGEST_RUN_PDOM)
  {
    unsigned long i;
    gt_file_xprintf(metadata_file,
                       "Protein domain models\t%lu (",
                       gt_str_array_size(pdom_opts->hmm_files));
    for (i=0;i<gt_str_array_size(pdom_opts->hmm_files);i++)
    {
      gt_file_xprintf(metadata_file, "%s",
                         gt_str_array_get(pdom_opts->hmm_files, i));
      if (i != gt_str_array_size(pdom_opts->hmm_files)-1)
        gt_file_xprintf(metadata_file, ", ");
    }
    gt_file_xprintf(metadata_file, ")\n");
    gt_file_xprintf(metadata_file,
                       "pHMM e-value cutoff \t%g\t%g\n",
                       pdom_opts->evalue_cutoff,
                       0.000001);
    gt_file_xprintf(metadata_file,
                       "maximal allowed gap length between fragments to chain"
                       " \t%u\t%u\n",
                       pdom_opts->chain_max_gap_length,
                       50);
  }
#endif
  gt_file_xprintf(metadata_file, "\n");
  gt_free(buffer);
}

void gt_ltrfileout_stream_free(GtNodeStream *ns)
{
  GtLTRFileOutStream *ls = gt_ltr_fileout_stream_cast(ns);
  if (ls->tabout_file != NULL)
    gt_file_delete(ls->tabout_file);
  if (ls->metadata_file != NULL)
    gt_file_delete(ls->metadata_file);
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
  gt_hashmap_delete(ls->pdomout_files);
  gt_hashmap_delete(ls->pdomali_files);
  gt_hashmap_delete(ls->pdomaa_files);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_node_stream_delete(ls->in_stream);
}

const GtNodeStreamClass* gt_ltr_fileout_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRFileOutStream),
                                   gt_ltrfileout_stream_free,
                                   gt_ltrfileout_stream_next );
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

void gt_ltr_fileout_stream_enable_pdom_alignment_output(GtNodeStream *ns)
{
  GtLTRFileOutStream *ls;
  gt_assert(ns);

  ls = gt_ltr_fileout_stream_cast(ns);
  ls->write_pdom_alignments = true;
}

void gt_ltr_fileout_stream_enable_aa_sequence_output(GtNodeStream *ns)
{
  GtLTRFileOutStream *ls;
  gt_assert(ns);

  ls = gt_ltr_fileout_stream_cast(ns);
  ls->write_pdom_aaseqs = true;
}

GtNodeStream* gt_ltr_fileout_stream_new(GtNodeStream *in_stream,
                                     int tests_to_run,
                                     GtEncseq *encseq,
                                     char *file_prefix,
                                     GtPPTOptions *ppt_opts,
                                     GtPBSOptions *pbs_opts,
#ifdef HAVE_HMMER
                                     GtPdomOptions *pdom_opts,
#endif
                                     const char *trnafilename,
                                     const char *seqfilename,
                                     const char *gfffilename,
                                     unsigned int seqnamelen,
                                     GtError* err)
{
  GtNodeStream *ns;
  GtLTRFileOutStream *ls;
  char fn[GT_MAXFILENAMELEN];
  gt_error_check(err);

  gt_assert(file_prefix && in_stream && encseq && ppt_opts && pbs_opts
#ifdef HAVE_HMMER
    && pdom_opts
#endif
    );

  ns = gt_node_stream_create(gt_ltr_fileout_stream_class(), true);
  ls = gt_ltr_fileout_stream_cast(ns);

  /* ref GFF input stream and sequences*/
  ls->in_stream = gt_node_stream_ref(in_stream);
  ls->encseq = encseq;
  ls->tests_to_run = tests_to_run;
  ls->seqnamelen = seqnamelen;
  ls->write_pdom_alignments = false;
  ls->write_pdom_aaseqs = false;

  /* open outfiles */
  ls->fileprefix = file_prefix;
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_tabout.csv", file_prefix);
  ls->tabout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  if (tests_to_run & GT_LTRDIGEST_RUN_PPT)
  {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_ppt.fas", file_prefix);
    ls->pptout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  }
  if (tests_to_run & GT_LTRDIGEST_RUN_PBS)
  {
    (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                    "%s_pbs.fas", file_prefix);
    ls->pbsout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  }
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                   "%s_5ltr.fas", file_prefix);
  ls->ltr5out_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_3ltr.fas", file_prefix);
  ls->ltr3out_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_complete.fas", file_prefix);
  ls->elemout_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);
  (void) snprintf(fn, (size_t) (GT_MAXFILENAMELEN-1),
                  "%s_conditions.csv", file_prefix);
  ls->metadata_file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, fn, "w+", err);

  /* create hashmaps to hold protein domain output files */
  ls->pdomout_files = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                     (GtFree) gt_file_delete);
  ls->pdomali_files = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                     (GtFree) gt_file_delete);
  ls->pdomaa_files  = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                     (GtFree) gt_file_delete);

  /* log run conditions in file */
  write_metadata(ls->metadata_file,
                 tests_to_run,
                 ppt_opts,
                 pbs_opts,
#ifdef HAVE_HMMER
                 pdom_opts,
#endif
                 trnafilename,
                 seqfilename,
                 gfffilename);

  /* print tabular outfile headline */
  gt_file_xprintf(ls->tabout_file,
              "element start\telement end\telement length\tsequence\t"
              "lLTR start\tlLTR end\tlLTR length\t"
              "rLTR start\trLTR end\trLTR length\t"
              "lTSD start\tlTSD end\tlTSD motif\t"
              "rTSD start\trTSD end\trTSD motif\t"
              "PPT start\tPPT end\tPPT motif\tPPT strand\tPPT offset");
  gt_file_xprintf(ls->tabout_file,
              "\tPBS start\tPBS end\tPBS strand\ttRNA\ttRNA motif\tPBS offset\t"
              "tRNA offset\tPBS/tRNA edist");
#ifdef HAVE_HMMER
  gt_file_xprintf(ls->tabout_file, "\tProtein domain hits");
#endif
  gt_file_xprintf(ls->tabout_file, "\n");

  /* create visitor */
  ls->lv = (GtLTRVisitor*) gt_ltr_visitor_new(&ls->element);
  return ns;
}
