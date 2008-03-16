/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "libgtcore/log.h"
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/reverse.h"
#include "libgtltr/ltrdigest_stream.h"
#include "libgtltr/ltr_visitor.h"
#include "libgtltr/pbs.h"
#include "libgtltr/ppt.h"
#include "libgtltr/pdom.h"

struct LTRdigestStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  Bioseq *bioseq;
  PBSOptions *pbs_opts;
  PPTOptions *ppt_opts;
  PdomOptions *pdom_opts;
  LTRVisitor *lv;
  Str *ltrdigest_tag;
  LTRElement element;
};

#define ltrdigest_stream_cast(GS)\
        genome_stream_cast(ltrdigest_stream_class(), GS)

int pdom_domain_attach_gff3(void *key, void *value, void *data,
                            UNUSED Error *err)
{
  struct plan7_s *model = (struct plan7_s *) key;
  LTRdigestStream  *ls = (LTRdigestStream *) data;
  PdomHit *hit = (PdomHit*) value;
  Range rng;
  Phase frame = phase_get(hit->best_hit->name[0]);
  rng.start = hit->best_hit->sqfrom;
  rng.end = hit->best_hit->sqto;
  pdom_convert_frame_position(&rng, frame);
  switch(strand_get(hit->best_hit->name[1]))
  {
    case STRAND_FORWARD:
      ltrelement_offset2pos_fwd(&ls->element, &rng, 0,
                                OFFSET_BEGIN_LEFT_LTR);
      break;
    case STRAND_REVERSE:
      ltrelement_offset2pos_rev(&ls->element, &rng, 0,
                                OFFSET_BEGIN_LEFT_LTR);
      break;
    default:
      break;
  }
  GenomeNode *gf = genome_feature_new(gft_protein_match,
                                      rng,
                                      strand_get(hit->best_hit->name[1]),
                                      NULL,
                                      UNDEF_ULONG);
  genome_feature_set_source(gf, ls->ltrdigest_tag);
  genome_feature_set_phase(gf, frame);
  genome_feature_add_attribute((GenomeFeature*) gf,"pfamname", model->name);
  genome_feature_add_attribute((GenomeFeature*) gf,"pfamid", model->acc);
  genome_node_set_seqid((GenomeNode*) gf,
                        genome_node_get_idstr(
                          (GenomeNode*) ls->element.mainnode));

  genome_node_is_part_of_genome_node((GenomeNode*) ls->element.mainnode, gf);
  return 0;
}

void pbs_attach_results_to_gff3(PBSResults *results, LTRElement *element,
                                Str *tag, unsigned int radius)
{
  Range pbs_range;
  GenomeNode *gf;
  char buffer[BUFSIZ];
  pbs_range.start = results->best_hit->start;
  pbs_range.end   = results->best_hit->end;
  switch(results->best_hit->strand)
  {
    case STRAND_FORWARD:
      ltrelement_offset2pos_fwd(element, &pbs_range, radius,
                                OFFSET_END_LEFT_LTR);
      break;
    case STRAND_REVERSE:
      ltrelement_offset2pos_rev(element, &pbs_range, radius,
                                OFFSET_END_LEFT_LTR);
      break;
    default:
      break;
  }
  results->best_hit->start = pbs_range.start;
  results->best_hit->end = pbs_range.end;
  gf = genome_feature_new(gft_primer_binding_site,
                          pbs_range,
                          results->best_hit->strand,
                          NULL,
                          UNDEF_ULONG);
  genome_feature_set_source(gf, tag);
  genome_node_set_seqid((GenomeNode*) gf,
                        genome_node_get_idstr(
                          (GenomeNode*) element->mainnode));
  genome_feature_add_attribute((GenomeFeature*) gf,"trna",
                                results->best_hit->trna);
  snprintf(buffer, BUFSIZ-1, "%lu", results->best_hit->tstart);
  genome_feature_add_attribute((GenomeFeature*) gf,"trnaoffset", buffer);
  snprintf(buffer, BUFSIZ-1, "%lu", results->best_hit->offset);
  genome_feature_add_attribute((GenomeFeature*) gf,"pbsoffset", buffer);
  snprintf(buffer, BUFSIZ-1, "%lu", results->best_hit->edist);
  genome_feature_add_attribute((GenomeFeature*) gf,"edist", buffer);
  genome_node_is_part_of_genome_node((GenomeNode*) element->mainnode, gf);
}

void ppt_attach_results_to_gff3(PPTResults *results, LTRElement *element,
                                Str *tag, unsigned int radius)
{
  Range ppt_range;
  GenomeNode *gf;
  ppt_range.start = results->best_hit->start;
  ppt_range.end   = results->best_hit->end;
  switch(results->best_hit->strand)
  {
    case STRAND_FORWARD:
      ltrelement_offset2pos_fwd(element, &ppt_range, radius,
                                OFFSET_BEGIN_RIGHT_LTR);
      break;
    case STRAND_REVERSE:
      ltrelement_offset2pos_rev(element, &ppt_range, radius,
                                OFFSET_BEGIN_RIGHT_LTR);
      break;
    default:
      break;
  }
  results->best_hit->start = ppt_range.start;
  results->best_hit->end = ppt_range.end;
  gf = genome_feature_new(gft_RR_tract,
                                      ppt_range,
                                      results->best_hit->strand,
                                      NULL,
                                      UNDEF_ULONG);
  genome_feature_set_source(gf, tag);
  genome_node_set_seqid((GenomeNode*) gf,
                        genome_node_get_idstr(
                          (GenomeNode*) element->mainnode));
  genome_node_is_part_of_genome_node((GenomeNode*) element->mainnode, gf);
}

void run_ltrdigest(LTRElement *element, Seq *seq, LTRdigestStream *ls,
                   Error *err)
{
  PPTResults ppt_results;
  PBSResults pbs_results;
  PdomResults pdom_results;
  char *rev_seq;
  const char *base_seq = seq_get_orig(seq)+element->leftLTR_5;
  unsigned long seqlen = ltrelement_length(element);

  /* create reverse strand sequence */
  rev_seq = ma_malloc(sizeof (char) * seqlen+1);
  memcpy(rev_seq, base_seq, sizeof (char) * seqlen);
  rev_seq[seqlen] = '\0';
  reverse_complement(rev_seq, seqlen, err);

  /* PPT finding
   * -----------*/
  ppt_find((const char*) base_seq, (const char*) rev_seq,
           element, &ppt_results, ls->ppt_opts);
  if (ppt_results.best_hit)
  {
    log_log("    PPT: \tscore: %f \t(%c)\n",
            ppt_results.best_hit->score,
            STRANDCHARS[ppt_results.best_hit->strand]);
    ppt_attach_results_to_gff3(&ppt_results, element,
                               ls->ltrdigest_tag, ls->pbs_opts->radius);
  }

  /* PBS finding
   * ----------- */
  pbs_find((const char*) base_seq, (const char*) rev_seq,
           element, &pbs_results, ls->pbs_opts, err);
   if (pbs_results.best_hit)
   {
    log_log("    PBS: \tscore: %f \t(%c, %s)\n",
            pbs_results.best_hit->score,
            STRANDCHARS[pbs_results.best_hit->strand],
            pbs_results.best_hit->trna);
    pbs_attach_results_to_gff3(&pbs_results, element,
                               ls->ltrdigest_tag, ls->pbs_opts->radius);
   }

  /* Protein domain finding
   * ----------------------*/
  pdom_results.domains = hashtable_new(HASH_STRING, NULL, pdom_clear_domain_hit);
  pdom_find((const char*) base_seq, (const char*) rev_seq,
            element, &pdom_results, ls->pdom_opts);
  if (!pdom_results.empty)
  {
    hashtable_foreach(pdom_results.domains,
                      (Hashiteratorfunc) pdom_domain_attach_gff3,
                      ls,
                      err);
  }

  ma_free(rev_seq);
  ppt_clear_results(&ppt_results);
  pbs_clear_results(&pbs_results);
  hashtable_delete(pdom_results.domains);
}

int ltrdigest_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                               Error *e)
{
  LTRdigestStream *ls;
  int had_err;

  error_check(e);
  ls = ltrdigest_stream_cast(gs);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (LTRElement));

  /* get annotations from parser */
  had_err = genome_stream_next_tree(ls->in_stream, gn, e);
  if (!had_err && *gn)
  {
    GenomeNodeIterator* gni;
    GenomeNode *mygn;
    /* fill LTRElement structure from GFF3 subgraph */
    gni = genome_node_iterator_new(*gn);
    for(mygn = *gn; mygn; mygn = genome_node_iterator_next(gni))
      genome_node_accept(mygn, (GenomeVisitor*) ls->lv, e);
    genome_node_iterator_delete(gni);
  }

  if(ls->element.mainnode)
  {
    unsigned long seqid;
    log_log("%lu %lu %lu %lu %lu %lu\n", ls->element.leftLTR_5,
                                         ls->element.rightLTR_3,
                                         ls->element.leftLTR_5,
                                         ls->element.leftLTR_3,
                                         ls->element.rightLTR_5,
                                         ls->element.rightLTR_3);

    /* TODO: use MD5 hashes to identify sequence */
    const char *sreg = str_get(genome_node_get_seqid((GenomeNode*) ls->element.mainnode));
    sscanf(sreg,"seq%lu", &seqid);
    Seq *seq = bioseq_get_seq(ls->bioseq, seqid);

    /* run LTRdigest core routine */
    run_ltrdigest(&ls->element, seq, ls, e);
  }
  return had_err;
}

void ltrdigest_stream_free(GenomeStream *gs)
{
  LTRdigestStream *ls = ltrdigest_stream_cast(gs);
  genome_visitor_delete((GenomeVisitor*) ls->lv);
  str_delete(ls->ltrdigest_tag);
  genome_stream_delete(ls->in_stream);
}

const GenomeStreamClass* ltrdigest_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (LTRdigestStream),
                                         ltrdigest_stream_next_tree,
                                         ltrdigest_stream_free };
  return &gsc;
}

GenomeStream* ltrdigest_stream_new(GenomeStream *in_stream,
                                   Bioseq *bioseq,
                                   PBSOptions *pbs_opts,
                                   PPTOptions *ppt_opts,
                                   PdomOptions *pdom_opts)
{
  GenomeStream *gs;
  LTRdigestStream *ls;
  gs = genome_stream_create(ltrdigest_stream_class(), true);
  ls = ltrdigest_stream_cast(gs);
  ls->in_stream = genome_stream_ref(in_stream);
  ls->ppt_opts = ppt_opts;
  ls->pbs_opts = pbs_opts;
  ls->pdom_opts = pdom_opts;
  ls->bioseq = bioseq;
  ls->ltrdigest_tag = str_new_cstr("LTRdigest");
  ls->lv = (LTRVisitor*) ltr_visitor_new(&ls->element);
  return gs;
}
