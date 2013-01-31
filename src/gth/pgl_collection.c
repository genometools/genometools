/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/parseutils.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/consensus_sa.h"
#include "gth/ags_build.h"
#include "gth/pgl_collection.h"
#include "gth/gthcluster.h"
#include "gth/gthsortags.h"
#include "gth/sa_cmp.h"

#define PGLS_DELIMITERCHAR  '-'

struct GthPGLCollection {
  GtArray *pgls;
};

/* necessary to call consensus_sa() */
static GtRange pgl_get_genomic_range(const void *sa)
{
  GthSACluster *sacluster = *(GthSACluster**) sa;
  return gth_sa_range_actual(sacluster->representative);
}

/* necessary to call consensus_sa() */
static GtStrand pgl_get_strand(const void *sa)
{
  GthSACluster *sacluster = *(GthSACluster**) sa;
  gt_assert(sacluster);
  return gth_sa_gen_strand(sacluster->representative);
}

/* necessary to call consensus_sa() */
static void get_exons_func(GtArray *exons, const void *sa)
{
  GthSACluster *sacluster = *(GthSACluster**) sa;
  GtArray *genomicexons;
  GtRange exon;
  unsigned long i;

  genomicexons = gt_array_new(sizeof (GtRange));
  gth_sa_get_exons(sacluster->representative, genomicexons);
  gt_assert(!gt_array_size(exons));
  for (i = 0; i < gt_array_size(genomicexons); i++) {
    exon.start = ((GtRange*) gt_array_get(genomicexons, i))->start;
    exon.end   = ((GtRange*) gt_array_get(genomicexons, i))->end;
    gt_array_add(exons, exon);
  }
  gt_array_delete(genomicexons);
}

/* necessary to call consensus_sa() */
static void process_splice_form_func(GtArray *spliced_alignments_in_form,
                                     GT_UNUSED const void *set_of_sas,
                                     GT_UNUSED unsigned long number_of_sas,
                                     GT_UNUSED size_t size_of_sa,
                                     void *userdata)
{
  GthPGL *pgl = (GthPGL*) userdata;
  GthAGS *ags;
  GtBittab *assemblytab;
  unsigned long i;

  ags = gth_ags_new(pgl);

  assemblytab = gt_bittab_new(gt_array_size(pgl->saclusters));
  for (i = 0; i < gt_array_size(spliced_alignments_in_form); i++) {
    gt_bittab_set_bit(assemblytab, *(unsigned long*)
                                gt_array_get(spliced_alignments_in_form, i));
  }

  gth_build_AGS_from_assembly(ags, assemblytab, pgl->saclusters);

  gt_bittab_delete(assemblytab);

  gt_array_add(pgl->assemblies, ags);
}

void assemble_cluster(GthPGL *pgl, bool disableclustersas)
{
  GthSACluster *sacluster;
  GthSA *sa;
  unsigned long i;

  sacluster = gt_malloc(sizeof (GthSACluster));
  sacluster->representative = *(GthSA**) gt_array_get_first(pgl->alignments);
  sacluster->members = gt_array_new(sizeof (GthSA*));

  for (i = 1; i < gt_array_size(pgl->alignments); i++) {
    sa = *(GthSA**) gt_array_get(pgl->alignments, i);
    if (disableclustersas ||
        gth_sa_cmp_genomic_actual(&sacluster->representative, &sa)) {
      /* spliced alignments differ -> create a new cluster */
      gt_array_add(pgl->saclusters, sacluster);
      sacluster = gt_malloc(sizeof (GthSACluster));
      sacluster->representative = sa;
      sacluster->members = gt_array_new(sizeof (GthSA*));
    }
    else {
      /* spliced alignments are equal -> store new sa also in current cluster */
      gt_array_add(sacluster->members, sa);
    }
  }

  /* store last cluster */
  gt_array_add(pgl->saclusters, sacluster);
}

GthPGLCollection* gth_pgl_collection_new(GthSACollection *sacollection,
                                         bool disableclustersas)
{
  GthPGLCollection *pgl_collection;
  GthPGL *pgl;
  unsigned long i;
  gt_assert(sacollection);

  /* init */
  pgl_collection = gt_malloc(sizeof *pgl_collection);
  pgl_collection->pgls = gt_array_new(sizeof (GthPGL*));

  /* cluster alignments */
  gthclusterSAstoPGLs(pgl_collection->pgls, sacollection);

  /* assemble (clustered) alignments */
  for (i = 0; i < gt_array_size(pgl_collection->pgls); i++) {
    pgl = *(GthPGL**) gt_array_get(pgl_collection->pgls, i);
    /* sort the spliced alignments */
    qsort(gt_array_get_space(pgl->alignments), gt_array_size(pgl->alignments),
          sizeof (GthSA*), gth_sa_cmp_genomic_actual);

    /* cluster spliced alignments which are equal on the genomic sequence.
       this way we only have to consider one spliced alignment for each cluster
       later on */
    assemble_cluster(pgl, disableclustersas);

    /* call consensus phase */
    gt_consensus_sa(gt_array_get_space(pgl->saclusters),
                    gt_array_size(pgl->saclusters),
                    sizeof (GthSACluster*), pgl_get_genomic_range,
                    pgl_get_strand, get_exons_func, process_splice_form_func,
                    pgl);
  }

  return pgl_collection;
}

void gth_pgl_collection_delete(GthPGLCollection *pgl_collection)
{
  unsigned long i;
  if (!pgl_collection) return;
  for (i = 0; i < gt_array_size(pgl_collection->pgls); i++)
    gth_pgl_delete(*(GthPGL**) gt_array_get(pgl_collection->pgls, i));
  gt_array_delete(pgl_collection->pgls);
  gt_free(pgl_collection);
}

void gth_pgl_collection_sortAGSs(GthPGLCollection *pgl_collection,
                                 double sortagswf)
{
  gt_assert(pgl_collection && pgl_collection->pgls);
  /* sort AGSs according to the weighted mean of the average exon score and
     the average splice site probability */
  gth_sortAGSs(pgl_collection->pgls, sortagswf);
}

void gth_pgl_collection_set_max_ags(GthPGLCollection *pgl_collection,
                                    unsigned int maxagsnum)
{
  unsigned long i;
  gt_assert(pgl_collection && maxagsnum && maxagsnum != GT_UNDEF_UINT);
  for (i = 0; i < gt_array_size(pgl_collection->pgls); i++) {
    gth_pgl_set_max_ags(*(GthPGL**) gt_array_get(pgl_collection->pgls, i),
                        maxagsnum);
  }
}

static GthPGL* gth_pgl_collection_get(const GthPGLCollection *pgl_collection,
                                      unsigned long i)
{
  gt_assert(pgl_collection && pgl_collection->pgls);
  gt_assert(i < gt_array_size(pgl_collection->pgls));
  return *(GthPGL**) gt_array_get(pgl_collection->pgls, i);
}

unsigned long gth_pgl_collection_size(const GthPGLCollection *pgl_collection)
{
  gt_assert(pgl_collection && pgl_collection->pgls);
  return gt_array_size(pgl_collection->pgls);
}

void gth_pgl_collection_traverse(const GthPGLCollection *pgl_collection,
                                 GthPGLVisitor *pgl_visitor, GthInput *input,
                                 bool use_desc_ranges)
{
  unsigned long i;
  gt_assert(pgl_collection && pgl_visitor && input);
  gth_pgl_visitor_preface(pgl_visitor, gth_pgl_collection_size(pgl_collection));
  for (i = 0; i < gth_pgl_collection_size(pgl_collection); i++) {
    GtRegionMapping *region_mapping;
    unsigned long file_num, seq_num, offset = 1;
    GthSeqCon *seq_con;
    GthPGL *pgl = gth_pgl_collection_get(pgl_collection, i);
    file_num = gth_pgl_filenum(pgl);
    gth_input_load_genomic_file(input, file_num, false);
    seq_con = gth_input_current_gen_seq_con(input);
    seq_num = gth_pgl_seqnum(pgl);
    if (use_desc_ranges) {
      GtRange descrange;
      GtStr *description = gt_str_new();
      gth_input_get_genomic_description(input, description, file_num, seq_num);
      if (!gt_parse_description_range(gt_str_get(description), &descrange))
        offset = descrange.start;
      gt_str_delete(description);
    }
    region_mapping =
      gt_region_mapping_new_rawseq((const char*)
                                   gth_seq_con_get_orig_seq(seq_con, seq_num),
                                   gth_seq_con_get_length(seq_con, seq_num),
                                   offset);
    gth_pgl_visitor_set_region_mapping(pgl_visitor, region_mapping);
    gth_pgl_visitor_visit_pgl(pgl_visitor, pgl, i);
    gt_region_mapping_delete(region_mapping);
  }
  gth_pgl_visitor_trailer(pgl_visitor);
}
