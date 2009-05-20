/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/fasta.h"
#include "core/string_distri.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_rep.h"
#include "extended/splice_site_info_visitor.h"
#include "extended/reverse.h"

struct GtSpliceSiteInfoVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *region_mapping;
  GtStringDistri *splicesites,
                 *donorsites,
                 *acceptorsites;
  bool show,
       intron_processed;
};

#define splice_site_info_visitor_cast(GV)\
        gt_node_visitor_cast(gt_splice_site_info_visitor_class(), GV)

static void splice_site_info_visitor_free(GtNodeVisitor *gv)
{
  GtSpliceSiteInfoVisitor *splice_site_info_visitor;
  gt_assert(gv);
  splice_site_info_visitor = splice_site_info_visitor_cast(gv);
  gt_region_mapping_delete(splice_site_info_visitor->region_mapping);
  gt_string_distri_delete(splice_site_info_visitor->splicesites);
  gt_string_distri_delete(splice_site_info_visitor->donorsites);
  gt_string_distri_delete(splice_site_info_visitor->acceptorsites);
}

static int process_intron(GtSpliceSiteInfoVisitor *ssiv, GtGenomeNode *intron,
                          GtError *err)
{
  const char *sequence;
  unsigned long seqlen;
  GtStrand strand;
  GtRange range;
  char site[5];
  GtStr *seqid;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(ssiv && intron);
  ssiv->intron_processed = true;
  range = gt_genome_node_get_range(intron);
  gt_assert(range.start); /* 1-based coordinates */
  if (gt_range_length(&range) >= 4) {
    seqid = gt_genome_node_get_seqid(intron);
    had_err = gt_region_mapping_get_raw_sequence(ssiv->region_mapping,
                                                 &sequence, seqid, err);
    if (!had_err) {
      had_err = gt_region_mapping_get_raw_sequence_length(ssiv->region_mapping,
                                                          &seqlen, seqid, err);
    }
    if (!had_err) {
      gt_assert(range.end <= seqlen);
      strand = gt_feature_node_get_strand((GtFeatureNode*) intron);
      if (strand == GT_STRAND_FORWARD || strand == GT_STRAND_REVERSE) {
        /* fill site */
        site[0] = tolower(sequence[range.start-1]);
        site[1] = tolower(sequence[range.start]);
        site[2] = tolower(sequence[range.end-2]);
        site[3] = tolower(sequence[range.end-1]);
        site[4] = '\0';
        if (strand == GT_STRAND_REVERSE)
          had_err = gt_reverse_complement(site, 4, err);
        if (!had_err) {
          /* add site to distributions */
          gt_string_distri_add(ssiv->splicesites, site);
          gt_string_distri_add(ssiv->acceptorsites, site + 2);
          site[2] = '\0';
          gt_string_distri_add(ssiv->donorsites, site);
          ssiv->show = true;
        }
      }
      else {
        gt_warning("skipping intron with unknown orientation "
                   "(file '%s', line %u)", gt_genome_node_get_filename(intron),
                   gt_genome_node_get_line_number(intron));
      }
    }
  }
  return had_err;
}

static int splice_site_info_visitor_genome_feature(GtNodeVisitor *gv,
                                                   GtFeatureNode *fn,
                                                   GtError *err)
{
  GtSpliceSiteInfoVisitor *ssiv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  int had_err = 0;
  gt_error_check(err);
  ssiv = splice_site_info_visitor_cast(gv);
  gt_assert(ssiv->region_mapping);
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_has_type((GtFeatureNode*) node, gt_ft_intron))
      had_err = process_intron(ssiv, (GtGenomeNode*) node, err);
  }
  gt_feature_node_iterator_delete(fni);
  return had_err;
}

const GtNodeVisitorClass* gt_splice_site_info_visitor_class()
{
  static const GtNodeVisitorClass *gvc = NULL;
  if (!gvc) {
   gvc = gt_node_visitor_class_new(sizeof (GtSpliceSiteInfoVisitor),
                                   splice_site_info_visitor_free,
                                   NULL,
                                   splice_site_info_visitor_genome_feature,
                                   NULL,
                                   NULL);
  }
  return gvc;
}

GtNodeVisitor* gt_splice_site_info_visitor_new(GtRegionMapping *rm)
{
  GtNodeVisitor *gv;
  GtSpliceSiteInfoVisitor *ssiv;
  gt_assert(rm);
  gv = gt_node_visitor_create(gt_splice_site_info_visitor_class());
  ssiv = splice_site_info_visitor_cast(gv);
  ssiv->region_mapping = rm;
  ssiv->splicesites = gt_string_distri_new();
  ssiv->acceptorsites = gt_string_distri_new();
  ssiv->donorsites = gt_string_distri_new();
  return gv;
}

static void showsplicesite(const char *string, unsigned long occurrences,
                           double probability, GT_UNUSED void *unused)
{
  gt_assert(string && strlen(string) == 4);
  gt_xputchar(string[0]);
  gt_xputchar(string[1]);
  gt_xputchar('-');
  gt_xputchar(string[2]);
  gt_xputchar(string[3]);
  printf(": %6.2f%% (n=%lu)\n", probability * 100.0, occurrences);
}

static void showsinglesite(const char *string, unsigned long occurrences,
                           double probability, GT_UNUSED void *unused)
{
  gt_assert(string && strlen(string) == 2);
  printf("%s: %6.2f%% (n=%lu)\n", string, probability * 100.0, occurrences);
}

bool gt_splice_site_info_visitor_show(GtNodeVisitor *gv)
{
  GtSpliceSiteInfoVisitor *ssiv;
  gt_assert(gv);
  ssiv = splice_site_info_visitor_cast(gv);

  if (ssiv->show) {
    /* show splice sites */
    printf("splice site distribution (for introns >= 4bp)\n");
    gt_string_distri_foreach(ssiv->splicesites, showsplicesite, NULL);
    gt_xputchar('\n');

    /* show donor sites */
    printf("donor site distribution (for introns >= 4bp)\n");
    gt_string_distri_foreach(ssiv->donorsites, showsinglesite, NULL);
    gt_xputchar('\n');

    /* show acceptor sites */
    printf("acceptor site distribution (for introns >= 4bp)\n");
    gt_string_distri_foreach(ssiv->acceptorsites, showsinglesite, NULL);
  }
  return ssiv->intron_processed;
}
