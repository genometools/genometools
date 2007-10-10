/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "libgtcore/fasta.h"
#include "libgtcore/stringdistri.h"
#include "libgtcore/translate.h"
#include "libgtcore/warning.h"
#include "libgtcore/xansi.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/splicesiteinfo_visitor.h"
#include "libgtext/reverse.h"

struct SpliceSiteInfoVisitor {
  const GenomeVisitor parent_instance;
  RegionMapping *regionmapping;
  StringDistri *splicesites,
               *donorsites,
               *acceptorsites;
  bool show;
};

#define splicesiteinfo_visitor_cast(GV)\
        genome_visitor_cast(splicesiteinfo_visitor_class(), GV)

static void splicesiteinfo_visitor_free(GenomeVisitor *gv, Env *env)
{
  SpliceSiteInfoVisitor *splicesiteinfo_visitor;
  assert(gv);
  splicesiteinfo_visitor = splicesiteinfo_visitor_cast(gv);
  regionmapping_delete(splicesiteinfo_visitor->regionmapping, env);
  stringdistri_delete(splicesiteinfo_visitor->splicesites, env);
  stringdistri_delete(splicesiteinfo_visitor->donorsites, env);
  stringdistri_delete(splicesiteinfo_visitor->acceptorsites, env);
}

static int process_intron(SpliceSiteInfoVisitor *ssiv, GenomeNode *intron,
                          Env *env)
{
  const char *sequence;
  unsigned long seqlen;
  Strand strand;
  Range range;
  char site[5];
  Str *seqid;
  int had_err = 0;
  env_error_check(env);
  assert(ssiv && intron);
  range = genome_node_get_range(intron);
  assert(range.start); /* 1-based coordinates */
  if (range_length(range) >= 4) {
    seqid = genome_node_get_seqid(intron);
    had_err = regionmapping_get_raw_sequence(ssiv->regionmapping, &sequence,
                                             seqid, env);
    if (!had_err) {
      had_err = regionmapping_get_raw_sequence_length(ssiv->regionmapping,
                                                      &seqlen, seqid, env);
    }
    if (!had_err) {
      assert(range.end <= seqlen);
      strand = genome_feature_get_strand((GenomeFeature*) intron);
      if (strand == STRAND_FORWARD || strand == STRAND_REVERSE) {
        /* fill site */
        site[0] = tolower(sequence[range.start-1]);
        site[1] = tolower(sequence[range.start]);
        site[2] = tolower(sequence[range.end-2]);
        site[3] = tolower(sequence[range.end-1]);
        site[4] = '\0';
        if (strand == STRAND_REVERSE)
          had_err = reverse_complement(site, 4, env);
        if (!had_err) {
          /* add site to distributions */
          stringdistri_add(ssiv->splicesites, site, env);
          stringdistri_add(ssiv->acceptorsites, site + 2, env);
          site[2] = '\0';
          stringdistri_add(ssiv->donorsites, site, env);
          ssiv->show = true;
        }
      }
      else {
        warning("skipping intron with unknown orientation "
                "(file '%s', line %lu)", genome_node_get_filename(intron),
                genome_node_get_line_number(intron));
      }
    }
  }
  return had_err;
}

static int splicesiteinfo_visitor_genome_feature(GenomeVisitor *gv,
                                                 GenomeFeature *gf, Env *env)
{
  SpliceSiteInfoVisitor *ssiv;
  GenomeNodeIterator *gni;
  GenomeNode *node;
  int had_err = 0;
  env_error_check(env);
  ssiv = splicesiteinfo_visitor_cast(gv);
  assert(ssiv->regionmapping);
  gni = genome_node_iterator_new((GenomeNode*) gf, env);
  while (!had_err && (node = genome_node_iterator_next(gni, env))) {
    if (genome_feature_get_type((GenomeFeature*) node) == gft_intron)
      had_err = process_intron(ssiv, node, env);
  }
  genome_node_iterator_delete(gni, env);
  return had_err;
}

const GenomeVisitorClass* splicesiteinfo_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (SpliceSiteInfoVisitor),
                                          splicesiteinfo_visitor_free,
                                          NULL,
                                          splicesiteinfo_visitor_genome_feature,
                                          NULL };
  return &gvc;
}

GenomeVisitor* splicesiteinfo_visitor_new(RegionMapping *rm, Env *env)
{
  GenomeVisitor *gv;
  SpliceSiteInfoVisitor *ssiv;
  assert(rm);
  gv = genome_visitor_create(splicesiteinfo_visitor_class(), env);
  ssiv = splicesiteinfo_visitor_cast(gv);
  ssiv->regionmapping = rm;
  ssiv->splicesites = stringdistri_new(env);
  ssiv->acceptorsites = stringdistri_new(env);
  ssiv->donorsites = stringdistri_new(env);
  return gv;
}

static void showsplicesite(const char *string, unsigned long occurrences,
                           double probability, void *unused)
{
  assert(string && strlen(string) == 4);
  xputchar(string[0]);
  xputchar(string[1]);
  xputchar('-');
  xputchar(string[2]);
  xputchar(string[3]);
  printf(": %6.2f%% (n=%lu)\n", probability * 100.0, occurrences);
}

static void showsinglesite(const char *string, unsigned long occurrences,
                           double probability, void *unused)
{
  assert(string && strlen(string) == 2);
  printf("%s: %6.2f%% (n=%lu)\n", string, probability * 100.0, occurrences);
}

void splicesiteinfo_visitor_show(GenomeVisitor *gv, Env *env)
{
  SpliceSiteInfoVisitor *ssiv;
  env_error_check(env);
  assert(gv);
  ssiv = splicesiteinfo_visitor_cast(gv);

  if (ssiv->show) {
    /* show splice sites */
    printf("splice site distribution (for introns >= 4bp)\n");
    stringdistri_foreach(ssiv->splicesites, showsplicesite, NULL, env);
    xputchar('\n');

    /* show donor sites */
    printf("donor site distribution (for introns >= 4bp)\n");
    stringdistri_foreach(ssiv->donorsites, showsinglesite, NULL, env);
    xputchar('\n');

    /* show acceptor sites */
    printf("acceptor site distribution (for introns >= 4bp)\n");
    stringdistri_foreach(ssiv->acceptorsites, showsinglesite, NULL, env);
  }
}
