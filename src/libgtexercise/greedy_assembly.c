/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/fasta_separator.h"
#include "libgtcore/ma.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtexercise/greedy_assembly.h"
#include "libgtext/union_find.h"

struct GreedyAssembly {
  unsigned long first_fragment,   /* the first fragment in assembly */
                num_of_fragments, /* number of fragments in assembly */
                *next_fragment,   /* maps current fragment to next fragment */
                *overlap;         /* maps current fragment to its overlap with
                                     the previous fragment */
};

static void try_to_select_edge(GreedyAssembly *ga, const Overlap *edge,
                               bool *inedges, bool *outedges, UnionFind *uf,
                               unsigned long *selected)
{
  if (!inedges[edge->end] && !outedges[edge->start] &&
      union_find_find(uf, edge->start) != union_find_find(uf, edge->end)) {
    /* valid edge found */
    ga->next_fragment[edge->start] = edge->end;
    ga->overlap[edge->end] = edge->weight;
    (*selected)++;
    inedges[edge->end] = true;
    outedges[edge->start] = true;
    union_find_union(uf, edge->start, edge->end);
  }
}

static void add_zero_weight_edges(GreedyAssembly *ga, bool *inedges,
                                  bool *outedges, UnionFind *uf,
                                  unsigned long *selected)
{
  Overlap zeroedge = { 0 };
  unsigned long i, j;
  for (i = 0; i < ga->num_of_fragments; i++) {
    for (j = i + 1; j < ga->num_of_fragments; j++) {
      /* forward */
      zeroedge.start = i;
      zeroedge.end = j;
      try_to_select_edge(ga, &zeroedge, inedges, outedges, uf, selected);
      if (*selected + 1 == ga->num_of_fragments)
        return;
      /* reverse */
      zeroedge.start = j;
      zeroedge.end = i;
      try_to_select_edge(ga, &zeroedge, inedges, outedges, uf, selected);
      if (*selected + 1 == ga->num_of_fragments)
        return;
    }
  }
}

static void assemble(GreedyAssembly *ga, FragmentOverlaps *sorted_overlaps)
{
  unsigned long i, current_edge, selected = 0;
  bool *inedges, *outedges;
  UnionFind *uf;
  assert(ga && sorted_overlaps);
  /* init */
  inedges = ma_calloc(sizeof *inedges, ga->num_of_fragments);
  outedges = ma_calloc(sizeof *inedges, ga->num_of_fragments);
  uf = union_find_new(ga->num_of_fragments);
  current_edge = fragment_overlaps_size(sorted_overlaps);
  /* process given edges */
  while (selected  + 1 < ga->num_of_fragments) {
    const Overlap *edge;
    if (!current_edge)
      break; /* no more edges to process */
    edge = fragment_overlaps_get(sorted_overlaps, current_edge - 1);
    try_to_select_edge(ga, edge, inedges, outedges, uf, &selected);
    current_edge--;
  }
  /* add missing edges with weight 0, if necessary */
  if (selected + 1 < ga->num_of_fragments)
    add_zero_weight_edges(ga, inedges, outedges, uf, &selected);
  assert(selected + 1 == ga->num_of_fragments);
  /* determine first fragment */
  for (i = 0; i < ga->num_of_fragments; i++) {
    if (!inedges[i]) {
      ga->first_fragment = i;
      break;
    }
  }
  assert(ga->first_fragment != UNDEF_ULONG);
  assert(selected + 1 == ga->num_of_fragments);
  union_find_delete(uf);
  ma_free(outedges);
  ma_free(inedges);
}

GreedyAssembly* greedy_assembly_new(Bioseq *fragments,
                                    FragmentOverlaps *sorted_overlaps)
{
  GreedyAssembly *ga;
  assert(fragments && sorted_overlaps);
  assert(fragment_overlaps_are_sorted(sorted_overlaps));
  ga = ma_malloc(sizeof *ga);
  ga->first_fragment = UNDEF_ULONG;
  ga->num_of_fragments = bioseq_number_of_sequences(fragments);
  ga->next_fragment = ma_calloc(ga->num_of_fragments, sizeof (unsigned long));
  ga->overlap = ma_calloc(ga->num_of_fragments, sizeof (unsigned long));
  assemble(ga, sorted_overlaps);
  return ga;
}

void greedy_assembly_delete(GreedyAssembly *ga)
{
  if (!ga) return;
  ma_free(ga->overlap);
  ma_free(ga->next_fragment);
  ma_free(ga);
}

typedef void (*ProcFragment)(unsigned long fragnum, unsigned long overlap,
                             void *data);

static void greedy_assembly_show_generic(const GreedyAssembly *ga,
                                         ProcFragment proc_fragment, void *data)
{
  unsigned long i, current_frag;
  assert(ga && proc_fragment);
  current_frag = ga->first_fragment;
  proc_fragment(current_frag, 0, data);
  for (i = 1; i < ga->num_of_fragments; i++) {
    current_frag = ga->next_fragment[current_frag];
    proc_fragment(current_frag, ga->overlap[current_frag], data);
  }
}

static void show_assembly_part(unsigned long fragnum, unsigned long overlap,
                               void *data)
{
  Bioseq *fragments = data;
  unsigned long i, fraglen;
  const char *frag;
  assert(fragments);
  frag = bioseq_get_sequence(fragments, fragnum);
  fraglen = bioseq_get_sequence_length(fragments, fragnum);
  for (i = overlap; i < fraglen; i++)
    xputchar(frag[i]);
}

void greedy_assembly_show(const GreedyAssembly *ga, Bioseq *fragments)
{
  printf("%cAssembled sequence\n", FASTA_SEPARATOR);
  greedy_assembly_show_generic(ga, show_assembly_part, fragments);
  xputchar('\n');
}

static void show_assebly_fragnum(unsigned long fragnum, unsigned long overlap,
                                 UNUSED void *data)
{
  printf("%lu (overlap=%lu)\n", fragnum, overlap);
}

void greedy_assembly_show_path(const GreedyAssembly *ga)
{
  greedy_assembly_show_generic(ga, show_assebly_fragnum, NULL);
}
