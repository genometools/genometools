/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include <stdint.h>
#include <inttypes.h>
#include "core/array_api.h"
#include "core/arraydef.h"
#include "core/disc_distri_api.h"
#include "core/ensure.h"
#include "core/fasta.h"
#include "core/fileutils.h"
#include "core/format64.h"
#include "core/hashmap-generic.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/progressbar.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/spacecalc.h"
#include "extended/assembly_stats_calculator.h"
#include "match/asqg_writer.h"
#include "match/reads_libraries_table.h"
#include "match/rdj-contigpaths.h"
#include "match/rdj-contig-info.h"
#include "match/rdj-contigs-writer.h"
#include "match/rdj-ensure-output.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-strgraph.h"

/* default representation: */

#ifndef GT_STRGRAPH_SHORT
#ifndef GT_STRGRAPH_BITFIELD
#define GT_STRGRAPH_BITPACK
#endif
#endif

/* Combinations of Edges/Vertices/Counts Representations */

#if defined(GT_STRGRAPH_BITFIELD)
#define GT_STRGRAPH_COUNTS_CHAR
#define GT_STRGRAPH_VERTICES_BITFIELD
#define GT_STRGRAPH_EDGES_BITFIELD
#elif defined(GT_STRGRAPH_SHORT)
#define GT_STRGRAPH_COUNTS_CHAR
#define GT_STRGRAPH_VERTICES_BITFIELD
#define GT_STRGRAPH_EDGES_SHORT
#elif defined(GT_STRGRAPH_BITPACK)
#define GT_STRGRAPH_COUNTS_CHAR
#define GT_STRGRAPH_VERTICES_BITPACK
#define GT_STRGRAPH_EDGES_SINGLE_BITPACK
#endif

/* Type For Vertex numbers */

typedef uint64_t GtStrgraphVnum;
#define FormatGtStrgraphVnum       Formatuint64_t
#define PRINTGtStrgraphVnumcast(X) PRINTuint64_tcast(X)
#define SCANGtStrgraphVnumcast(X)  SCANuint64_tcast(X)
#define GT_STRGRAPH_VNUM_MAX (GtStrgraphVnum)UINT64_MAX

/* Common Helper Macros */

#define GT_STRGRAPH_SERIALIZE_DATA(FP, NOFELEMS, SPACE)\
  gt_file_xwrite((FP), (void*)(SPACE), (sizeof (*(SPACE)) * (NOFELEMS)))

#define GT_STRGRAPH_DESERIALIZE_DATA(FP, NOFELEMS, SPACE)\
  {\
    GT_UNUSED int nofreadbytes;\
    nofreadbytes = gt_file_xread((FP), (void*)(SPACE), (sizeof (*(SPACE))) *\
        (NOFELEMS));\
    gt_assert(nofreadbytes == (int)(sizeof (*(SPACE)) * (NOFELEMS)));\
  }

/* Counts Representation */

/*
 * counts representations define:
 *
 * GtStrgraphCount          a type large at least as the largest count
 * GT_STRGRAPH_COUNT_MAX
 * FormatGtStrgraphCount
 * [PRINT|SCAN]GtStrgraphCountcast
 *
 * GT_STRGRAPH...           type:
 * _COUNTS_REPRESENTATION
 * _DECLARE_COUNTS
 * _INIT_COUNTS             (necessary if counts are not used)
 * _(ALLOC|FREE)_COUNTS
 * _(GET|INC)_COUNT         GtStrgraphCount
 * _[DE]SERIALIZE_COUNTS
 *
 */

#include "match/rdj-strgraph-counts-char-def.h"

/* Vertices Representation */

/* GtStrgraphVmark / GT_STRGRAPH_VMARK_BITS  */
#include "match/rdj-strgraph-vmark-def.h"

/*
 * vertices representations define:
 *
 * GtStrgraphVEdgenum        degree / edge number in a vertex
 * GT_STRGRAPH_V_EDGENUM_MAX
 * FormatGtStrgraphVEdgenum
 * [PRINT|SCAN]GtStrgraphVEdgenumcast
 *
 * GtStrgraphEdgenum         edge number in graph, offset value in edges tables
 * GT_STRGRAPH_EDGENUM_MAX
 * FormatGtStrgraphEdgenum
 * [PRINT|SCAN]GtStrgraphEdgenumcast
 *
 * GT_STRGRAPH_..            type:
 * VERTICES_REPRESENTATION   string literal
 * [SET_]NOFVERTICES         GtStrgraphVnum
 * [ALLOC|FREE]_VERTICES
 * SIZEOF_VERTICES           size_t
 * [DE]SERIALIZE_VERTICES
 *
 * GT_STRGRAPH_V_...
 * [INC_|DEC_]OUTDEG         GtStrgraphVEdgenum
 * [SET_]OFFSET              GtStrgraphEdgenum
 * [SET_]MARK                GtStrgraphVmark
 *
 */

#ifdef GT_STRGRAPH_VERTICES_BITPACK
#include "match/rdj-strgraph-vertices-bitpack-def.h"
#else
#include "match/rdj-strgraph-vertices-bitfield-def.h"
#endif

/* macros depending on vertices definitions,
 *  but common to all representations:
 *
 * GT_STRGRAPH...
 * NOFSPM_MAX                GtStrgraphEdgenum
 * NOFREADS                  GtStrgraphVnum
 *
 * GT_STRGRAPH_V_...
 * NOFEDGES                  GtStrgraphVEdgenum
 * NTH_EDGE_OFFSET           GtStrgraphEdgenum
 * INDEG                     GtStrgraphVEdgenum
 * IS_INTERNAL               bool
 * IS_JUNCTION               bool
 * [B|E]                     GtStrgraphVnum
 * READNUM                   GtUword
 * IS_[E|B]                  bool
 * OTHER                     GtStrgraphVnum
 * CHAR                      char
 */
#include "match/rdj-strgraph-vertices-common-def.h"

/* Edges Representation */

/*
 * edges representations define:
 *
 * GtStrgraphLength           length of an edge
 * FormatGtStrgraphLength
 * GT_STRGRAPH_LENGTH_MAX
 *
 * GT_STRGRAPH_N_READS_MAX     GtUword
 *
 * GT_STRGRAPH_..
 * EDGES_REPRESENTATION        string literal
 * [SET_]NOFEDGES              GtStrgraphEdgenum
 * SIZEOF_EDGES                size_t
 * [ALLOC|SHRINK|FREE]_EDGES
 * COPY_EDGE
 * [DE_]SERIALIZE_EDGES
 * FIND_LONGEST_EDGE           GtStrgraphLength
 * SORT_EDGES_BY_LENGTH_FROM_VERTEX XXX:rename
 *
 * GT_STRGRAPH_EDGE_...
 * [SET_]DEST                  GtStrgraphVnum
 * [SET_]LEN                   GtStrgraphLength
 * INIT
 * [SET|HAS]_MARK              bool
 * REDUCE
 * IS_REDUCED                  bool
 */

#ifdef GT_STRGRAPH_EDGES_SHORT
#include "match/rdj-strgraph-edges-short-def.h"
#elif defined(GT_STRGRAPH_EDGES_BITPACK)
#include "match/rdj-strgraph-edges-bitpack-def.h"
#elif defined(GT_STRGRAPH_EDGES_SINGLE_BITPACK)
#include "match/rdj-strgraph-edges-single-bitpack-def.h"
#else
#include "match/rdj-strgraph-edges-bitfield-def.h"
#endif

/* seqnum where to read label for edge/vertex in mirrored encseq */

#define GT_STRGRAPH_V_MIRROR_SEQNUM(NOFV, V) (GT_STRGRAPH_V_IS_E(V) \
    ? GT_STRGRAPH_V_READNUM(V) \
    : ((GtUword)NOFV - GT_STRGRAPH_V_READNUM(V) - 1))

void gt_strgraph_show_limits_debug_log(void)
{
  gt_log_log("string graph representation: vertices=%s edges=%s",
      GT_STRGRAPH_VERTICES_REPRESENTATION, GT_STRGRAPH_EDGES_REPRESENTATION);
}

void gt_strgraph_show_limits(void)
{
  printf("# max number of reads: "GT_WU"\n", (GtUword) GT_STRGRAPH_N_READS_MAX);
  printf("# max read length: "FormatGtStrgraphLength"\n",
      PRINTGtStrgraphLengthcast(GT_STRGRAPH_LENGTH_MAX));
  printf("# max degree of a vertex: "FormatGtStrgraphVEdgenum"\n",
      PRINTGtStrgraphVEdgenumcast(GT_STRGRAPH_V_EDGENUM_MAX));
  printf("# max number of spm: "FormatGtStrgraphEdgenum"\n",
      PRINTGtStrgraphEdgenumcast(GT_STRGRAPH_NOFSPM_MAX));
  printf("# compile-time chosen representations\n");
  printf("# - counts: %s\n",   GT_STRGRAPH_COUNTS_REPRESENTATION);
  printf("# - vertices: %s\n", GT_STRGRAPH_VERTICES_REPRESENTATION);
  printf("# - edges: %s\n",    GT_STRGRAPH_EDGES_REPRESENTATION);
}

GT_UNUSED
static inline void gt_strgraph_show_current_space(const char *label)
{
  GtUword m, f;
  if (gt_ma_bookkeeping_enabled())
  {
    m = gt_ma_get_space_current();
    f = gt_fa_get_space_current();
    gt_log_log("used space %s: %.2f MB (ma: %.2f MB; fa: %.2f MB)",
        label == NULL ? "" : label, GT_MEGABYTES(m + f), GT_MEGABYTES(m),
        GT_MEGABYTES(f));
  }
}

typedef enum {
  GT_STRGRAPH_PREPARATION,
  GT_STRGRAPH_CONSTRUCTION,
  GT_STRGRAPH_SORTED_BY_L,
  GT_STRGRAPH_LOADED_FROM_FILE,
} GtStrgraphState;

struct GtStrgraph {
  const GtEncseq        *encseq;
  GtStrgraphLength      fixlen;
  GtStrgraphState       state;
  FILE                  *spmfile;
  char                  *spmfile_buffer;
  bool                  load_self_spm;
  bool                  binary_spmlist;
  GtStrgraphLength      minmatchlen;
  GtReadsLibrariesTable *rlt;
  GT_STRGRAPH_DECLARE_COUNTS;
  GT_STRGRAPH_DECLARE_VERTICES;
  GT_STRGRAPH_DECLARE_EDGES;
};

#ifdef GT_STRGRAPH_RUNTIME_CHECKS
#define GT_STRGRAPH_CHECK_NOFREADS(NOFREADS)\
  if ((NOFREADS) > GT_STRGRAPH_N_READS_MAX)\
  {\
    fprintf(stderr, "fatal: overflow\n");\
    fprintf(stderr, "more than "GT_WU" reads ("GT_WU" found)\n",\
        GT_STRGRAPH_N_READS_MAX, (NOFREADS));\
    exit(EXIT_FAILURE);\
  }

#define GT_STRGRAPH_CHECK_LEN(FROM, TO, LEN)\
  if ((LEN) > GT_STRGRAPH_LENGTH_MAX)\
  {\
    fprintf(stderr, "fatal: overflow\n");\
    fprintf(stderr, \
            "edge "GT_WU"%c -> "GT_WU"%c has length > "FormatGtStrgraphLength\
            " ("FormatGtStrgraphLength" found)\n",\
            GT_STRGRAPH_V_READNUM(FROM),\
            GT_STRGRAPH_V_IS_E(FROM) ? 'E' : 'B',\
            GT_STRGRAPH_V_READNUM(TO),\
            GT_STRGRAPH_V_IS_E(TO) ? 'E' : 'B',\
            PRINTGtStrgraphLengthcast(GT_STRGRAPH_LENGTH_MAX),\
            PRINTGtStrgraphLengthcast(LEN));\
    exit(EXIT_FAILURE);\
  }

#define GT_STRGRAPH_CHECK_NOFEDGES(NOFEDGES)\
  if ((NOFEDGES) > GT_STRGRAPH_EDGENUM_MAX)\
  {\
    fprintf(stderr, "fatal: overflow\n");\
    fprintf(stderr, "more than "FormatGtStrgraphEdgenum" spm ("\
        FormatGtStrgraphEdgenum" found)\n",\
        PRINTGtStrgraphEdgenumcast(GT_STRGRAPH_EDGENUM_MAX),\
        PRINTGtStrgraphEdgenumcast(NOFEDGES));\
    exit(EXIT_FAILURE);\
  }

#define GT_STRGRAPH_CHECK_OUTDEG(VNUM, OUTDEG)\
  if ((OUTDEG) > GT_STRGRAPH_V_EDGENUM_MAX)\
  {\
    fprintf(stderr, "fatal: overflow\n");\
    fprintf(stderr, "vertex "GT_WU"%c has more than "\
        FormatGtStrgraphVEdgenum" outgoing edges ("FormatGtStrgraphVEdgenum\
        " found)\n", GT_STRGRAPH_V_READNUM(VNUM),\
        GT_STRGRAPH_V_IS_E(VNUM) ? 'E' : 'B',\
        PRINTGtStrgraphVEdgenumcast(GT_STRGRAPH_V_EDGENUM_MAX),\
        PRINTGtStrgraphVEdgenumcast(OUTDEG));\
    exit(EXIT_FAILURE);\
  }
#else
#define GT_STRGRAPH_CHECK_NOFREADS(NOFREADS)
#define GT_STRGRAPH_CHECK_LEN(FROM, TO, LEN)
#define GT_STRGRAPH_CHECK_NOFEDGES(NOFEDGES)
#define GT_STRGRAPH_CHECK_OUTDEG(VNUM, OUTDEG)
#endif

#define GT_STRGRAPH_SEQLEN(STRGRAPH, READNUM)\
  (((STRGRAPH)->fixlen > 0) ? ((STRGRAPH)->fixlen)\
   : (GtStrgraphLength)gt_encseq_seqlength((STRGRAPH)->encseq,\
     (GtUword)(READNUM)))

GtStrgraphLength gt_strgraph_longest_read(GtStrgraph *strgraph)
{
  if (strgraph->fixlen > 0)
  {
    return strgraph->fixlen;
  }
  else
  {
    GtStrgraphLength seqlen, maxseqlen;
    GtStrgraphVnum i;

    gt_assert(strgraph->encseq != NULL);
    maxseqlen = 0;
    for (i = 0; i < GT_STRGRAPH_NOFREADS(strgraph); i++)
    {
      seqlen = (GtStrgraphLength)gt_encseq_seqlength(strgraph->encseq,
          (GtUword)i);
      if (seqlen > maxseqlen)
        maxseqlen = seqlen;
    }
    gt_assert(maxseqlen > 0);
    gt_assert(sizeof (GtStrgraphLength) >= sizeof (GtUword) ||
        maxseqlen <= GT_STRGRAPH_LENGTH_MAX);
    return maxseqlen;
  }
}

GtStrgraphCount gt_strgraph_largest_count(GtStrgraph *strgraph)
{
  GtStrgraphCount count, maxcount;
  GtStrgraphVnum i;
  maxcount = 0;

  gt_assert(strgraph != NULL &&\
            strgraph->__small_counts != NULL &&
            strgraph->__large_counts != NULL);
  for (i = 0; i < (GT_STRGRAPH_NOFVERTICES(strgraph)); i++)
  {
    GT_STRGRAPH_GET_COUNT(strgraph, count, i);
    if (count > maxcount)
      maxcount = count;
  }
  return maxcount;
}

GtStrgraphEdgenum gt_strgraph_counts_sum(GtStrgraph *strgraph)
{
  GtStrgraphCount count;
  GtStrgraphEdgenum sum;
  GtStrgraphVnum i;
  sum = 0;
  for (i = 0; i < (GT_STRGRAPH_NOFVERTICES(strgraph)); i++)
  {
    gt_assert(sizeof (GtStrgraphEdgenum) >= sizeof (GtStrgraphVEdgenum));
    GT_STRGRAPH_GET_COUNT(strgraph, count, i);
    sum += count;
  }
  return sum;
}

/* --- initialization / preparation --- */

GtStrgraph* gt_strgraph_new(GtUword nofreads)
{
  GtStrgraph *strgraph;
  strgraph = gt_calloc((size_t)1, sizeof (GtStrgraph));
  strgraph->state = GT_STRGRAPH_PREPARATION;
  strgraph->load_self_spm = false;
  strgraph->minmatchlen = GT_STRGRAPH_LENGTH_MAX;
  strgraph->rlt = NULL;
  gt_strgraph_show_limits_debug_log();
  GT_STRGRAPH_CHECK_NOFREADS(nofreads);
  GT_STRGRAPH_SET_NOFVERTICES(strgraph, (GtStrgraphVnum)nofreads << 1);
  GT_STRGRAPH_ALLOC_COUNTS(strgraph, GT_STRGRAPH_NOFVERTICES(strgraph));
  return strgraph;
}

static void gt_strgraph_create_vertices(GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtStrgraphCount c;
  GtStrgraphEdgenum offset;

  gt_assert(strgraph != NULL);
  gt_assert(strgraph->state == GT_STRGRAPH_PREPARATION);

  GT_STRGRAPH_ALLOC_VERTICES(strgraph);

  offset = 0;
  for (i = (GtStrgraphVnum)1; i <= GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    GT_STRGRAPH_GET_COUNT(strgraph, c, i - 1);
    gt_assert(sizeof (GtStrgraphVEdgenum) >= sizeof (GtStrgraphCount));
    GT_STRGRAPH_CHECK_OUTDEG(i - 1, (GtStrgraphVEdgenum)c);
    gt_assert(sizeof (GtStrgraphEdgenum) >= sizeof (GtStrgraphCount));
    offset += (GtStrgraphCount)c;
    GT_STRGRAPH_V_SET_OFFSET(strgraph, i, offset);
  }
  GT_STRGRAPH_CHECK_NOFEDGES(offset);
  GT_STRGRAPH_SET_NOFEDGES(strgraph, offset);
  GT_STRGRAPH_FREE_COUNTS(strgraph);
}

void gt_strgraph_allocate_graph(GtStrgraph *strgraph, GtUword fixlen,
    const GtEncseq *encseq)
{
  gt_assert(strgraph != NULL);
  gt_assert((fixlen == 0 && encseq != NULL)||
            (fixlen > 0 && encseq == NULL));
  gt_assert(sizeof (GtStrgraphLength) >= sizeof (GtUword) ||
      (fixlen > 0 && fixlen <= (GtUword)GT_STRGRAPH_LENGTH_MAX) ||
      (encseq != NULL && gt_encseq_max_seq_length(encseq)));
  strgraph->fixlen = (GtStrgraphLength)fixlen;
  strgraph->encseq = encseq;
  gt_log_log("minmatchlen = "FormatGtStrgraphLength, strgraph->minmatchlen);
  gt_strgraph_create_vertices(strgraph);
  GT_STRGRAPH_ALLOC_EDGES(strgraph);
  strgraph->state = GT_STRGRAPH_CONSTRUCTION;
}

void gt_strgraph_delete(GtStrgraph *strgraph)
{
  if (strgraph != NULL)
  {
    GT_STRGRAPH_FREE_VERTICES(strgraph);
    GT_STRGRAPH_FREE_EDGES(strgraph);
    GT_STRGRAPH_FREE_COUNTS(strgraph);
    gt_reads_libraries_table_delete(strgraph->rlt);
    gt_free(strgraph);
  }
}

int gt_strgraph_load_reads_library_table(GtStrgraph *strgraph,
    FILE *rlt_fp, GtError *err)
{
  gt_assert(strgraph != NULL);
  gt_assert(strgraph->rlt == NULL);
  strgraph->rlt = gt_reads_libraries_table_load(rlt_fp, err);
  return (strgraph->rlt == NULL) ? -1 : 0;
}

static void gt_strgraph_save(const GtStrgraph *strgraph, GtFile *outfp)
{
  gt_assert(strgraph != NULL);
  GT_STRGRAPH_SERIALIZE_VERTICES(strgraph, outfp);
  GT_STRGRAPH_SERIALIZE_EDGES(strgraph, outfp);
}

static void gt_strgraph_load(GtStrgraph *strgraph, GtFile *infp)
{
  gt_assert(strgraph != NULL);
  GT_STRGRAPH_DESERIALIZE_VERTICES(strgraph, infp);
  GT_STRGRAPH_DESERIALIZE_EDGES(strgraph, infp);
}

void gt_strgraph_compact(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphVnum i;
  GtStrgraphEdgenum old_offset, next_free_edge, new_offset;
  GtStrgraphVEdgenum old_nofedges, j;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  next_free_edge = 0;
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    new_offset = next_free_edge;
    old_offset = GT_STRGRAPH_V_OFFSET(strgraph, i);
    old_nofedges = GT_STRGRAPH_V_NOFEDGES(strgraph, i);
    for (j = 0; j < old_nofedges; j++)
    {
      if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
      {
        GT_STRGRAPH_COPY_EDGE(strgraph, old_offset + j, next_free_edge);
        next_free_edge++;
      }
    }
    GT_STRGRAPH_V_SET_OFFSET(strgraph, i, new_offset);
    gt_assert(next_free_edge - GT_STRGRAPH_V_OFFSET(strgraph, i) ==
        (GtStrgraphEdgenum)GT_STRGRAPH_V_OUTDEG(strgraph, i));
    if (show_progressbar)
      progress++;
  }
  GT_STRGRAPH_V_SET_OFFSET(strgraph, GT_STRGRAPH_NOFVERTICES(strgraph),
      next_free_edge);
  if (next_free_edge < GT_STRGRAPH_NOFEDGES(strgraph))
  {
    GT_STRGRAPH_SHRINK_EDGES(strgraph, next_free_edge);
  }

  if (show_progressbar)
    gt_progressbar_stop();
}

static GtFile* gt_strgraph_get_file(const char *indexname, const char *suffix,
    bool write, bool gzipped)
{
  GtFile *file;
  GtStr *filename;
  GtError *err;

  err = gt_error_new();
  filename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(filename, suffix);
  if (!write && !gt_file_exists(gt_str_get(filename)))
  {
    fprintf(stderr, "file %s does not exist\n", gt_str_get(filename));
    exit(EXIT_FAILURE);
  }
  file = gt_file_open(gzipped ? GT_FILE_MODE_GZIP : GT_FILE_MODE_UNCOMPRESSED,
      gt_str_get(filename), write ? "wb" : "r", err);
  if (file == NULL)
  {
    fprintf(stderr, "%s", gt_error_get(err));
    exit(EXIT_FAILURE);
  }
  gt_str_delete(filename);
  gt_error_delete(err);
  return file;
}

GtStrgraph* gt_strgraph_new_from_file(const GtEncseq *encseq,
    GtUword fixlen, const char *indexname, const char *suffix)
{
  GtStrgraph *strgraph;
  GtFile *infp;

  gt_assert(encseq != NULL || fixlen > 0);
  gt_assert(sizeof (GtStrgraphLength) >= sizeof (GtUword) ||
     fixlen <= (GtUword)GT_STRGRAPH_LENGTH_MAX);
  gt_strgraph_show_limits_debug_log();

  strgraph = gt_calloc((size_t)1, sizeof (GtStrgraph));
  strgraph->state = GT_STRGRAPH_LOADED_FROM_FILE;
  strgraph->load_self_spm = false;
  strgraph->minmatchlen = GT_STRGRAPH_LENGTH_MAX;
  strgraph->rlt = NULL;
  strgraph->encseq = encseq;
  strgraph->fixlen = (GtStrgraphLength)fixlen;
  GT_STRGRAPH_INIT_COUNTS(strgraph);
  infp = gt_strgraph_get_file(indexname, suffix, false, false);
  gt_strgraph_load(strgraph, infp);
  gt_file_delete(infp);
  return strgraph;
}

void gt_spmproc_strgraph_count(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *strgraph)
{
  GtStrgraph* g;
  GtStrgraphVnum position;

  gt_assert(strgraph != NULL);
  g = (GtStrgraph*)strgraph;

  gt_assert(g != NULL && g->state == GT_STRGRAPH_PREPARATION);
  position = suffixseq_direct ? GT_STRGRAPH_V_E(suffix_readnum)
                              : GT_STRGRAPH_V_B(suffix_readnum);
  GT_STRGRAPH_INC_COUNT(g, position);

  position = prefixseq_direct ? GT_STRGRAPH_V_B(prefix_readnum)
                              : GT_STRGRAPH_V_E(prefix_readnum);
  GT_STRGRAPH_INC_COUNT(g, position);

  if (g->minmatchlen > (GtStrgraphLength)length)
    g->minmatchlen = (GtStrgraphLength)length;
}

int gt_strgraph_open_spmlist_file(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, bool binary, GtUword bufsize, GtError *err)
{
  strgraph->binary_spmlist = binary;
  strgraph->spmfile = gt_fa_fopen_with_suffix(indexname, suffix,
      binary ? "wb" : "w", err);
  if (strgraph->spmfile == NULL)
    return -1;
  if (bufsize > 0)
  {
    strgraph->spmfile_buffer = NULL;
    (void)setvbuf(strgraph->spmfile, strgraph->spmfile_buffer, (int)_IOFBF,
        (size_t)bufsize << 20);
  }
  if (binary)
  {
    if (GT_STRGRAPH_NOFREADS(strgraph) > UINT32_MAX)
      gt_spmlist_write_header_bin64(strgraph->spmfile);
    else
      gt_spmlist_write_header_bin32(strgraph->spmfile);
  }
  return 0;
}

void gt_spmproc_strgraph_count_and_save(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *strgraph)
{
  GtStrgraph* g;

  gt_assert(strgraph != NULL);
  g = (GtStrgraph*)strgraph;
  gt_spmproc_strgraph_count(suffix_readnum, prefix_readnum, length,
      suffixseq_direct, prefixseq_direct, strgraph);
  if (g->binary_spmlist)
  {
    if (GT_STRGRAPH_NOFREADS(g) > (GtStrgraphVnum)UINT32_MAX)
    {
      /*@ignore@*/
      gt_spmproc_show_bin64(suffix_readnum, prefix_readnum, length,
          suffixseq_direct, prefixseq_direct, g->spmfile);
      /*@end@*/
    }
    else
    {
      /*@ignore@*/
      gt_spmproc_show_bin32(suffix_readnum, prefix_readnum, length,
          suffixseq_direct, prefixseq_direct, g->spmfile);
      /*@end@*/
    }
  }
  else
  {
    GtFile *gtfile = gt_file_new_from_fileptr(g->spmfile);
    gt_spmproc_show_ascii(suffix_readnum, prefix_readnum, length,
        suffixseq_direct, prefixseq_direct, gtfile);
    gt_file_delete_without_handle(gtfile);
  }
}

int gt_strgraph_save_counts(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, GT_UNUSED GtError *err)
{
  GtFile *outfp;

  gt_assert(strgraph != NULL);
  gt_assert(strgraph->state == GT_STRGRAPH_PREPARATION);
  outfp = gt_strgraph_get_file(indexname, suffix, true, false);
  gt_assert(outfp != NULL);
  GT_STRGRAPH_SERIALIZE_COUNTS(strgraph, outfp);
  gt_file_delete(outfp);
  return 0;
}

int gt_strgraph_load_counts(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, GT_UNUSED GtError *err)
{
  GtFile *infp;

  gt_assert(strgraph != NULL);
  gt_assert(strgraph->state == GT_STRGRAPH_PREPARATION);
  infp = gt_strgraph_get_file(indexname, suffix, false, false);
  gt_assert(infp != NULL);
  GT_STRGRAPH_DESERIALIZE_COUNTS(strgraph, infp);
  gt_file_delete(infp);
  return 0;
}

void gt_strgraph_close_spmlist_file(GtStrgraph *strgraph)
{
  gt_fa_fclose(strgraph->spmfile);
  gt_free(strgraph->spmfile_buffer);
}

static void gt_strgraph_mark_empty_edges(GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j, n_empty;

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    gt_assert(GT_STRGRAPH_V_OUTDEG(strgraph, i)
        <= GT_STRGRAPH_V_NOFEDGES(strgraph, i));
    n_empty = GT_STRGRAPH_V_NOFEDGES(strgraph, i) -
      GT_STRGRAPH_V_OUTDEG(strgraph, i);
    if (n_empty > 0)
    {
      for (j = GT_STRGRAPH_V_OUTDEG(strgraph, i);
           j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        GT_STRGRAPH_EDGE_REDUCE(strgraph, i, j);
      }
    }
  }
}

int gt_strgraph_load_spm_from_file(GtStrgraph *strgraph,
    GtUword min_length, bool load_self_spm, GtBitsequence *contained,
    const char *indexname, unsigned int nspmfiles, const char *suffix,
    GtError *err)
{
  int had_err = 0;
  GtStr *filename = gt_str_new();
  GtSpmprocSkipData skipdata;
  unsigned int i;

  gt_assert(strgraph != NULL);
  if (contained != NULL)
  {
    skipdata.out.e.proc = gt_spmproc_strgraph_add;
    skipdata.to_skip = contained;
    skipdata.out.e.data = strgraph;
  }
  strgraph->load_self_spm = load_self_spm;
  for (i = 0; i < nspmfiles; i++)
  {
    gt_str_append_cstr(filename, indexname);
    gt_str_append_char(filename, '.');
    gt_str_append_uint(filename, i);
    gt_str_append_cstr(filename, suffix);
    had_err = gt_spmlist_parse(gt_str_get(filename), min_length,
        contained != NULL ? gt_spmproc_skip : gt_spmproc_strgraph_add,
        contained != NULL ? (void*)&skipdata : (void*)strgraph, err);
    gt_str_reset(filename);
  }
  gt_str_delete(filename);
  if (!had_err)
    gt_strgraph_mark_empty_edges(strgraph);
  return had_err;
}

/* --- construction --- */

void gt_strgraph_set_encseq(GtStrgraph *strgraph, const GtEncseq *encseq)
{
  gt_assert(strgraph != NULL);
  strgraph->encseq = encseq;
}

static inline void gt_strgraph_add_edge(GtStrgraph *strgraph,
    GtStrgraphVnum from, GtStrgraphVnum to, GtStrgraphLength spmlen)
{
  GtStrgraphLength edgelen;
  GtStrgraphVEdgenum next_free_edge;

  gt_assert(strgraph != NULL);

  edgelen = GT_STRGRAPH_SEQLEN(strgraph, GT_STRGRAPH_V_READNUM(to)) - spmlen;

  GT_STRGRAPH_CHECK_LEN(from, to, edgelen);
  next_free_edge = GT_STRGRAPH_V_OUTDEG(strgraph, from);
  GT_STRGRAPH_EDGE_INIT(strgraph, from, next_free_edge);
  GT_STRGRAPH_EDGE_SET_DEST(strgraph, from, next_free_edge, to);
  GT_STRGRAPH_EDGE_SET_LEN(strgraph, from, next_free_edge, edgelen);
  GT_STRGRAPH_V_INC_OUTDEG(strgraph, from);
}

void gt_spmproc_strgraph_add(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *graph)
{
  GtStrgraph *strgraph = graph;
  GtStrgraphLength edgelen = (GtStrgraphLength)length;
  gt_assert(strgraph != NULL);
  gt_assert(suffixseq_direct || prefixseq_direct);

  if (suffix_readnum == prefix_readnum && !strgraph->load_self_spm)
    return;

  if (suffixseq_direct) {
    if (prefixseq_direct) {
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_E(suffix_readnum),
                         GT_STRGRAPH_V_E(prefix_readnum), edgelen);
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_B(prefix_readnum),
                         GT_STRGRAPH_V_B(suffix_readnum), edgelen);
    } else {
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_E(suffix_readnum),
                         GT_STRGRAPH_V_B(prefix_readnum), edgelen);
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_E(prefix_readnum),
                         GT_STRGRAPH_V_B(suffix_readnum), edgelen);
    }
  } else {
    if (prefixseq_direct) {
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_B(suffix_readnum),
                         GT_STRGRAPH_V_E(prefix_readnum), edgelen);
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_B(prefix_readnum),
                         GT_STRGRAPH_V_E(suffix_readnum), edgelen);
    } else {
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_B(suffix_readnum),
                         GT_STRGRAPH_V_B(prefix_readnum), edgelen);
      gt_strgraph_add_edge(strgraph, GT_STRGRAPH_V_E(prefix_readnum),
                         GT_STRGRAPH_V_E(suffix_readnum), edgelen);
    }
  }
}

void gt_strgraph_sort_edges_by_len(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphVnum i;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    GT_STRGRAPH_SORT_V_EDGES(strgraph, i);
    if (show_progressbar)
      progress++;
  }

  strgraph->state = GT_STRGRAPH_SORTED_BY_L;

  if (show_progressbar)
    gt_progressbar_stop();
}

#ifndef NDEBUG
static void gt_strgraph_check_outdegs(GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j, outdeg;
  for (i = 0UL; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    outdeg = 0;
    for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        outdeg++;
    gt_assert(outdeg == GT_STRGRAPH_V_OUTDEG(strgraph, i));
  }
}
#endif

static GtUword gt_strgraph_reduce_marked_edges(GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtUword counter = 0;
  GtStrgraphVEdgenum j;

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) == 0)
      continue;
    for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
    {
      if (GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        continue;
      if (GT_STRGRAPH_EDGE_HAS_MARK(strgraph, i, j))
      {
        GT_STRGRAPH_EDGE_REDUCE(strgraph, i, j);
        GT_STRGRAPH_V_DEC_OUTDEG(strgraph, i);
        counter++;
      }
    }
  }
  return counter;
}

/* return value: number of self edges */
GtUword gt_strgraph_redself(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphVnum vnum;
  GtStrgraphVEdgenum edgenum;
  GtUword counter = 0;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  counter = 0;
  if (show_progressbar)
    gt_progressbar_start(&(progress),
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));
  for (vnum = 0; vnum < GT_STRGRAPH_NOFVERTICES(strgraph); vnum++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, vnum) > 0)
    {
      for (edgenum = 0; edgenum < GT_STRGRAPH_V_NOFEDGES(strgraph, vnum);
          edgenum++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, vnum, edgenum) &&
            GT_STRGRAPH_EDGE_DEST(strgraph, vnum, edgenum) == vnum)
        {
          counter++;
          GT_STRGRAPH_EDGE_REDUCE(strgraph, vnum, edgenum);
          GT_STRGRAPH_V_DEC_OUTDEG(strgraph, vnum);
        }
      }
    }
    if (show_progressbar)
      progress++;
  }
  if (show_progressbar)
    gt_progressbar_stop();

  gt_log_log("self-matches counter: "GT_WU"", counter);
  /* self matches shoud be found twice, check number is even */
  gt_assert((counter & 1) == 0);
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return (counter >> 1);
}

/* return value: number of with-rc edges */
GtUword gt_strgraph_redwithrc(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphVnum vnum;
  GtStrgraphVEdgenum edgenum;
  GtUword counter = 0;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  counter = 0;
  if (show_progressbar)
    gt_progressbar_start(&(progress),
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));
  for (vnum = 0; vnum < GT_STRGRAPH_NOFVERTICES(strgraph); vnum++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, vnum) > 0)
    {
      for (edgenum = 0; edgenum < GT_STRGRAPH_V_NOFEDGES(strgraph, vnum);
          edgenum++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, vnum, edgenum) &&
            GT_STRGRAPH_EDGE_DEST(strgraph, vnum, edgenum) ==
            GT_STRGRAPH_V_OTHER(vnum))
        {
          counter++;
          GT_STRGRAPH_EDGE_REDUCE(strgraph, vnum, edgenum);
          GT_STRGRAPH_V_DEC_OUTDEG(strgraph, vnum);
        }
      }
    }
    if (show_progressbar)
      progress++;
  }
  if (show_progressbar)
    gt_progressbar_stop();

  gt_log_log("withrc-matches counter: "GT_WU"", counter);
  /* withrc matches shoud be found twice, check number is even */
  gt_assert((counter & 1) == 0);
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return (counter >> 1);
}

/* return value: number of transitive edges */
GtUword gt_strgraph_redtrans(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphLength jlen, klen, longest;
  GtStrgraphVEdgenum j, k, l;
  GtStrgraphVnum i, jdest, kdest;
  GtUword counter;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);
  gt_assert(strgraph->state == GT_STRGRAPH_SORTED_BY_L);

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_EDGE_DEST(strgraph, i, j),
            GT_STRGRAPH_V_INPLAY);
      }
      GT_STRGRAPH_FIND_LONGEST_EDGE(strgraph, i, longest);
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        jdest = GT_STRGRAPH_EDGE_DEST(strgraph, i, j);
        jlen = GT_STRGRAPH_EDGE_LEN(strgraph, i, j);
        for (k = 0; k < GT_STRGRAPH_V_NOFEDGES(strgraph, jdest) &&
            GT_STRGRAPH_EDGE_LEN(strgraph, jdest, k) + jlen <= longest; k++)
        {
          kdest = GT_STRGRAPH_EDGE_DEST(strgraph, jdest, k);
          klen = GT_STRGRAPH_EDGE_LEN(strgraph, jdest, k);
          if (GT_STRGRAPH_V_MARK(strgraph, kdest) == GT_STRGRAPH_V_INPLAY)
          {
            for (l = 0; l < GT_STRGRAPH_V_NOFEDGES(strgraph, i); l++)
            {
              if (GT_STRGRAPH_EDGE_DEST(strgraph, i, l) == kdest &&
                  GT_STRGRAPH_EDGE_LEN(strgraph, i, l) == jlen + klen)
              {
                GT_STRGRAPH_EDGE_SET_MARK(strgraph, i, l);
              }
            }
          }
        }
      }
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_EDGE_DEST(strgraph, i, j),
            GT_STRGRAPH_V_VACANT);
      }
    }
    if (show_progressbar)
      progress++;
  }
  if (show_progressbar)
    gt_progressbar_stop();

  counter = gt_strgraph_reduce_marked_edges(strgraph);
  gt_log_log("transitive counter: "GT_WU"", counter);
  /* trans spm should be found twice, check number is even */
  /*gt_assert((counter & 1) == 0);*/
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return (counter >> 1);
}

GtUword gt_strgraph_redsubmax(GtStrgraph *strgraph, bool show_progressbar)
{
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j;
  GtUword counter = 0;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);
  gt_assert(strgraph->state == GT_STRGRAPH_SORTED_BY_L);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
          continue;
        if (GT_STRGRAPH_V_MARK(strgraph, GT_STRGRAPH_EDGE_DEST(strgraph, i, j))
              == GT_STRGRAPH_V_INPLAY)
        {
          counter++;
          GT_STRGRAPH_EDGE_REDUCE(strgraph, i, j);
          GT_STRGRAPH_V_DEC_OUTDEG(strgraph, i);
        }
        GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_INPLAY);
      }
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_EDGE_DEST(strgraph, i, j),
            GT_STRGRAPH_V_VACANT);
      }
    }
    if (show_progressbar)
      progress++;
  }
  if (show_progressbar)
    gt_progressbar_stop();
  gt_log_log("submaximal counter: "GT_WU"", counter);
  /* nonmax spm shoud be found twice, check number is even */
  gt_assert((counter & 1) == 0);
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return (counter >> 1);
}

static inline GtStrgraphVEdgenum gt_strgraph_find_only_edge(
    GtStrgraph *strgraph, GtStrgraphVnum from)
{
  GtStrgraphVEdgenum j;
  gt_assert(GT_STRGRAPH_V_OUTDEG(strgraph, from) == (GtStrgraphVEdgenum)1);
  for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, from); j++)
  {
    if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, from, j))
      return j;
  }
  gt_assert(false); /* outdeg error */
  return 0; /* to avoid warnings */
}

typedef struct {
  GtStrgraphVnum vnum;
  GtStrgraphVEdgenum edgenum;
} GtStrgraphEdgeID;

GtUword gt_strgraph_reddepaths(GtStrgraph *strgraph,
    GtUword maxdepth, bool show_progressbar)
{
  GtStrgraphVnum i, from, to;
  GtStrgraphVEdgenum j, from_to;
  GtUword depth, d;
  GtUword counter = 0, nofdepaths = 0;
  GtUint64 progress = 0;
  bool i_branching;
  GtStrgraphEdgeID *edges;

  gt_assert(strgraph != NULL);

  edges = gt_malloc(sizeof (GtStrgraphEdgeID) * (maxdepth + 1));

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      if (GT_STRGRAPH_V_IS_INTERNAL(strgraph, i))
        continue;
      i_branching =
        (GT_STRGRAPH_V_OUTDEG(strgraph, i) > (GtStrgraphVEdgenum)1 &&
         GT_STRGRAPH_V_INDEG(strgraph, i) > 0) ||
        (GT_STRGRAPH_V_OUTDEG(strgraph, i) == (GtStrgraphVEdgenum)1 &&
         GT_STRGRAPH_V_INDEG(strgraph, i) > (GtStrgraphVEdgenum)1);
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j) &&
            !GT_STRGRAPH_EDGE_HAS_MARK(strgraph, i, j))
        {
          from = i;
          from_to = j;
          to = GT_STRGRAPH_EDGE_DEST(strgraph, from, from_to);
          edges->vnum = from;
          edges->edgenum = from_to;
          depth = 1UL;
          while (GT_STRGRAPH_V_IS_INTERNAL(strgraph, to) &&
              depth <= maxdepth)
          {
            depth++;
            from = to;
            from_to = gt_strgraph_find_only_edge(strgraph, from);
            to = GT_STRGRAPH_EDGE_DEST(strgraph, from, from_to);
            gt_assert(depth >= 1UL);
            gt_assert(depth - 1UL <= maxdepth);
            edges[depth - 1UL].vnum = from;
            edges[depth - 1UL].edgenum = from_to;
          }
          if (depth <= maxdepth &&
              (!i_branching || GT_STRGRAPH_V_OUTDEG(strgraph, to) == 0))
          {
            nofdepaths++;
            for (d = 0; d < depth; d++)
            {
              GT_STRGRAPH_EDGE_SET_MARK(strgraph, edges[d].vnum,
                  edges[d].edgenum);
            }
          }
        }
      }
    }
    if (show_progressbar)
      progress++;
  }
  gt_free(edges);
  counter = gt_strgraph_reduce_marked_edges(strgraph);
  if (show_progressbar)
    gt_progressbar_stop();
  gt_log_log("dead-paths = "GT_WU"", nofdepaths);
  gt_log_log("dead-path edges = "GT_WU"", counter);
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return counter;
}

typedef struct {
  GtStrgraphVEdgenum edgenum;
  GtStrgraphVnum dest;
  GtUword depth;
  GtUword width;
} GtStrgraphPathInfo;

static int gt_strgraph_path_info_compare(const void *pi_a, const void *pi_b)
{
  int retv;
  const GtStrgraphPathInfo *a = pi_a, *b = pi_b;
  retv = (int)(a->dest > b->dest) -
    (int)(a->dest < b->dest);
  if (retv == 0)
    retv = (int)(a->width > b->width) - (int)(a->width < b->width);
  return retv;
}

GtUword gt_strgraph_redpbubbles(GtStrgraph *strgraph,
    GtUword maxwidth, const GtUword maxdiff,
    bool show_progressbar)
{
  GtStrgraphVnum i, from, to;
  GtStrgraphVEdgenum j, from_to, p, nofpaths;
  GtStrgraphLength len;
  GtUword depth, width, counter = 0, nofpbubbles = 0;
  GtUint64 progress = 0;
  GtStrgraphPathInfo *info, *prev;

  gt_assert(strgraph != NULL);

  if (maxwidth == 0)
    maxwidth = (GtUword)(gt_strgraph_longest_read(strgraph) << 2) -
        (strgraph->minmatchlen << 1) - 1;
  gt_log_log("redpbubbles(maxwidth="GT_WU", maxdiff="GT_WU")", maxwidth,
             maxdiff);

  /* allocate info and set all marks to VACANT */
  {
    GtUword maxoutdeg = 0;
    for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    {
      GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);
      if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > (GtStrgraphVEdgenum)maxoutdeg)
        maxoutdeg = (GtUword)GT_STRGRAPH_V_OUTDEG(strgraph, i);
    }
    gt_log_log("maxoutdeg = "GT_WU"", maxoutdeg);
    info = gt_malloc(sizeof (GtStrgraphPathInfo) * maxoutdeg);
  }

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      if (GT_STRGRAPH_V_IS_INTERNAL(strgraph, i))
        continue;
      nofpaths = 0;
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        {
          to = GT_STRGRAPH_EDGE_DEST(strgraph, i, j);
          depth = 1UL;
          len = GT_STRGRAPH_EDGE_LEN(strgraph, i, j);
          gt_assert(sizeof (GtUword) >= sizeof (GtStrgraphLength) ||
              len <= (GtStrgraphLength)ULONG_MAX);
          width = (GtUword)len;
          while (GT_STRGRAPH_V_IS_INTERNAL(strgraph, to) && width <= maxwidth)
          {
            depth++;
            from = to;
            from_to = gt_strgraph_find_only_edge(strgraph, from);
            len = GT_STRGRAPH_EDGE_LEN(strgraph, from, from_to);
            gt_assert((sizeof (GtUword) >= sizeof (GtStrgraphLength) &&
                width <= ULONG_MAX - (GtUword)len) ||
                (GtStrgraphLength)width + len < (GtStrgraphLength)ULONG_MAX);
            width += (GtUword)len;
            to = GT_STRGRAPH_EDGE_DEST(strgraph, from, from_to);
          }
          if (width <= maxwidth && depth > 1UL)
          {
            info[nofpaths].edgenum = j;
            info[nofpaths].dest = to;
            info[nofpaths].depth = depth;
            info[nofpaths].width = width;
            nofpaths++;
          }
        }
      }
      if (nofpaths > 0)
      {
        qsort(info, (size_t)nofpaths, sizeof (*info),
            gt_strgraph_path_info_compare);
        prev = info;
        for (p = (GtStrgraphVEdgenum)1; p < nofpaths; p++)
        {
          if (info[p].dest == prev->dest &&
              (info[p].width - prev->width <= maxdiff))
          {
            nofpbubbles++;
            if (info[p].depth <= prev->depth)
            {
              from_to = info[p].edgenum;
            }
            else
            {
              from_to = prev->edgenum;
              prev = info + p;
            }
            GT_STRGRAPH_EDGE_SET_MARK(strgraph, i, from_to);
            to = GT_STRGRAPH_EDGE_DEST(strgraph, i, from_to);
            while (GT_STRGRAPH_V_IS_INTERNAL(strgraph, to))
            {
              from = to;
              from_to = gt_strgraph_find_only_edge(strgraph, from);
              GT_STRGRAPH_EDGE_SET_MARK(strgraph, from, from_to);
              to = GT_STRGRAPH_EDGE_DEST(strgraph, from, from_to);
            }
          }
          else
          {
            prev = info + p;
          }
        }
      }
    }
    if (show_progressbar)
      progress++;
  }
  counter = gt_strgraph_reduce_marked_edges(strgraph);

  if (show_progressbar)
    gt_progressbar_stop();
  gt_free(info);
  gt_log_log("p-bubbles = "GT_WU"", nofpbubbles);
  gt_log_log("removed p-bubble edges = "GT_WU"", counter);
#ifndef NDEBUG
  gt_strgraph_check_outdegs(strgraph);
#endif
  return counter;
}

#define GT_STRGRAPH_DOT_HEADER     "digraph StringGraph {\n"
#define GT_STRGRAPH_DOT_BI_HEADER  "graph StringGraph {\n"
#define GT_STRGRAPH_DOT_BC_HEADER  GT_STRGRAPH_DOT_BI_HEADER\
                                   "concentrate=true\n"
#define GT_STRGRAPH_DOT_FOOTER     "}\n"

#define GT_STRGRAPH_DOT_VSHAPE_INTERNAL    "ellipse"
#define GT_STRGRAPH_DOT_VSHAPE_JUNCTION    "box"
#define GT_STRGRAPH_DOT_VSHAPE_END         "triangle"

#define GT_STRGRAPH_DOT_SELECT_VSHAPE(STRGRAPH, V) \
  (GT_STRGRAPH_V_IS_INTERNAL(STRGRAPH, V) ? \
    GT_STRGRAPH_DOT_VSHAPE_INTERNAL : \
    (GT_STRGRAPH_V_IS_JUNCTION(STRGRAPH, V) ? \
      GT_STRGRAPH_DOT_VSHAPE_JUNCTION : \
      GT_STRGRAPH_DOT_VSHAPE_END))

#define GT_STRGRAPH_DOT_VCOLOR_GROUP_1   "palegreen"
#define GT_STRGRAPH_DOT_VCOLOR_GROUP_2   "lightskyblue"
#define GT_STRGRAPH_DOT_VCOLOR_MARKED    "moccasin"

#define GT_STRGRAPH_DOT_VCOLOR_ATTRS(C) "style=filled,fillcolor="C

/* color selection rules:
 * if depth > 1: generic color,
 * otherwise if marked: marked color,
 * otherwise if group is 1: group 1 color,
 * otherwise group 2 color */
#define GT_STRGRAPH_DOT_SELECT_VCOLOR(STRGRAPH, V, DEPTH, GROUP) \
  ((DEPTH) > 1UL ? "" : \
   (GT_STRGRAPH_V_MARK(STRGRAPH, V) == GT_STRGRAPH_V_MARKED ? \
    GT_STRGRAPH_DOT_VCOLOR_ATTRS(GT_STRGRAPH_DOT_VCOLOR_MARKED) : \
    ((GROUP) == 1UL ? \
     GT_STRGRAPH_DOT_VCOLOR_ATTRS(GT_STRGRAPH_DOT_VCOLOR_GROUP_1) : \
     GT_STRGRAPH_DOT_VCOLOR_ATTRS(GT_STRGRAPH_DOT_VCOLOR_GROUP_2))))

static inline void gt_strgraph_dot_show_vertex(GtFile *outfp,
    GtStrgraphFormat format, GtStrgraph *strgraph, GtStrgraphVnum v,
    GtUword depth, GtUword group)
{
  const char *shape = GT_STRGRAPH_DOT_SELECT_VSHAPE(strgraph, v),
             *color = GT_STRGRAPH_DOT_SELECT_VCOLOR(strgraph, v, depth, group);
  if (format == GT_STRGRAPH_DOT)
  {
    gt_file_xprintf(outfp, " \""GT_WU"%c\" [shape=%s,%s]\n",
        GT_STRGRAPH_V_READNUM(v), GT_STRGRAPH_V_CHAR(v), shape, color);
  }
  else
  {
    gt_assert(format == GT_STRGRAPH_DOT_BI);
    gt_file_xprintf(outfp, " "GT_WU" [shape=%s,%s]\n",
        GT_STRGRAPH_V_READNUM(v), shape, color);
  }
}

static inline void gt_strgraph_dot_show_edge(GtFile *outfp,
    GtStrgraphVnum from, GtStrgraphVnum to, GtStrgraphLength length)
{
  gt_file_xprintf(outfp,
      " \""GT_WU"%c\" -> \""GT_WU"%c\" "
      "[label="FormatGtStrgraphLength"];\n",
      GT_STRGRAPH_V_READNUM(from),
      GT_STRGRAPH_V_CHAR(from),
      GT_STRGRAPH_V_READNUM(to),
      GT_STRGRAPH_V_CHAR(to),
      PRINTGtStrgraphLengthcast(length));
}

static inline void gt_strgraph_dot_bi_show_edge(GtFile *outfp,
    GtUword sn1, bool towards1, GtUword sn2,
    bool towards2)
{
  gt_file_xprintf(outfp, " "GT_WU" -- "GT_WU
                  " [arrowtail=%s,arrowhead=%s,dir=both];\n", sn1, sn2,
                  towards1 ? "normal" : "inv", towards2 ? "normal" : "inv");
}

static void gt_strgraph_dot_show(const GtStrgraph *strgraph, GtFile *outfp,
    bool show_progressbar)
{
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  gt_file_xprintf(outfp, GT_STRGRAPH_DOT_HEADER);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++) {
    for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
    {
      if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
      {
        gt_strgraph_dot_show_edge(outfp, i,
            GT_STRGRAPH_EDGE_DEST(strgraph, i, j),
            GT_STRGRAPH_EDGE_LEN(strgraph, i, j));
      }
    }
    if (show_progressbar)
      progress++;
  }
  gt_file_xprintf(outfp, GT_STRGRAPH_DOT_FOOTER);

  if (show_progressbar)
    gt_progressbar_stop();
}

static void gt_strgraph_dot_bi_show(const GtStrgraph *strgraph, GtFile *outfp,
    bool show_progressbar)
{
  GtUword sn1, sn2;
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j;
  bool is_e1, is_e2;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  gt_file_xprintf(outfp, GT_STRGRAPH_DOT_BI_HEADER);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++) {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      sn1 = GT_STRGRAPH_V_READNUM(i);
      is_e1 = GT_STRGRAPH_V_IS_E(i);
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        {
          sn2 = GT_STRGRAPH_V_READNUM(GT_STRGRAPH_EDGE_DEST(strgraph, i, j));
          is_e2 = GT_STRGRAPH_V_IS_E(GT_STRGRAPH_EDGE_DEST(strgraph, i, j));
          if (is_e1 && is_e2) {
            gt_strgraph_dot_bi_show_edge(outfp, sn1, false, sn2, true);
          }
          else if ((is_e1 && !is_e2) && (sn1 < sn2)) {
            gt_strgraph_dot_bi_show_edge(outfp, sn1, false, sn2, false);
          }
          else if ((!is_e1 && is_e2) && (sn1 > sn2)) {
            gt_strgraph_dot_bi_show_edge(outfp, sn1, true, sn2, true);
          }
          /* other E->B / B->E cases, as well as all B-B
             are not considered to avoid double output */
        }
      }
    }
    if (show_progressbar)
      progress++;
  }
  gt_file_xprintf(outfp, GT_STRGRAPH_DOT_FOOTER);

  if (show_progressbar)
    gt_progressbar_stop();
}

typedef struct {
  GtStrgraphVnum v;
  GtUword d;
} GtStrgraphVnumAndDepth;

#define GT_STRGRAPH_DOT_VISITED(V) \
  (GT_STRGRAPH_V_MARK(strgraph, V) == GT_STRGRAPH_V_ELIMINATED ||\
   GT_STRGRAPH_V_MARK(strgraph, V) == GT_STRGRAPH_V_MARKED)

#define GT_STRGRAPH_DOT_SET_VISITED(V, D) \
  ((D) == 1UL ? \
    GT_STRGRAPH_V_SET_MARK(strgraph, V, GT_STRGRAPH_V_MARKED) : \
    GT_STRGRAPH_V_SET_MARK(strgraph, V, GT_STRGRAPH_V_ELIMINATED))

static int gt_strgraph_dot_show_context(GtStrgraph *strgraph,
    GtStrgraphFormat format, GtFile *outfp, GtUword *readnums,
    GtUword nofreadnums, GtUword maxdepth, bool extend,
    GtUword group, GtError *err)
{
  int had_err = 0;
  GtArray *stack;
  GtStrgraphVnumAndDepth to_add, *current;
  GtStrgraphVnum v, w;
  GtStrgraphVEdgenum j;
  GtUword d;
  GtUword sn1, sn2, readnum, i;
  bool is_e1, is_e2;

  gt_assert(format == GT_STRGRAPH_DOT || format == GT_STRGRAPH_DOT_BI);
  stack = gt_array_new(sizeof (GtStrgraphVnumAndDepth));
  for (i = 0; i < nofreadnums; i++)
  {
    readnum = readnums[i];
    if (readnum >= (GtUword)GT_STRGRAPH_NOFREADS(strgraph))
    {
      had_err = -1;
      gt_error_set(err, "Can't show context of read "GT_WU" "
          "because the readset has "FormatGtStrgraphVnum" reads", readnum,
          GT_STRGRAPH_NOFREADS(strgraph));
      break;
    }
    to_add.v = GT_STRGRAPH_V_B(readnum);
    to_add.d = 1UL;
    gt_array_add(stack, to_add);
    gt_strgraph_dot_show_vertex(outfp, format, strgraph, to_add.v, 1UL, group);
    to_add.v = GT_STRGRAPH_V_E(readnum);
    to_add.d = 1UL;
    gt_array_add(stack, to_add);
    gt_strgraph_dot_show_vertex(outfp, format, strgraph, to_add.v, 1UL, group);
  }
  if (had_err)
  {
    gt_array_delete(stack);
    return had_err;
  }
  while (gt_array_size(stack) > 0)
  {
    current = gt_array_pop(stack);
    v = current->v;
    d = current->d;
    if (!GT_STRGRAPH_DOT_VISITED(v))
    {
      sn1 = GT_STRGRAPH_V_READNUM(v);
      if (d > 1UL)
      {
        gt_strgraph_dot_show_vertex(outfp, format, strgraph, v, d, 0);
      }
      if (GT_STRGRAPH_V_OUTDEG(strgraph, v) > 0)
      {
        is_e1 = GT_STRGRAPH_V_IS_E(v);
        for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, v); j++)
        {
          if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, v, j))
          {
            w = GT_STRGRAPH_EDGE_DEST(strgraph, v, j);
            sn2 = GT_STRGRAPH_V_READNUM(w);
            is_e2 = GT_STRGRAPH_V_IS_E(w);
            if (format == GT_STRGRAPH_DOT)
            {
              gt_strgraph_dot_show_edge(outfp, v, w,
                  GT_STRGRAPH_EDGE_LEN(strgraph, v, j));
            }
            else
            {
              gt_strgraph_dot_bi_show_edge(outfp, sn1, !is_e1, sn2, is_e2);
            }
            if (d < maxdepth ||
                (extend && GT_STRGRAPH_V_IS_INTERNAL(strgraph, v)))
            {
              if (!GT_STRGRAPH_DOT_VISITED(w))
              {
                to_add.v = w;
                to_add.d = d + 1;
                gt_array_add(stack, to_add);
              }
              w = GT_STRGRAPH_V_OTHER(w);
              if (!GT_STRGRAPH_DOT_VISITED(w))
              {
                to_add.v = w;
                to_add.d = d + 1;
                gt_array_add(stack, to_add);
              }
            }
          }
        }
      }
    }
    GT_STRGRAPH_DOT_SET_VISITED(v, d);
  }
  gt_array_delete(stack);
  return had_err;
}

static void gt_strgraph_output_connecting_path(GtStrgraph *strgraph,
    GtArray *path, GtContigsWriter *cw)
{
  GtUword path_size, i;
  GtStrgraphVnumAndDepth *vd;
  path_size = gt_array_size(path);
  vd = gt_array_get(path, 0);
  gt_contigs_writer_start(cw,
      GT_STRGRAPH_V_MIRROR_SEQNUM(GT_STRGRAPH_NOFVERTICES(strgraph), vd->v));
  for (i = (GtUword)1UL; i < path_size; i++)
  {
    vd = gt_array_get(path, i);
    gt_contigs_writer_append(cw, GT_STRGRAPH_V_MIRROR_SEQNUM(
          GT_STRGRAPH_NOFVERTICES(strgraph), vd->v), (GtUword)vd->d);
  }
  gt_contigs_writer_write(cw);
}

/* return code: 0 if goal_sn was found, -1 otherwise */
static int gt_strgraph_add_edge_to_connecting_path(GtStrgraph *strgraph,
    GtStrgraphVnum from, GtStrgraphVEdgenum j, GtUword goal_sn,
    GtStrgraphVtype goal_vt, GtUword current_length, GtUword minlen,
    GtUword maxlen, GtArray *path, bool first_path_only,
    GtUword *count_too_long, GtUword *count_too_short, GtUword *count_found,
    GtUword *count_circular, GtUword *count_redundant, GtContigsWriter *cw,
    GtLogger *logger)
{
  GtStrgraphVnumAndDepth wl;
  GtStrgraphVEdgenum e;
  wl.v = GT_STRGRAPH_EDGE_DEST(strgraph, from, j);
  wl.d = GT_STRGRAPH_EDGE_LEN(strgraph, from, j);
  current_length += wl.d;
  if (current_length > maxlen)
  {
    (*count_too_long)++;
    return -1;
  }
  if ((GT_STRGRAPH_V_READNUM(wl.v) == goal_sn) &&
      ((goal_vt == GT_STRGRAPH_VTYPE_A) ||
       ((goal_vt == GT_STRGRAPH_VTYPE_E) && GT_STRGRAPH_V_IS_E(wl.v)) ||
       ((goal_vt == GT_STRGRAPH_VTYPE_B) && GT_STRGRAPH_V_IS_B(wl.v))))
  {
    if (current_length <= minlen)
    {
      (*count_too_short)++;
      return -1;
    }
    else
    {
      (*count_found)++;
      gt_array_add(path, wl);
      gt_logger_log(logger, "Path found, length: "GT_WU, current_length);
      gt_strgraph_output_connecting_path(strgraph, path, cw);
      gt_array_pop(path);
      return 0;
    }
  }
  if (GT_STRGRAPH_V_MARK(strgraph, wl.v) == GT_STRGRAPH_V_MARKED)
  {
    (*count_circular)++;
    return -1;
  }
  else if (GT_STRGRAPH_V_MARK(strgraph, wl.v) == GT_STRGRAPH_V_INPLAY)
  {
    (*count_redundant)++;
    return -1;
  }
  else
  {
    int found = -1;
    GT_STRGRAPH_V_SET_MARK(strgraph, wl.v, GT_STRGRAPH_V_MARKED);
    for (e = 0; e < GT_STRGRAPH_V_NOFEDGES(strgraph, wl.v); e++)
    {
      if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, wl.v, e))
      {
        int retcode;
        gt_array_add(path, wl);
        retcode = gt_strgraph_add_edge_to_connecting_path(strgraph, wl.v, e,
            goal_sn, goal_vt, current_length, minlen, maxlen, path,
            first_path_only, count_too_long, count_too_short, count_found,
            count_circular, count_redundant, cw, logger);
        gt_array_pop(path);
        if (retcode == 0)
        {
          if (first_path_only)
            return retcode;
          found = retcode;
        }
      }
    }
    GT_STRGRAPH_V_SET_MARK(strgraph, wl.v, GT_STRGRAPH_V_INPLAY);
    return found;
  }
}

/* return code: 0 if goal_sn was found, -1 otherwise */
static int gt_strgraph_find_connecting_path_from_vertex(GtStrgraph *strgraph,
    GtStrgraphVnumAndDepth vd, GtStrgraphVEdgenum nofedges, GtUword to,
    GtStrgraphVtype to_vt, GtUword minlen, GtUword maxlen, bool first_path_only,
    GtContigsWriter *cw, GtLogger *logger)
{
  GtArray *path;
  GtStrgraphVEdgenum j;
  int found = -1;
  GtUword count_too_long = 0, count_too_short = 0,
          count_found = 0, count_circular = 0,
          count_redundant = 0;

  path = gt_array_new(sizeof (GtStrgraphVnumAndDepth));
  gt_array_add(path, vd);
  for (j = 0; j < nofedges; j++)
  {
    if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, vd.v, j))
    {
      int retcode;
      retcode = gt_strgraph_add_edge_to_connecting_path(strgraph, vd.v, j,
          to, to_vt, vd.d, minlen, maxlen, path, first_path_only,
          &count_too_long, &count_too_short, &count_found, &count_circular,
          &count_redundant, cw, logger);
      if (retcode == 0)
      {
        found = retcode;
        if (first_path_only)
          break;
      }
    }
  }
  gt_array_delete(path);
  gt_logger_log(logger, "Paths interrupted as too long:  "GT_WU,
      count_too_long);
  gt_logger_log(logger, "Paths interrupted as circular:  "GT_WU,
      count_circular);
  gt_logger_log(logger, "Paths interrupted as redundant: "GT_WU,
      count_redundant);
  gt_logger_log(logger, "Paths discarded as too short:   "GT_WU,
      count_too_short);
  gt_logger_log(logger, "Paths output:                   "GT_WU,
      count_found);
  return found;
}

/* note: the algorithm is heuristic and generally will not find all
   connecting paths; e.g. if a vertex vA is reached by some path P1
   of length L1 and all paths from vA to the goal have a length
   > maxlen - L1, the vertex vA will be marked and not visited anymore;
   however, it could be that arriving later to vA through some other
   path P2 of length L2 < L1, the same paths would not be too long
   anymore - but it's too late, as vA has been marked; to fix this
   one would have to use, instead of a mark, a length value (which would
   make the algorithm similar to Dijekstra shortest paths)
*/
int gt_strgraph_find_connecting_path(GtStrgraph *strgraph, GtUword from,
    GtStrgraphVtype from_vt, GtUword to, GtStrgraphVtype to_vt,
    GtUword minlen, GtUword maxlen, bool first_path_only,
    const char *indexname, const char *suffix,
    GtLogger *logger, GtError *err)
{
  int had_err = 0;

  if (from >= (GtUword)GT_STRGRAPH_NOFREADS(strgraph))
  {
    had_err = -1;
    gt_error_set(err, "Can't search path from read "GT_WU" "
        "because the readset has "FormatGtStrgraphVnum" reads", from,
        GT_STRGRAPH_NOFREADS(strgraph));
  }
  if (to >= (GtUword)GT_STRGRAPH_NOFREADS(strgraph))
  {
    had_err = -1;
    gt_error_set(err, "Can't search path from read "GT_WU" "
        "because the readset has "FormatGtStrgraphVnum" reads", to,
        GT_STRGRAPH_NOFREADS(strgraph));
  }
  if (!had_err)
  {
    GtStrgraphVnumAndDepth vd;
    GtStrgraphVEdgenum nofedges, nofedges_sum = 0;
    GtStrgraphVnum i;
    GtContigsWriter *cw;
    GtFile *outfp;
    GtStr *complete_suffix;
    bool unreachable = false;
    int found = -1;

    /* check if "to" is reachable */
    if (to_vt == GT_STRGRAPH_VTYPE_B || to_vt == GT_STRGRAPH_VTYPE_A)
    {
      nofedges = GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_B(to));
      nofedges_sum = nofedges;
      if (nofedges == 0)
      {
        gt_logger_log(logger, "Destination read has no edges from B vertex");
        if (to_vt == GT_STRGRAPH_VTYPE_B)
          unreachable = true;
      }
    }
    if (to_vt == GT_STRGRAPH_VTYPE_E || to_vt == GT_STRGRAPH_VTYPE_A)
    {
      nofedges = GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_E(to));
      nofedges_sum += nofedges;
      if (nofedges == 0)
      {
        gt_logger_log(logger, "Destination read has no edges from E vertex");
        if (to_vt == GT_STRGRAPH_VTYPE_E || nofedges_sum == 0)
          unreachable = true;
      }
    }
    if (unreachable)
    {
      gt_logger_log(logger, "Destination read unreachable");
      return 0;
    }

    gt_assert(strgraph->encseq != NULL);
    gt_assert(gt_encseq_is_mirrored(strgraph->encseq));
    complete_suffix = gt_str_new();
    gt_str_append_char(complete_suffix, '.');
    gt_str_append_uword(complete_suffix, from);
    gt_str_append_char(complete_suffix, '.');
    gt_str_append_uword(complete_suffix, to);
    gt_str_append_cstr(complete_suffix, suffix);
    outfp = gt_strgraph_get_file(indexname, gt_str_get(complete_suffix),
        true, false);
    gt_logger_log(logger, "Connecting paths output filename: %s%s", indexname,
        gt_str_get(complete_suffix));
    cw = gt_contigs_writer_new(strgraph->encseq, outfp);
    gt_contigs_writer_enable_complete_path_output(cw);

    /* reset vertex marks */
    for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
      GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);

    vd.d = GT_STRGRAPH_SEQLEN(strgraph, from);

    /* disallow returning to origin read in any orientation
       (however, from == to should still work) */
    GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_V_B(from),
        GT_STRGRAPH_V_MARKED);
    GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_V_E(from),
        GT_STRGRAPH_V_MARKED);

    /* search path from B vertex */
    if (from_vt == GT_STRGRAPH_VTYPE_B || from_vt == GT_STRGRAPH_VTYPE_A)
    {
      vd.v = GT_STRGRAPH_V_B(from);
      nofedges = GT_STRGRAPH_V_NOFEDGES(strgraph, vd.v);
      if (nofedges == 0)
        gt_logger_log(logger, "Origin read has no edges from B vertex");
      else
      {
        gt_logger_log(logger, "Computing paths starting from B vertex...");
        found = gt_strgraph_find_connecting_path_from_vertex(strgraph, vd,
            nofedges, to, to_vt, minlen, maxlen, first_path_only, cw, logger);
      }
    }

    if ((from_vt == GT_STRGRAPH_VTYPE_E || from_vt == GT_STRGRAPH_VTYPE_A) &&
        ((found != 0) || !first_path_only))
    {
      /* search path from E vertex */
      vd.v = GT_STRGRAPH_V_E(from);
      nofedges = GT_STRGRAPH_V_NOFEDGES(strgraph, vd.v);
      if (nofedges == 0)
        gt_logger_log(logger, "Origin read has no edges from E vertex");
      else
      {
        gt_logger_log(logger, "Computing paths starting from E vertex...");
        gt_strgraph_find_connecting_path_from_vertex(strgraph, vd, nofedges,
            to, to_vt, minlen, maxlen, first_path_only, cw, logger);
      }
    }
    gt_str_delete(complete_suffix);
    gt_contigs_writer_delete(cw);
    gt_file_delete(outfp);
  }
  return had_err;
}

/* format: read[E|B] [numofedges]: neighbor1[E|B], neighbor2[E|B]... */
static void gt_strgraph_adjlist_show(const GtStrgraph *strgraph, GtFile *outfp)
{
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j;

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    gt_file_xprintf(outfp, ""GT_WU"%c "
        "["FormatGtStrgraphVEdgenum"]",
        GT_STRGRAPH_V_READNUM(i),
        GT_STRGRAPH_V_CHAR(i),
        PRINTGtStrgraphVEdgenumcast(GT_STRGRAPH_V_OUTDEG(strgraph, i)));

    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      gt_file_xprintf(outfp, " --> ");
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        {
          gt_file_xprintf(outfp, ""GT_WU"%c ("FormatGtStrgraphLength");",
              GT_STRGRAPH_V_READNUM(GT_STRGRAPH_EDGE_DEST(strgraph, i, j)),
              GT_STRGRAPH_V_CHAR(GT_STRGRAPH_EDGE_DEST(strgraph, i, j)),
              PRINTGtStrgraphVEdgenumcast(
                GT_STRGRAPH_EDGE_LEN(strgraph, i, j)));
        }
      }
    }
    gt_file_xfputc('\n', outfp);
  }
}

static void gt_strgraph_spm_show(const GtStrgraph *strgraph, GtFile *outfp)
{
  GtUword sn1, sn2;
  GtStrgraphVnum i, v2;
  GtStrgraphLength spm_len;
  GtStrgraphVEdgenum j;
  bool is_e1, is_e2;

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      sn1 = GT_STRGRAPH_V_READNUM(i);
      is_e1 = GT_STRGRAPH_V_IS_E(i);
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        {
          v2 = GT_STRGRAPH_EDGE_DEST(strgraph, i, j);
          sn2 = GT_STRGRAPH_V_READNUM(v2);
          spm_len = GT_STRGRAPH_SEQLEN(strgraph, sn2) -
            GT_STRGRAPH_EDGE_LEN(strgraph, i, j);
          is_e2 = GT_STRGRAPH_V_IS_E(v2);
          if ((is_e1 && is_e2) || (is_e1 && !is_e2 && (sn1 <= sn2)) ||
              (!is_e1 &&  is_e2 && (sn1 >= sn2)))
          {
            /* other E->B / B->E cases, as well as all B-B
               are not considered to avoid double output */
            gt_file_xprintf(outfp, GT_WU" %s " GT_WU " %s "
                            FormatGtStrgraphLength"\n", sn1, is_e1 ? "+" : "-",
                            sn2, is_e2 ? "+" : "-",
                PRINTGtStrgraphLengthcast(spm_len));
          }
        }
      }
    }
  }
}

static void gt_strgraph_asqg_show(const GtStrgraph *strgraph,
    const char *indexname, GtFile *outfp)
{
  GtUword sn1, sn2;
  GtStrgraphVnum i, v2;
  GtUword sl2, spm_len;
  GtStrgraphVEdgenum j;
  bool is_e1, is_e2;
  GtError *err;
  int had_err = 0;
  GtAsqgWriter *aw = NULL;

  err = gt_error_new();
  gt_assert(strgraph->encseq != NULL);
  aw = gt_asqg_writer_new(outfp, strgraph->encseq);
  had_err = gt_asqg_writer_show_header(aw, 0.0, (GtUword)
      strgraph->minmatchlen, indexname, false, false, err);
  if (!had_err)
    had_err = gt_asqg_writer_show_vertices(aw, err);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph) && !had_err; i++)
  {
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      sn1 = GT_STRGRAPH_V_READNUM(i);
      is_e1 = GT_STRGRAPH_V_IS_E(i);
      for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i) && !had_err; j++)
      {
        if (!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
        {
          v2 = GT_STRGRAPH_EDGE_DEST(strgraph, i, j);
          sn2 = GT_STRGRAPH_V_READNUM(v2);
          sl2 = (GtUword)GT_STRGRAPH_SEQLEN(strgraph, sn2);
          spm_len = sl2 - (GtUword)GT_STRGRAPH_EDGE_LEN(strgraph, i, j);
          is_e2 = GT_STRGRAPH_V_IS_E(v2);
          if ((is_e1 && is_e2 && (sn1 >= sn2)) ||
              (!is_e1 && !is_e2 && (sn1 > sn2)) ||
              (is_e1 && !is_e2 && (sn1 >= sn2)) ||
              (!is_e1 &&  is_e2 && (sn1 > sn2)))
          {
            /* other E->B / B->E cases, as well as all B-B
               are not considered to avoid double output */
            gt_spmproc_show_asgq(sn1, sn2, spm_len, is_e1, is_e2, aw);
          }
        }
      }
    }
  }
  if (had_err)
  {
    fprintf(stderr, "%s", gt_error_get(err));
    exit(EXIT_FAILURE);
  }
  gt_asqg_writer_delete(aw);
  gt_error_delete(err);
}

void gt_strgraph_show(const GtStrgraph *strgraph, GtStrgraphFormat format,
    const char *indexname, const char *suffix, bool show_progressbar)
{
  GtFile *outfp = NULL;

  gt_assert(strgraph != NULL);
  outfp = gt_strgraph_get_file(indexname, suffix, true,
      format == GT_STRGRAPH_ASQG_GZ ? true : false);
  switch (format)
  {
    case GT_STRGRAPH_DOT:
      gt_strgraph_dot_show(strgraph, outfp, show_progressbar);
      break;
    case GT_STRGRAPH_DOT_BI:
      gt_strgraph_dot_bi_show(strgraph, outfp, show_progressbar);
      break;
    case GT_STRGRAPH_ADJLIST:
      gt_strgraph_adjlist_show(strgraph, outfp);
      break;
    case GT_STRGRAPH_SPM:
      gt_strgraph_spm_show(strgraph, outfp);
      break;
    case GT_STRGRAPH_ASQG_GZ: /*@ fallthrough @*/
    case GT_STRGRAPH_ASQG:
      gt_strgraph_asqg_show(strgraph, indexname, outfp);
      break;
    case GT_STRGRAPH_BIN:
      gt_strgraph_save(strgraph, outfp);
      break;
    default:
      gt_assert(false);
      break;
  }
  gt_file_delete(outfp);
}

static int gt_strgraph_show_disc_distri_datapoint(GtUword key,
    GtUint64 value, GtFile *outfile)
{
  gt_file_xprintf(outfile, ""GT_WU" "Formatuint64_t"\n", key,
     PRINTuint64_tcast((uint64_t)value));
  return 0;
}

void gt_strgraph_show_edge_lengths_distribution(const GtStrgraph *strgraph,
    const char *indexname, const char *suffix)
{
  GtDiscDistri *d;
  GtStrgraphVnum i;
  GtStrgraphVEdgenum j;
  GtFile *outfp;

  gt_assert(strgraph != NULL);
  gt_assert(sizeof (GtUword) >= sizeof (GtStrgraphLength));
  d = gt_disc_distri_new();
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
    {
      gt_disc_distri_add(d,
          (GtUword)GT_STRGRAPH_EDGE_LEN(strgraph, i, j));
    }
  }
  outfp = gt_strgraph_get_file(indexname, suffix, true, false);
  gt_file_xprintf(outfp, "# length nofedges\n");
  gt_disc_distri_foreach(d,
     (GtDiscDistriIterFunc) gt_strgraph_show_disc_distri_datapoint, outfp);
  gt_disc_distri_delete(d);
  gt_file_delete(outfp);
}

void gt_strgraph_show_counts_distribution(const GtStrgraph *strgraph,
    const char *indexname, const char *suffix)
{
  GtDiscDistri *d;
  GtFile *outfp;
  GtStrgraphVnum i;

  gt_assert(sizeof (GtUword) >= sizeof (GtStrgraphCount));
  d = gt_disc_distri_new();
  gt_assert(strgraph != NULL &&\
            strgraph->__small_counts != NULL &&
            strgraph->__large_counts != NULL);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    GtStrgraphCount c;
    GT_STRGRAPH_GET_COUNT(strgraph, c, i);
    gt_disc_distri_add(d, (GtUword)c);
  }
  outfp = gt_strgraph_get_file(indexname, suffix, true, false);
  gt_file_xprintf(outfp, "# count nofvertices\n");
  gt_disc_distri_foreach(d,
     (GtDiscDistriIterFunc) gt_strgraph_show_disc_distri_datapoint, outfp);
  gt_disc_distri_delete(d);
  gt_file_delete(outfp);
}

int gt_strgraph_show_context(GtStrgraph *strgraph, GtStrgraphFormat format,
    const char *indexname, const char *suffix, GtUword *readnums,
    GtUword nofreadnums, GtUword *otherreadnums,
    GtUword nofotherreadnums, GtUword maxdepth, bool extend,
    GtError *err)
{
  int had_err = 0;
  GtFile *outfp = NULL;
  GtStrgraphVnum i;

  gt_assert(strgraph != NULL);
  gt_assert(format == GT_STRGRAPH_DOT || format == GT_STRGRAPH_DOT_BI);

  outfp = gt_strgraph_get_file(indexname, suffix, true, false);
  gt_file_xprintf(outfp, format == GT_STRGRAPH_DOT ?
      GT_STRGRAPH_DOT_HEADER : GT_STRGRAPH_DOT_BC_HEADER);
  gt_assert(sizeof (GtStrgraphVnum) >= sizeof (GtUword));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);

  had_err = gt_strgraph_dot_show_context(strgraph, format, outfp, readnums,
      nofreadnums, maxdepth, extend, 1UL, err);
  if (!had_err && nofotherreadnums > 0)
    had_err = gt_strgraph_dot_show_context(strgraph, format, outfp,
        otherreadnums, nofotherreadnums, maxdepth, extend, 2UL, err);
  gt_file_xprintf(outfp, GT_STRGRAPH_DOT_FOOTER);
  gt_file_delete(outfp);
  return had_err;
}

/* log using for large values */
static void gt_strgraph_log_space_value(const char *prefix, size_t value)
{
  char unit_prefix;
  float float_value;
  if (value < (size_t)1024) {
    gt_log_log("%s = "GT_WU" bytes", prefix, (GtUword)value);
  } else {
    if (value < (size_t)1048576) {
      unit_prefix = 'k';
      float_value = (float)value / 1024.F;
    } else if (value < (size_t)1073741824UL) {
      unit_prefix = 'M';
      float_value = (float)value / 1048576.F;
    } else {
      unit_prefix = 'G';
      float_value = (float)value / 1073741824.F;
    }
    gt_log_log("%s = "GT_WU" bytes (%.2f %cb)", prefix, (GtUword)value,
        float_value, unit_prefix);
  }
}

void gt_strgraph_log_space(const GtStrgraph *strgraph)
{
  size_t space_for_edges, space_for_vertices, total_space;
  gt_assert(strgraph != NULL);
  space_for_edges = GT_STRGRAPH_SIZEOF_EDGES(strgraph);
  space_for_vertices = GT_STRGRAPH_SIZEOF_VERTICES(strgraph);
  total_space = sizeof (GtStrgraph) + space_for_edges + space_for_vertices;
  gt_strgraph_log_space_value("space graph", total_space);
  gt_strgraph_log_space_value("- edges", space_for_edges);
  gt_strgraph_log_space_value("- vertices", space_for_vertices);
}

static inline GtStrgraphVnum gt_strgraph_nofvertices_Vnum(
    const GtStrgraph *strgraph)
{
  GtStrgraphVnum i, connected_vertices;

  gt_assert(strgraph != NULL);

  connected_vertices = 0;
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if ((GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0) ||
        (GT_STRGRAPH_V_INDEG(strgraph, i) > 0))
    {
      connected_vertices++;
    }
  }
  return connected_vertices;
}

GtUword gt_strgraph_nofvertices(const GtStrgraph *strgraph)
{
  GtStrgraphVnum retval;
  retval = gt_strgraph_nofvertices_Vnum(strgraph);
  gt_assert(sizeof (GtUword) >= sizeof (GtStrgraphVnum) ||
      retval <= (GtStrgraphVnum)ULONG_MAX);
  return (GtUword)retval;
}

GtUword gt_strgraph_nofreads(const GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtUword connected_reads = 0;
  bool already_counted = false;

  gt_assert(strgraph != NULL);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_IS_B(i))
      already_counted = false;
      /* if IS_E, then the read may already have been counted,
         in the case B was connected*/
    if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0)
    {
      if (already_counted == false)
      {
        connected_reads++;
        already_counted = true;
      }
    }
  }
  return connected_reads;
}

static inline GtStrgraphEdgenum gt_strgraph_nofedges_Edgenum(
    const GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtStrgraphEdgenum total_degree = 0;
  gt_assert(strgraph != NULL);
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    total_degree += GT_STRGRAPH_V_OUTDEG(strgraph, i);
  return total_degree;
}

GtUword gt_strgraph_nofedges(const GtStrgraph *strgraph)
{
  GtStrgraphEdgenum retval;
  retval = gt_strgraph_nofedges_Edgenum(strgraph);
  gt_assert(sizeof (GtUword) >= sizeof (GtStrgraphEdgenum) ||
      retval <= (GtStrgraphEdgenum)ULONG_MAX);
  return (GtUword)retval;
}

GtUword gt_strgraph_nofspm(const GtStrgraph *strgraph)
{
  GtUword total_degree;
  gt_assert(strgraph != NULL);
  /* nofspm is half of total degree */
  total_degree = gt_strgraph_nofedges(strgraph);
  /*  gt_assert(!(total_degree & 1)); */
  return ((total_degree >> 1));
}

void gt_strgraph_log_stats(const GtStrgraph *strgraph, GtLogger *logger)
{
  GtStrgraphVnum v = 0;
  GtStrgraphEdgenum e = 0;
  float avg_deg;

  gt_assert(strgraph != NULL);

  v = gt_strgraph_nofvertices_Vnum(strgraph);
  e = gt_strgraph_nofedges_Edgenum(strgraph);

  avg_deg = (v == 0UL) ? (float)0 : (float)(e) / (float)(v);

  gt_logger_log(logger, "vertices: "FormatGtStrgraphVnum
      " (reads: "FormatGtStrgraphVnum") -- edges: "FormatGtStrgraphEdgenum" "
      "(spm: "FormatGtStrgraphEdgenum") -- e/v: %.4f", v, v >> 1, e, e >> 1,
      avg_deg);
}

static inline void gt_strgraph_traverse_simple_path(GtStrgraph *strgraph,
    GtStrgraphVnum i, GtStrgraphVEdgenum j, void(*process_edge)
    (GtStrgraphVnum, GtStrgraphLength, void*), void *data)
{
  GtStrgraphVnum from, to;
  GtStrgraphVEdgenum from_to;

  from = i;
  from_to = j;
  to = GT_STRGRAPH_EDGE_DEST(strgraph, i, j);

  while (GT_STRGRAPH_V_IS_INTERNAL(strgraph, to) && i != to &&
      GT_STRGRAPH_V_MARK(strgraph, to) != GT_STRGRAPH_V_ELIMINATED)
  {
    if (process_edge != NULL)
      process_edge(to, GT_STRGRAPH_EDGE_LEN(strgraph, from, from_to), data);
    GT_STRGRAPH_V_SET_MARK(strgraph, to, GT_STRGRAPH_V_ELIMINATED);
    GT_STRGRAPH_V_SET_MARK(strgraph, GT_STRGRAPH_V_OTHER(to),
        GT_STRGRAPH_V_ELIMINATED);
    from = to;
    from_to = gt_strgraph_find_only_edge(strgraph, from);
    to = GT_STRGRAPH_EDGE_DEST(strgraph, from, from_to);
  }
  if (process_edge != NULL)
    process_edge(to, GT_STRGRAPH_EDGE_LEN(strgraph, from, from_to), data);
}

static inline void gt_strgraph_traverse_from_vertex(GtStrgraph *strgraph,
    GtStrgraphVnum i, void(*process_start) (GtStrgraphVnum, void*),
    void(*process_edge) (GtStrgraphVnum, GtStrgraphLength, void*), void *data)
{
  GtStrgraphVEdgenum j;
  gt_assert(GT_STRGRAPH_V_OUTDEG(strgraph, i) > 0);
  for (j = 0; j < GT_STRGRAPH_V_NOFEDGES(strgraph, i); j++)
  {
    if (GT_STRGRAPH_EDGE_IS_REDUCED(strgraph, i, j))
      continue;
    if (GT_STRGRAPH_V_MARK(strgraph, GT_STRGRAPH_EDGE_DEST(strgraph, i, j))
        == GT_STRGRAPH_V_ELIMINATED)
      continue;
    if (process_start != NULL)
      process_start(i, data);
    gt_strgraph_traverse_simple_path(strgraph, i, j, process_edge, data);
  }
}

#ifdef GG_DEBUG
static void gt_strgraph_count_junctions(GtStrgraph *strgraph)
{
  GtStrgraphVnum i;
  GtUword nof_out1_junctions  = 0,
                nof_in1_junctions   = 0,
                nof_multi_junctions = 0;
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i+= (GtStrgraphVnum)2)
  {
    if (GT_STRGRAPH_V_IS_JUNCTION(strgraph, i))
    {
      if (GT_STRGRAPH_V_OUTDEG(strgraph, i) > (GtStrgraphVEdgenum)1 && \
          GT_STRGRAPH_V_INDEG(strgraph, i) == (GtStrgraphVEdgenum)1)
      {
        nof_in1_junctions++;
      }
      else if (GT_STRGRAPH_V_INDEG(strgraph, i) > (GtStrgraphVEdgenum)1 && \
          GT_STRGRAPH_V_OUTDEG(strgraph, i) == (GtStrgraphVEdgenum)1)
      {
        nof_out1_junctions++;
      }
      else
      {
        nof_multi_junctions++;
      }
    }
  }
  gt_log_log("junctions: in-1:"GT_WU" out-1:"GT_WU" multi:"GT_WU,
             nof_in1_junctions, nof_out1_junctions, nof_multi_junctions);
}
#endif

static void gt_strgraph_traverse(GtStrgraph *strgraph,
    void(*process_start) (GtStrgraphVnum, void*),
    void(*process_edge) (GtStrgraphVnum, GtStrgraphLength, void*),
    void *data, bool show_progressbar)
{
  GtStrgraphVnum i;
  GtUint64 progress = 0;

  gt_assert(strgraph != NULL);

#ifdef GG_DEBUG
  gt_strgraph_count_junctions(strgraph);
#endif

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
    GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_VACANT);

  if (show_progressbar)
    gt_progressbar_start(&progress,
        (GtUint64)GT_STRGRAPH_NOFVERTICES(strgraph));

  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_MARK(strgraph, i) != GT_STRGRAPH_V_ELIMINATED)
    {
      if (GT_STRGRAPH_V_OUTDEG(strgraph, i) == 0)
      {
        GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_ELIMINATED);
      }
      else if (!GT_STRGRAPH_V_IS_INTERNAL(strgraph, i))
      {
        gt_strgraph_traverse_from_vertex(strgraph, i, process_start,
            process_edge, data);
        GT_STRGRAPH_V_SET_MARK(strgraph, i, GT_STRGRAPH_V_ELIMINATED);
      }
    }
    if (show_progressbar)
      progress++;
  }

  if (show_progressbar)
    gt_progressbar_stop();

  /* handle circles of internal vertices only */
  for (i = 0; i < GT_STRGRAPH_NOFVERTICES(strgraph); i++)
  {
    if (GT_STRGRAPH_V_IS_INTERNAL(strgraph, i) &&
        GT_STRGRAPH_V_MARK(strgraph, i) != GT_STRGRAPH_V_ELIMINATED)
      gt_strgraph_traverse_from_vertex(strgraph, i, process_start,
          process_edge, data);
  }
}

#define GT_STRGRAPH_CONTIG_INC 16384UL

/* --- Contig Paths Output --- */

GT_DECLAREARRAYSTRUCT(GtContigpathElem);

typedef struct {
  GtUword                 total_depth, current_depth, contignum,
                          min_depth, jnum;
  GtStrgraphVnum          nof_v, firstnode, lastnode;
  FILE                    *p_file, *cjl_i_file, *cjl_o_file, *ji_file;
  GtArrayGtContigpathElem contig;
  GtStrgraph              *strgraph;
  bool                    show_contigs_info;
} GtStrgraphContigpathsData;

static void gt_strgraph_show_contigpath_edge(GtStrgraphVnum v,
    GtStrgraphLength len, void *data)
{
  GtStrgraphContigpathsData *pdata = data;
  GtContigpathElem *value;
  GtUword seqnum;

  (pdata->current_depth)++;

  /* first element: length */
  GT_GETNEXTFREEINARRAY(value, &pdata->contig, GtContigpathElem,
      GT_STRGRAPH_CONTIG_INC);
  gt_assert(sizeof (GtContigpathElem) >= sizeof (GtStrgraphLength) ||
      len <= (GtStrgraphLength)GT_CONTIGPATH_ELEM_MAX);
  *value = (GtContigpathElem)len;

  /* second element: vertex number */
  GT_GETNEXTFREEINARRAY(value, &pdata->contig, GtContigpathElem,
      GT_STRGRAPH_CONTIG_INC);
  seqnum = GT_STRGRAPH_V_MIRROR_SEQNUM(pdata->nof_v, v);
  gt_assert(sizeof (GtContigpathElem) >= sizeof (GtStrgraphVnum) ||
      seqnum <= (GtUword)GT_CONTIGPATH_ELEM_MAX);
  *value = (GtContigpathElem)seqnum;

  pdata->lastnode = v;
}

static inline void gt_strgraph_show_contiginfo(GtStrgraphContigpathsData *pdata)
{
  GtContigEdgesLink cjl;
  GtContigJunctionInfo ji;

  cjl.deg = (uint64_t)GT_STRGRAPH_V_INDEG(pdata->strgraph, pdata->firstnode);
  cjl.ptr = (uint64_t)GT_STRGRAPH_V_OTHER(pdata->firstnode);
  (void)gt_xfwrite(&cjl, sizeof (cjl), (size_t)1, pdata->cjl_i_file);

  cjl.deg = (uint64_t)GT_STRGRAPH_V_OUTDEG(pdata->strgraph, pdata->lastnode);
  cjl.ptr = (uint64_t)pdata->lastnode;
  (void)gt_xfwrite(&cjl, sizeof (cjl), (size_t)1, pdata->cjl_o_file);

  ji.contig_num = (uint32_t)pdata->contignum;
  if (GT_STRGRAPH_V_IS_JUNCTION(pdata->strgraph, pdata->firstnode))
  {
    ji.junction_num = (GtUword)pdata->firstnode;
    ji.firstnode = 1U;
    ji.length = (unsigned int)GT_STRGRAPH_SEQLEN(pdata->strgraph,
          GT_STRGRAPH_V_READNUM(pdata->firstnode));
    (void)gt_xfwrite(&ji, sizeof (ji), (size_t)1, pdata->ji_file);
    pdata->jnum++;
  }
  if (GT_STRGRAPH_V_IS_JUNCTION(pdata->strgraph, pdata->lastnode))
  {
    ji.junction_num = (GtUword)(GT_STRGRAPH_V_OTHER(pdata->lastnode));
    ji.firstnode = 0;
    ji.length = (unsigned int)GT_STRGRAPH_SEQLEN(pdata->strgraph,
          GT_STRGRAPH_V_READNUM(pdata->lastnode));
    (void)gt_xfwrite(&ji, sizeof (ji), (size_t)1, pdata->ji_file);
    pdata->jnum++;
  }
}

static inline void gt_strgraph_end_contigpath(GtStrgraphContigpathsData *pdata)
{
  (void)gt_xfwrite(pdata->contig.spaceGtContigpathElem,
    sizeof (GtContigpathElem), (size_t)pdata->contig.nextfreeGtContigpathElem,
    pdata->p_file);
  pdata->total_depth += pdata->current_depth;

  if (pdata->show_contigs_info)
    gt_strgraph_show_contiginfo(pdata);
  (pdata->contignum)++;
}

static inline void gt_strgraph_begin_contigpath(
    GtStrgraphContigpathsData *pdata, GtStrgraphVnum firstvertex)
{
  GtContigpathElem *value;
  GtUword seqnum;

  pdata->current_depth = 1UL;
  pdata->firstnode = firstvertex;

  /* 0 to start a new contig */
  pdata->contig.nextfreeGtContigpathElem = 0;
  GT_GETNEXTFREEINARRAY(value, &pdata->contig, GtContigpathElem,
      GT_STRGRAPH_CONTIG_INC);
  *value = 0;

  /* seqnum of origin vertex */
  GT_GETNEXTFREEINARRAY(value, &pdata->contig, GtContigpathElem,
      GT_STRGRAPH_CONTIG_INC);
  seqnum = GT_STRGRAPH_V_MIRROR_SEQNUM(pdata->nof_v, firstvertex);
  gt_assert(sizeof (GtContigpathElem) >= sizeof (GtStrgraphVnum) ||
      seqnum <= (GtUword)GT_CONTIGPATH_ELEM_MAX);
  *value = (GtContigpathElem)seqnum;
}

static void gt_strgraph_show_contigpath_vertex(GtStrgraphVnum firstvertex,
    void *data)
{
  GtStrgraphContigpathsData *pdata = data;

  if (pdata->current_depth >= pdata->min_depth)
    gt_strgraph_end_contigpath(pdata);

  gt_strgraph_begin_contigpath(pdata, firstvertex);
}

static void gt_strgraph_show_contigpaths(GtStrgraph *strgraph,
    GtUword min_path_depth, FILE *p_file, FILE *cjl_i_file,
    FILE *cjl_o_file, FILE *ji_file, bool show_contigs_info,
    bool show_progressbar)
{
  GtStrgraphContigpathsData pdata;

  pdata.total_depth = 1UL;
  pdata.current_depth = 1UL;
  pdata.contignum = 0;
  pdata.jnum = 0;

  pdata.lastnode = 0;
  GT_INITARRAY(&pdata.contig, GtContigpathElem);
  pdata.p_file = p_file;
  pdata.show_contigs_info = show_contigs_info;
  pdata.cjl_i_file = cjl_i_file;
  pdata.cjl_o_file = cjl_o_file;
  pdata.ji_file = ji_file;
  pdata.min_depth = min_path_depth;
  pdata.nof_v = GT_STRGRAPH_NOFVERTICES(strgraph);
  pdata.strgraph = strgraph;

  /* leave space for header */
  if (show_contigs_info)
  {
    gt_xfseek(ji_file, (GtWord) sizeof (pdata.jnum), SEEK_SET);
    gt_xfseek(cjl_i_file, (GtWord) sizeof (pdata.contignum), SEEK_SET);
    gt_xfseek(cjl_o_file, (GtWord) sizeof (pdata.contignum), SEEK_SET);
  }

  gt_strgraph_traverse(strgraph, gt_strgraph_show_contigpath_vertex,
      gt_strgraph_show_contigpath_edge, &pdata, show_progressbar);

  /* show last contig path */
  if (pdata.current_depth >= pdata.min_depth)
    gt_strgraph_end_contigpath(&pdata);

  /* write header */
  if (show_contigs_info)
  {
    gt_xfseek(ji_file, 0, SEEK_SET);
    (void)gt_xfwrite(&pdata.jnum, sizeof (pdata.jnum), (size_t)1, ji_file);
    gt_xfseek(cjl_i_file, 0, SEEK_SET);
    (void)gt_xfwrite(&pdata.contignum, sizeof (pdata.contignum),
        (size_t)1, cjl_i_file);
    gt_xfseek(cjl_o_file, 0, SEEK_SET);
    (void)gt_xfwrite(&pdata.contignum, sizeof (pdata.contignum),
        (size_t)1, cjl_o_file);
  }

  gt_log_log("traversed edges = "GT_WU"", pdata.total_depth);
  gt_log_log("numofcontigs = "GT_WU"", pdata.contignum);

  GT_FREEARRAY(&pdata.contig, GtContigpathElem);
}

/* --- Direct Contig Output --- */

typedef struct {
  GtUword total_depth, current_depth, min_depth,
                current_length, min_length, contignum;
  GtStrgraph *strgraph;
  GtContigsWriter *cw;
} GtStrgraphSpellData;

static void gt_strgraph_spell_edge(GtStrgraphVnum v, GtStrgraphLength len,
    void *data)
{

  GtStrgraphSpellData *sdata = data;
  gt_contigs_writer_append(sdata->cw, GT_STRGRAPH_V_MIRROR_SEQNUM(
        GT_STRGRAPH_NOFVERTICES(sdata->strgraph), v), (GtUword)len);
  (sdata->current_depth)++;
  sdata->current_length += len;
}

static void gt_strgraph_spell_vertex(GtStrgraphVnum firstvertex, void *data)
{
  GtStrgraphSpellData *sdata = data;

  if ((sdata->current_depth >= sdata->min_depth) &&
      (sdata->current_length >= sdata->min_length))
  {
    sdata->total_depth += sdata->current_depth;
    gt_contigs_writer_write(sdata->cw);
    (sdata->contignum)++;
  }
  else
    gt_contigs_writer_abort(sdata->cw);

  gt_contigs_writer_start(sdata->cw,
      GT_STRGRAPH_V_MIRROR_SEQNUM(GT_STRGRAPH_NOFVERTICES(sdata->strgraph),
        firstvertex));
  sdata->current_length = (GtUword)GT_STRGRAPH_SEQLEN(sdata->strgraph,
      firstvertex);
  sdata->current_depth = 1UL;
}

static void gt_strgraph_show_contigs(GtStrgraph *strgraph,
    GtUword min_path_depth, GtUword min_contig_length,
    bool showpaths, GtFile *outfp, const GtEncseq *encseq,
    bool show_progressbar, GtLogger *logger)
{
  GtStrgraphSpellData sdata;

  gt_strgraph_set_encseq(strgraph, encseq);
  sdata.cw = gt_contigs_writer_new(encseq, outfp);
  if (showpaths)
    gt_contigs_writer_enable_complete_path_output(sdata.cw);
  sdata.strgraph = strgraph;
  sdata.total_depth = 1UL;
  sdata.current_depth = 1UL;
  sdata.current_length = 0;
  sdata.contignum = 0;
  sdata.min_depth = min_path_depth;
  sdata.min_length = min_contig_length;

  gt_strgraph_traverse(strgraph, gt_strgraph_spell_vertex,
      gt_strgraph_spell_edge, &sdata, show_progressbar);

  /* show last contig */
  if (sdata.current_depth >= min_path_depth &&
      sdata.current_length >= min_contig_length)
    gt_contigs_writer_write(sdata.cw);
  else
    gt_contigs_writer_abort(sdata.cw);

  gt_log_log("traversed edges = "GT_WU"", sdata.total_depth);
  gt_log_log("numofcontigs = "GT_WU"", sdata.contignum);

  if (sdata.contignum > 0)
    gt_contigs_writer_show_stats(sdata.cw, logger);
  else
    gt_logger_log(logger, "no contigs respect the given cutoff parameters");

  gt_contigs_writer_delete(sdata.cw);
}

void gt_strgraph_spell(GtStrgraph *strgraph, GtUword min_path_depth,
    GtUword min_contig_length, bool showpaths, const char *indexname,
    const char *suffix, const GtEncseq *encseq, bool delay_reads_mapping,
    bool show_contigs_info, bool show_progressbar, GtLogger *logger)
{
  FILE *main_file = NULL;
  GtStr *filename;

  gt_assert(strgraph != NULL);
  filename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(filename, suffix);

  main_file = gt_fa_xfopen(gt_str_get(filename), "w");

  if (!delay_reads_mapping)
  {
    GtFile *gt_outfp = gt_file_new_from_fileptr(main_file);
    gt_strgraph_show_contigs(strgraph, min_path_depth, min_contig_length,
        showpaths, gt_outfp, encseq, show_progressbar, logger);
    gt_file_delete_without_handle(gt_outfp);
  }
  else if (show_contigs_info)
  {
    FILE *cjl_i_file, *cjl_o_file, *ji_file;
    gt_str_set(filename, indexname);
    gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_CJ_I_LINKS);
    cjl_i_file = gt_fa_xfopen(gt_str_get(filename), "w");
    gt_str_set(filename, indexname);
    gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_CJ_O_LINKS);
    cjl_o_file = gt_fa_xfopen(gt_str_get(filename), "w");
    gt_str_set(filename, indexname);
    gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_JUNCTIONS);
    ji_file = gt_fa_xfopen(gt_str_get(filename), "w");
    gt_strgraph_show_contigpaths(strgraph, min_path_depth, main_file,
       cjl_i_file, cjl_o_file, ji_file, true, show_progressbar);
    gt_fa_xfclose(cjl_i_file);
    gt_fa_xfclose(cjl_o_file);
    gt_fa_xfclose(ji_file);
  }
  else
  {
    gt_strgraph_show_contigpaths(strgraph, min_path_depth, main_file,
       NULL, NULL, NULL, false, show_progressbar);
  }

  gt_str_delete(filename);
  gt_fa_xfclose(main_file);
}

/* --- UNIT TESTS --- */

#define GT_STRGRAPH_UTEST(TEST)  \
  if (!had_err)        \
    had_err = gt_strgraph_ ## TEST ## _unit_test(err)

static int gt_strgraph_new_unit_test(GtError *err)
{
  int had_err = 0;
  GtStrgraph *strgraph = NULL;

  gt_error_check(err);

  /* ordinary number of reads */
  strgraph = gt_strgraph_new(100UL);
  gt_ensure(strgraph != NULL);
  gt_strgraph_delete(strgraph);

  return had_err;
}

static int gt_strgraph_delete_unit_test(GT_UNUSED /* when compiled without
  assertions */ GtError *err)
{
  GtStrgraph *strgraph;

  gt_error_check(err);

  strgraph = NULL;
  gt_strgraph_delete(strgraph);

  strgraph = gt_strgraph_new(100UL);
  gt_strgraph_delete(strgraph);

  return 0;
}

static int gt_strgraph_creation_unit_test(GtError *err)
{
  int had_err = 0;
  GtStrgraph *strgraph;
  GtStrgraphCount c;
  GtUword nofreads = 2UL;

  gt_error_check(err);
  strgraph = gt_strgraph_new(nofreads);
  gt_spmproc_strgraph_count(0, 1UL, 10UL, true, true, strgraph);
  GT_STRGRAPH_GET_COUNT(strgraph, c, GT_STRGRAPH_V_B(0));
  gt_ensure(c == 0);
  GT_STRGRAPH_GET_COUNT(strgraph, c, GT_STRGRAPH_V_E(0));
  gt_ensure(c == (GtStrgraphCount)1);
  GT_STRGRAPH_GET_COUNT(strgraph, c, GT_STRGRAPH_V_B(1));
  gt_ensure(c == (GtStrgraphCount)1);
  GT_STRGRAPH_GET_COUNT(strgraph, c, GT_STRGRAPH_V_E(1));
  gt_ensure(c == 0);
  gt_strgraph_allocate_graph(strgraph, 100UL, NULL);
  gt_ensure(GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_B(0)) == 0);
  gt_ensure(GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_E(0)) ==
      (GtStrgraphVEdgenum)1);
  gt_ensure(GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_B(1)) ==
      (GtStrgraphVEdgenum)1);
  gt_ensure(GT_STRGRAPH_V_NOFEDGES(strgraph, GT_STRGRAPH_V_E(1)) == 0);
  gt_strgraph_delete(strgraph);
  return had_err;
}

static int gt_strgraph_add_spm_unit_test(GtError *err)
{
  int had_err = 0;
  GtStrgraph *strgraph;
  GtUword nofreads = 2UL;

  gt_error_check(err);
  strgraph = gt_strgraph_new(nofreads);
  gt_spmproc_strgraph_count(0, 1UL, 10UL, true, true, strgraph);
  gt_strgraph_allocate_graph(strgraph, 22UL, NULL);
  gt_spmproc_strgraph_add(0, 1UL, 10UL, true, true, strgraph);
  gt_ensure(GT_STRGRAPH_NOFVERTICES(strgraph) == (GtStrgraphVnum)4);
  gt_ensure(GT_STRGRAPH_V_OUTDEG(strgraph, GT_STRGRAPH_V_B(0)) == 0);
  gt_ensure(GT_STRGRAPH_V_OUTDEG(strgraph, GT_STRGRAPH_V_E(0)) ==
      (GtStrgraphVEdgenum)1);
  gt_ensure(GT_STRGRAPH_EDGE_DEST(strgraph, GT_STRGRAPH_V_E(0), 0) ==
      GT_STRGRAPH_V_E(1));
  gt_ensure(GT_STRGRAPH_EDGE_LEN(strgraph, GT_STRGRAPH_V_E(0), 0) ==
      (GtStrgraphLength)12);
  gt_ensure(!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph,
      GT_STRGRAPH_V_E(0), 0));
  gt_ensure(GT_STRGRAPH_V_OUTDEG(strgraph, GT_STRGRAPH_V_B(1)) ==
      (GtStrgraphVEdgenum)1);
  gt_ensure(GT_STRGRAPH_V_OUTDEG(strgraph, GT_STRGRAPH_V_E(1)) == 0);
  gt_ensure(GT_STRGRAPH_EDGE_DEST(strgraph, GT_STRGRAPH_V_B(1), 0) ==
      (GtStrgraphVnum)GT_STRGRAPH_V_B(0));
  gt_ensure(GT_STRGRAPH_EDGE_LEN(strgraph, GT_STRGRAPH_V_B(1), 0) ==
      (GtStrgraphLength)12);
  gt_ensure(!GT_STRGRAPH_EDGE_IS_REDUCED(strgraph,
      GT_STRGRAPH_V_B(1), 0));
  gt_strgraph_delete(strgraph);
  return had_err;
}

static int gt_strgraph_redtrans_unit_test(GtError *err)
{
  int had_err = 0;
  GtStrgraph *strgraph;
  GtUword nofreads = 5UL;
  GT_ENSURE_OUTPUT_DECLARE(2000);

  gt_error_check(err);

  /*
  test case:
                          --------2-------------
                   ----------0-----------
             --------RC(1)---------
          ---------3------------
       -------RC(4)----------
  */

  strgraph = gt_strgraph_new(nofreads);
  gt_spmproc_strgraph_count(4UL, 3UL, 19UL, false, true, strgraph);
  gt_spmproc_strgraph_count(1UL, 4UL, 16UL, true, true, strgraph);
  gt_spmproc_strgraph_count(4UL, 0UL, 10UL, false, true, strgraph);
  gt_spmproc_strgraph_count(4UL, 2UL, 3UL, false, true, strgraph);
  gt_spmproc_strgraph_count(3UL, 1UL, 19UL, true, false, strgraph);
  gt_spmproc_strgraph_count(3UL, 0UL, 13UL, true, true, strgraph);
  gt_spmproc_strgraph_count(3UL, 2UL, 6UL, true, true, strgraph);
  gt_spmproc_strgraph_count(1UL, 0UL, 16UL, false, true, strgraph);
  gt_spmproc_strgraph_count(1UL, 2UL, 9UL, false, true, strgraph);
  gt_spmproc_strgraph_count(0UL, 2UL, 15UL, true, true, strgraph);
  gt_strgraph_allocate_graph(strgraph, 22UL, NULL);
  gt_spmproc_strgraph_add(4UL, 3UL, 19UL, false, true, strgraph);
  gt_spmproc_strgraph_add(1UL, 4UL, 16UL, true, true, strgraph);
  gt_spmproc_strgraph_add(4UL, 0UL, 10UL, false, true, strgraph);
  gt_spmproc_strgraph_add(4UL, 2UL, 3UL, false, true, strgraph);
  gt_spmproc_strgraph_add(3UL, 1UL, 19UL, true, false, strgraph);
  gt_spmproc_strgraph_add(3UL, 0UL, 13UL, true, true, strgraph);
  gt_spmproc_strgraph_add(3UL, 2UL, 6UL, true, true, strgraph);
  gt_spmproc_strgraph_add(1UL, 0UL, 16UL, false, true, strgraph);
  gt_spmproc_strgraph_add(1UL, 2UL, 9UL, false, true, strgraph);
  gt_spmproc_strgraph_add(0UL, 2UL, 15UL, true, true, strgraph);
  GT_ENSURE_OUTPUT(gt_strgraph_dot_show(strgraph, outfp, false),
      "digraph StringGraph {\n"
      " \"0B\" -> \"4E\" [label=12];\n"
      " \"0B\" -> \"3B\" [label=9];\n"
      " \"0B\" -> \"1E\" [label=6];\n"
      " \"0E\" -> \"2E\" [label=7];\n"
      " \"1B\" -> \"0E\" [label=6];\n"
      " \"1B\" -> \"2E\" [label=13];\n"
      " \"1E\" -> \"4E\" [label=6];\n"
      " \"1E\" -> \"3B\" [label=3];\n"
      " \"2B\" -> \"4E\" [label=19];\n"
      " \"2B\" -> \"3B\" [label=16];\n"
      " \"2B\" -> \"1E\" [label=13];\n"
      " \"2B\" -> \"0B\" [label=7];\n"
      " \"3B\" -> \"4E\" [label=3];\n"
      " \"3E\" -> \"1B\" [label=3];\n"
      " \"3E\" -> \"0E\" [label=9];\n"
      " \"3E\" -> \"2E\" [label=16];\n"
      " \"4B\" -> \"3E\" [label=3];\n"
      " \"4B\" -> \"1B\" [label=6];\n"
      " \"4B\" -> \"0E\" [label=12];\n"
      " \"4B\" -> \"2E\" [label=19];\n"
      "}\n"
    );
  gt_strgraph_sort_edges_by_len(strgraph, false);
  GT_ENSURE_OUTPUT(gt_strgraph_dot_show(strgraph, outfp, false),
      "digraph StringGraph {\n"
      " \"0B\" -> \"1E\" [label=6];\n"
      " \"0B\" -> \"3B\" [label=9];\n"
      " \"0B\" -> \"4E\" [label=12];\n"
      " \"0E\" -> \"2E\" [label=7];\n"
      " \"1B\" -> \"0E\" [label=6];\n"
      " \"1B\" -> \"2E\" [label=13];\n"
      " \"1E\" -> \"3B\" [label=3];\n"
      " \"1E\" -> \"4E\" [label=6];\n"
      " \"2B\" -> \"0B\" [label=7];\n"
      " \"2B\" -> \"1E\" [label=13];\n"
      " \"2B\" -> \"3B\" [label=16];\n"
      " \"2B\" -> \"4E\" [label=19];\n"
      " \"3B\" -> \"4E\" [label=3];\n"
      " \"3E\" -> \"1B\" [label=3];\n"
      " \"3E\" -> \"0E\" [label=9];\n"
      " \"3E\" -> \"2E\" [label=16];\n"
      " \"4B\" -> \"3E\" [label=3];\n"
      " \"4B\" -> \"1B\" [label=6];\n"
      " \"4B\" -> \"0E\" [label=12];\n"
      " \"4B\" -> \"2E\" [label=19];\n"
      "}\n"
    );
  (void)gt_strgraph_redtrans(strgraph, false);
  GT_ENSURE_OUTPUT(gt_strgraph_dot_show(strgraph, outfp, false),
      "digraph StringGraph {\n"
      " \"0B\" -> \"1E\" [label=6];\n"
      " \"0E\" -> \"2E\" [label=7];\n"
      " \"1B\" -> \"0E\" [label=6];\n"
      " \"1E\" -> \"3B\" [label=3];\n"
      " \"2B\" -> \"0B\" [label=7];\n"
      " \"3B\" -> \"4E\" [label=3];\n"
      " \"3E\" -> \"1B\" [label=3];\n"
      " \"4B\" -> \"3E\" [label=3];\n"
      "}\n"
    );

  gt_strgraph_delete(strgraph);
  return had_err;
}

int gt_strgraph_unit_test(GtError *err)
{
  int had_err = 0;

  gt_error_check(err);
  GT_STRGRAPH_UTEST(new);
  GT_STRGRAPH_UTEST(delete);
  GT_STRGRAPH_UTEST(creation);
  GT_STRGRAPH_UTEST(add_spm);
  GT_STRGRAPH_UTEST(redtrans);
  return had_err;
}
