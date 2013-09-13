/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/array_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "core/xansi_api.h"
#include "match/rdj-contig-info.h"
#include "match/rdj-contigs-graph.h"
#include "match/rdj-contigpaths.h"
#include "match/reads_libraries_table.h"

#define GT_CONTIGS_GRAPH_ASTAT_SINGLECOPY 17.0

#define GT_CONTIGS_GRAPH_IS_SINGLE_COPY(CG, CNUM) \
  ((CG)->v_d[CNUM].astat >= GT_CONTIGS_GRAPH_ASTAT_SINGLECOPY)

#define GT_CONTIGS_GRAPH_IS_OPTIONAL(CG, CNUM) \
  ((CG)->v_m[CNUM].optional)

#define GT_CONTIGS_GRAPH_EACH_SPM_EDGE(CG, CNUM, INCOMING, EDGE) \
  (EDGE) = (CG)->e_spm[INCOMING] + (CG)->v_spm[(INCOMING)][(CNUM)].ptr;\
  (EDGE) < (CG)->e_spm[INCOMING] + (CG)->v_spm[(INCOMING)][(CNUM)+1UL].ptr;\
  (EDGE)++

#define GT_CONTIGS_GRAPH_VALID_PATH_DIR(INCOMING, REVERSE)\
  (((bool)(INCOMING) == (REVERSE)) ? 0 : 1U)

#define GT_CONTIGS_GRAPH_OTHER_DIR(INCOMING)\
  ((INCOMING) == 0 ? 1U : 0)

/* pairing edges
   (draft) */
typedef struct {
  GtUword from,
                to,
                librarynum;
} GtContigsGraphScfEdge;

typedef struct {
  /* contig properties */
  bool deleted   : 1;
  bool optional  : 1;
  /* traversal markers (must be reset each time before traversal) */
  bool selected  : 1;
  bool processed : 1;
  bool visited   : 1;
  bool marked    : 1;
  bool mark0     : 1;
  bool mark1     : 1;
} GtContigsGraphMarks;

typedef struct {
  GtUword seqnum  : GT_INTWORDSIZE - 1;
  bool          revcomp : 1;
  GtUword offset;
} GtContigsGraphSeqUnit;

/* junction edge;

   implicit information coded by array position:
   - source contig number;
   - orientation at source (incoming edge if a prefix of source matches,
                            outgoing edge if a suffix of source matches)
*/
typedef struct {
  uint32_t      dest;          /* destination contig number */
  unsigned int  ovlen   : 29;  /* length of contig overlap */
  bool          reverse : 1;   /* edge orientation at destination:
                                  true if the destination contig has reverse
                                  orientation with respect to source contig */
  bool          deleted  : 1;
} GtContigsGraphSpmEdge;

typedef enum {
  GT_CONTIGS_GRAPH_JUNCTION,
  GT_CONTIGS_GRAPH_SINGLE_END,
  GT_CONTIGS_GRAPH_MULTIPLE_END,
  GT_CONTIGS_GRAPH_CIRCULAR,
} GtContigsGraphPathEndType;

typedef struct {
  GtUword dest;
  unsigned int dir;
  bool extended;
} GtContigsGraphPathElem;

GT_DECLAREARRAYSTRUCT(GtContigsGraphPathElem);

typedef struct {
  GtUword cnum, dest, depth;
  unsigned int  dir;
  GtContigsGraphPathEndType t;
} GtContigsGraphPathEndInfo;

struct GtContigsGraph
{
  /* vertices */
  GtUword            nof_v;
  GtUword            alloc_v;
  GtUword            nof_simple_v;
  GtContigEdgesLink        *v_spm[2]; /* link to spm edges:
                                         [0] = outgoing,
                                         [1] = incoming */
  GtContigEdgesLink        *v_scf; /* link to scaffolding edges */
  GtContigDepthInfo        *v_d;   /* depth information */
  GtContigsGraphMarks      *v_m;   /* marks */
  GtContigEdgesLink        *v_cmp; /* link to composition */
  /* composition units */
  GtUword            nof_units;
  GtUword            alloc_units;
  GtContigsGraphSeqUnit    *units;
  /* spm edges */
  GtUword            nof_spm_edges[2];
  GtUword            alloc_spm_edges[2];
  GtContigsGraphSpmEdge    *e_spm[2];
  /* scf edges */
  GtUword            nofscfedges;
  GtContigsGraphScfEdge    *e_scf;
  /* further information */
  GtReadsLibrariesTable    *rlt;
  /* misc */
  bool                     dot_show_deleted;
};

static int gt_contigs_graph_jinfo_cmp(const void *a, const void *b)
{
  return (int)
         (((GtContigJunctionInfo*)a)->junction_num -
          ((GtContigJunctionInfo*)b)->junction_num);
}

static GtContigsGraph* gt_contigs_graph_init(void)
{
  GtContigsGraph *cg;
  cg = gt_malloc(sizeof (GtContigsGraph));
  cg->nof_v = 0;
  cg->alloc_v = 0;
  cg->nof_simple_v = 0;
  cg->v_spm[0] = NULL;
  cg->v_spm[1] = NULL;
  cg->v_scf = NULL;
  cg->v_d   = NULL;
  cg->v_m   = NULL;
  cg->v_cmp = NULL;
  cg->nof_units = 0;
  cg->alloc_units = 0;
  cg->units = NULL;
  cg->nof_spm_edges[0] = 0;
  cg->nof_spm_edges[1] = 0;
  cg->alloc_spm_edges[0] = 0;
  cg->alloc_spm_edges[1] = 0;
  cg->e_spm[0] = NULL;
  cg->e_spm[1] = NULL;
  cg->nofscfedges = 0;
  cg->e_scf = NULL;
  cg->rlt = NULL;
  gt_log_log("sizeof (GtContigEdgesLink) = "GT_WU"",
      (GtUword) sizeof (GtContigEdgesLink));
  gt_log_log("sizeof (GtContigJunctionInfo) = "GT_WU"",
      (GtUword) sizeof (GtContigJunctionInfo));
  gt_log_log("sizeof (GtContigDepthInfo) = "GT_WU"",
      (GtUword) sizeof (GtContigDepthInfo));
  gt_log_log("sizeof (GtContigsGraphSpmEdge) = "GT_WU"",
      (GtUword) sizeof (GtContigsGraphSpmEdge));
  gt_log_log("sizeof (GtContigsGraphScfEdge) = "GT_WU"",
      (GtUword) sizeof (GtContigsGraphScfEdge));
  gt_log_log("sizeof (GtContigsGraphMarks) = "GT_WU"",
      (GtUword) sizeof (GtContigsGraphMarks));
  gt_log_log("sizeof (GtContigsGraphSeqUnit) = "GT_WU"",
      (GtUword) sizeof (GtContigsGraphSeqUnit));
  cg->dot_show_deleted = false;
  return cg;
}

void gt_contigs_graph_output_paths(GtContigsGraph *cg, FILE *fp)
{
  GtContigpathElem cnum;
  uint64_t unum;
  GtContigpathElem elems[2];
  GtContigsGraphSeqUnit *unit;
  gt_assert(cg != NULL);
  gt_assert(fp != NULL);
  elems[0] = 0;
  gt_log_log("nof_simple_v: "GT_WU"", cg->nof_simple_v);
  gt_log_log("nof_v: "GT_WU"", cg->nof_v);
  for (cnum = 0; cnum < (uint32_t)cg->nof_simple_v; cnum++)
  {
    if (cg->v_m[cnum].deleted || cg->v_m[cnum].optional)
      continue;
    elems[1] = cnum;
    (void)gt_xfwrite(&elems, sizeof (*elems), (size_t)2, fp);
  }
  for (cnum = (uint32_t)cg->nof_simple_v; cnum < (uint32_t)cg->nof_v; cnum++)
  {
    if (cg->v_m[cnum].deleted || cg->v_m[cnum].optional)
      continue;
    for (unum = 0; unum < cg->v_cmp[cnum - cg->nof_simple_v].deg; unum++)
    {
      unit = cg->units + cg->v_cmp[cnum - cg->nof_simple_v].ptr + unum;
      gt_assert(cg->v_d[unit->seqnum].length >= unit->offset);
      elems[0] = (unum == 0)
        ? 0
        : (uint32_t) cg->v_d[unit->seqnum].length - unit->offset;
      elems[1] = unit->revcomp
        ? (uint32_t) GT_MULT2(cg->nof_simple_v) - 1 - unit->seqnum
        : (uint32_t) unit->seqnum;
      (void)gt_xfwrite(&elems, sizeof (*elems), (size_t)2, fp);
    }
  }
}

void gt_contigs_graph_enable_dot_show_deleted(GtContigsGraph *cg)
{
  gt_assert(cg != NULL);
  cg->dot_show_deleted = true;
}

#define GT_CONTIGS_GRAPH_NEW_FREAD(PTR, COUNT, STREAM, ERRSTR1, ERRSTR2)\
  {\
    size_t count;\
    count = fread(PTR, sizeof (*(PTR)), (size_t)COUNT, STREAM);\
    if (count != (size_t)(COUNT))\
    {\
      gt_error_set(err, "error reading " ERRSTR1 " from " ERRSTR2 " file "\
          "(exp:"GT_WU"; found="GT_WU")", (GtUword)COUNT, (GtUword)count);\
      had_err = -1;\
    }\
  }

#ifdef GG_DEBUG
static void gt_contigs_graph_show_junctions(GtContigJunctionInfo *j_info,
    GtUword nofjinfos, const char *path)
{
  GtUword jnum;
  FILE *fp = fopen(path, "w");
  for (jnum = 0; jnum < nofjinfos; jnum++)
  {
    fprintf(fp, "junction="GT_WU" contig=%u\n",
        j_info[jnum].junction_num,
        j_info[jnum].contig_num);
  }
  (void)fclose(fp);
}

static void gt_contigs_graph_show_e_links(GtContigEdgesLink *el,
    GtUword nof_v, const char *path)
{
  GtUword vnum;
  FILE *fp = fopen(path, "w");
  for (vnum = 0; vnum < nof_v; vnum++)
  {
    fprintf(fp, "contignum="GT_WU" deg="GT_WU" ptr="GT_WU"\n", vnum,
        (GtUword)el[vnum].deg, (GtUword)el[vnum].ptr);
  }
  (void)fclose(fp);
}
#endif

#define GT_CONTIGS_GRAPH_V_INC 128UL
#define GT_CONTIGS_GRAPH_E_INC 256UL
#define GT_CONTIGS_GRAPH_U_INC 256UL

GtUword gt_contigs_graph_append_vertex(GtContigsGraph *cg,
    GtUword nof_spm_o, GtUword nof_spm_i,
    GtUword nof_units)
{
  unsigned int incoming;
  gt_assert(cg->nof_v <= cg->alloc_v);
  gt_log_log("append_vertex(nof_spm_o="GT_WU", nof_spm_i="GT_WU", nof_units="
             GT_WU")", nof_spm_o, nof_spm_i, nof_units);
  if (cg->nof_v == cg->alloc_v)
  {
    cg->alloc_v += GT_CONTIGS_GRAPH_V_INC;
    cg->v_spm[0] = gt_realloc(cg->v_spm[0],
        sizeof (*(cg->v_spm[0])) * (cg->alloc_v + 1UL));
    cg->v_spm[1] = gt_realloc(cg->v_spm[1],
        sizeof (*(cg->v_spm[1])) * (cg->alloc_v + 1UL));
    cg->v_cmp = gt_realloc(cg->v_cmp,
        sizeof (*cg->v_cmp) * (cg->alloc_v - cg->nof_simple_v));
    /* cg->v_scf = gt_realloc(cg->v_scf,
        sizeof (*cg->v_scf) * (cg->alloc_v + 1UL)); */
    cg->v_d = gt_realloc(cg->v_d, sizeof (*cg->v_d) * (cg->alloc_v));
    cg->v_m = gt_realloc(cg->v_m, sizeof (*cg->v_m) * (cg->alloc_v));
  }
  for (incoming = 0; incoming < 2U; incoming++)
  {
    GtUword nof_e = incoming ? nof_spm_i : nof_spm_o;
    GtContigsGraphSpmEdge *edge;
    cg->v_spm[incoming][cg->nof_v].deg = (uint64_t)nof_e;
    cg->v_spm[incoming][cg->nof_v + 1UL].ptr =
      cg->v_spm[incoming][cg->nof_v].ptr + nof_e;
    if (cg->nof_spm_edges[incoming] + nof_e > cg->alloc_spm_edges[incoming])
    {
      cg->alloc_spm_edges[incoming] += MAX(GT_CONTIGS_GRAPH_E_INC, nof_e);
      cg->e_spm[incoming] = gt_realloc(cg->e_spm[incoming],
          sizeof (*(cg->e_spm[incoming])) * cg->alloc_spm_edges[incoming]);
    }
    for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cg->nof_v, incoming, edge))
    {
      edge->deleted = false;
    }
  }
  cg->v_cmp[cg->nof_v - cg->nof_simple_v].deg = (uint64_t) nof_units;
  cg->v_cmp[cg->nof_v - cg->nof_simple_v].ptr =
    cg->nof_v == cg->nof_simple_v ? 0 :
    cg->v_cmp[cg->nof_v - 1 - cg->nof_simple_v].ptr +
    cg->v_cmp[cg->nof_v - 1 - cg->nof_simple_v].deg;
  if (nof_units > 0)
  {
    cg->nof_units += nof_units;
    if (cg->nof_units > cg->alloc_units)
    {
      cg->alloc_units += MAX(GT_CONTIGS_GRAPH_U_INC, nof_units);
      cg->units = gt_realloc(cg->units,
          sizeof (*cg->units) * (cg->alloc_units));
    }
    gt_assert(cg->nof_units <= cg->alloc_units);
  }
  cg->v_m[cg->nof_v].deleted = false;
  cg->v_m[cg->nof_v].optional = false;
  cg->nof_v++;
  return cg->nof_v - 1UL;
}

static int gt_contigs_graph_create_vertices(GtContigsGraph *cg,
    FILE *cjl_i_fp, FILE *cjl_o_fp, FILE *depthinfo_fp, GtError *err)
{
  int had_err = 0;
  GtUword vnum, nofv;
  gt_assert(cg != NULL);
  gt_assert(cjl_i_fp != NULL);
  gt_assert(cjl_o_fp != NULL);
  gt_assert(depthinfo_fp != NULL);
  gt_error_check(err);
  GT_CONTIGS_GRAPH_NEW_FREAD(&cg->nof_v, 1, cjl_i_fp,
      "number of contigs", "incoming contigs-junctions links file");
  if (!had_err)
  {
    GT_CONTIGS_GRAPH_NEW_FREAD(&nofv, 1, cjl_o_fp,
        "number of contigs", "outgoing contigs-junctions links file");
    if (nofv != cg->nof_v)
    {
      gt_error_set(err, "number of contigs differ in "
      "incoming ("GT_WU") and outgoing ("GT_WU") contigs-junctions links files",
      (GtUword)cg->nof_v, (GtUword)nofv);
    }
  }
  cg->alloc_v = cg->nof_v;
  cg->nof_simple_v = cg->nof_v;
  if (!had_err)
  {
    cg->v_spm[0] = gt_malloc(sizeof (*cg->v_spm[0]) * (cg->alloc_v + 1UL));
    GT_CONTIGS_GRAPH_NEW_FREAD(cg->v_spm[0], cg->nof_v, cjl_o_fp,
        "links", "outgoing contigs-junctions links");
  }
  if (!had_err)
  {
    cg->v_spm[1] = gt_malloc(sizeof (*cg->v_spm[1]) * (cg->alloc_v + 1UL));
    GT_CONTIGS_GRAPH_NEW_FREAD(cg->v_spm[1], cg->nof_v, cjl_i_fp,
        "links", "incoming contigs-junctions links file");
  }
  if (!had_err)
  {
    cg->v_d = gt_malloc(sizeof (*cg->v_d) * (cg->alloc_v));
    GT_CONTIGS_GRAPH_NEW_FREAD(cg->v_d, cg->nof_v, depthinfo_fp,
        "depth information", "depth info");
  }
  if (!had_err)
  {
    /* cg->v_scf = gt_malloc(sizeof (*cg->v_scf) * (cg->alloc_v + 1UL)); */
    cg->v_m = gt_malloc(sizeof (*cg->v_m) * (cg->alloc_v));
    for (vnum = 0; vnum < cg->nof_v; vnum++)
    {
      cg->v_m[vnum].deleted = false;
      cg->v_m[vnum].optional = (cg->v_d[vnum].depth == 2UL);
    }
#ifdef GG_DEBUG
    gt_contigs_graph_show_e_links(cg->v_spm[1], cg->nof_v, "debug.cjl_i");
    gt_contigs_graph_show_e_links(cg->v_spm[0], cg->nof_v, "debug.cjl_o");
#endif
  }
  return had_err;
}

static int gt_contigs_graph_read_junctions(GtContigJunctionInfo **j_info,
    GtUword *nofjinfos, FILE *junctions_fp, GtError *err)
{
  int had_err = 0;
  gt_assert(nofjinfos != NULL);
  gt_assert(j_info != NULL);
  gt_assert(junctions_fp != NULL);
  gt_error_check(err);
  GT_CONTIGS_GRAPH_NEW_FREAD(nofjinfos, 1, junctions_fp,
      "number of junctions", "junctions info");
  *j_info = gt_malloc(sizeof (**j_info) * (*nofjinfos));
  GT_CONTIGS_GRAPH_NEW_FREAD(*j_info, *nofjinfos, junctions_fp,
      "junctions", "junctions info");
  qsort(*j_info, (size_t)*nofjinfos, sizeof (**j_info),
      gt_contigs_graph_jinfo_cmp);
#ifdef GG_DEBUG
  gt_contigs_graph_show_junctions(*j_info, *nofjinfos,
      "debug.sorted_junctions");
#endif
  return had_err;
}

static void gt_contigs_graph_create_spm_edges_for_vertex(GtContigEdgesLink *v,
    GtContigsGraphSpmEdge *e, GtUword cnum, GtUword *nextfree_edge,
    GtContigJunctionInfo *j_info, GtUword nofjinfos, bool inwards)
{
  GtContigJunctionInfo *j_elem, key;
  uint32_t degree, j;
  key.junction_num = (GtUword)v[cnum].ptr;
  v[cnum].ptr = (uint64_t)(*nextfree_edge);
  degree = v[cnum].deg;
  if (degree > 0)
  {
    j_elem = bsearch(&key, j_info, (size_t)nofjinfos,
        sizeof (*j_info), gt_contigs_graph_jinfo_cmp);
    gt_assert(j_elem != NULL);
    while (j_elem != j_info && (j_elem - 1)->junction_num == key.junction_num)
      j_elem--;
    for (j = 0; j < degree; j++, (*nextfree_edge)++, j_elem++)
    {
      gt_assert(j_elem->junction_num == key.junction_num);
      e[*nextfree_edge].dest = j_elem->contig_num;
      e[*nextfree_edge].ovlen = (unsigned int) j_elem->length;
      e[*nextfree_edge].reverse = inwards
        ? (j_elem->firstnode == 1U ? true : false)
        : (j_elem->firstnode == 1U ? false : true);
      e[*nextfree_edge].deleted = false;
    }
  }
}

static void gt_contigs_graph_create_e_spm(GtUword nof_v,
    GtContigEdgesLink *v, GtContigsGraphSpmEdge **e,
    GtUword *nof_e, GtUword *alloc_e, GtContigJunctionInfo *j_info,
    GtUword nofjinfos, bool incoming)
{
  GtUword cnum, nextfree_edge;
  *nof_e = 0;
  for (cnum = 0; cnum < nof_v; cnum++)
  {
    *nof_e += v[cnum].deg;
  }
  *alloc_e = *nof_e;

  *e = gt_malloc(sizeof (**e) * (*alloc_e));
  nextfree_edge = 0;
  for (cnum = 0; cnum < nof_v; cnum++)
  {
    gt_contigs_graph_create_spm_edges_for_vertex(v, *e,
        cnum, &nextfree_edge, j_info, nofjinfos, incoming);
  }
  v[nof_v].ptr = (uint64_t)nextfree_edge;
}

static int gt_contigs_graph_load_reads_libraries_table(GtContigsGraph *cg,
    FILE *rlt_fp, GtError *err)
{
  gt_assert(cg != NULL);
  gt_assert(rlt_fp != NULL);
  cg->rlt = gt_reads_libraries_table_load(rlt_fp, err);
  return (cg->rlt == NULL) ? -1 : 0;
}

GtContigsGraph* gt_contigs_graph_new(FILE *cjl_i_fp, FILE *cjl_o_fp,
    FILE *junctions_fp, FILE *rlt_fp, FILE *depthinfo_fp, GtError *err)
{
  GtContigsGraph *cg;
  GtUword nofjinfos;
  GtContigJunctionInfo *j_info = NULL;
  int had_err = 0;

  cg = gt_contigs_graph_init();

  had_err = gt_contigs_graph_load_reads_libraries_table(cg, rlt_fp, err);

  if (!had_err)
    had_err = gt_contigs_graph_create_vertices(cg, cjl_i_fp, cjl_o_fp,
        depthinfo_fp, err);

  if (!had_err)
    had_err = gt_contigs_graph_read_junctions(&j_info, &nofjinfos,
        junctions_fp, err);

  if (!had_err)
  {
    unsigned int incoming;
    for (incoming = 0; incoming < 2U; incoming++)
    {
      gt_contigs_graph_create_e_spm(cg->nof_v, cg->v_spm[incoming],
          &(cg->e_spm[incoming]), &(cg->nof_spm_edges[incoming]),
          &(cg->alloc_spm_edges[incoming]), j_info, nofjinfos, (bool)incoming);
    }
  }

  gt_free(j_info);
  if (had_err)
  {
    gt_contigs_graph_delete(cg);
    return NULL;
  }
  return cg;
}

#define GT_CONTIGS_GRAPH_DOT_HEADER "graph ContigsGraph {\n"
#define GT_CONTIGS_GRAPH_DOT_FOOTER "}\n"

#define GT_CONTIGS_GRAPH_DOT_SELECT_V_SHAPE(CG, CNUM)\
  (GT_CONTIGS_GRAPH_IS_OPTIONAL(CG, CNUM) ? "diamond" : "ellipse")

#define GT_CONTIGS_GRAPH_DOT_SELECT_V_STYLE(CG, CNUM)\
  ((CG)->v_m[CNUM].deleted ? "dotted" : "solid")

#define GT_CONTIGS_GRAPH_DOT_SELECT_V_COLOR(CG, CNUM)\
  (GT_CONTIGS_GRAPH_IS_SINGLE_COPY(CG, CNUM) ? "red" : "black")

static void gt_contigs_graph_show_dot_vertex(GtContigsGraph *cg,
    GtUword cnum, GtFile *outfp)
{
  const char *shape = GT_CONTIGS_GRAPH_DOT_SELECT_V_SHAPE(cg, cnum),
             *style = GT_CONTIGS_GRAPH_DOT_SELECT_V_STYLE(cg, cnum),
             *color = GT_CONTIGS_GRAPH_DOT_SELECT_V_COLOR(cg, cnum);

  gt_file_xprintf(outfp, "  "GT_WU" [style=%s,color=%s,shape=%s];\n", cnum,
        style, color, shape);
}

static void gt_contigs_graph_show_dot_for_contig(GtContigsGraph *cg,
    GtUword cnum, GtFile *outfp)
{
  unsigned int incoming;
  GtContigsGraphSpmEdge *edge;
  for (incoming = 0; incoming < 2U; incoming++)
  {
    for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
    {
      if (!cg->v_m[edge->dest].processed)
      {
        if (cg->dot_show_deleted || !edge->deleted)
          gt_file_xprintf(outfp, "  "GT_WU" -- "GT_WU" "
              "[dir=both,arrowtail=%s,arrowhead=%s%s];\n", cnum,
              (GtUword)edge->dest,
              incoming ? "normal" : "inv",
              ((incoming && edge->reverse) ||
               (!incoming && !edge->reverse)) ? "normal" : "inv",
              edge->deleted ? ",style=dotted" : "");
      }
    }
  }
  if (cg->dot_show_deleted || !cg->v_m[cnum].deleted)
    gt_contigs_graph_show_dot_vertex(cg, cnum, outfp);
}

void gt_contigs_graph_show_dot(GtContigsGraph *cg, GtFile *outfp)
{
  GtUword cnum;
  gt_assert(cg != NULL);
  gt_file_xprintf(outfp, GT_CONTIGS_GRAPH_DOT_HEADER);
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    cg->v_m[cnum].processed = false;
  }
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    gt_contigs_graph_show_dot_for_contig(cg, cnum, outfp);
    cg->v_m[cnum].processed = true;
  }
  gt_file_xprintf(outfp, GT_CONTIGS_GRAPH_DOT_FOOTER);
}

GtContigsGraphSpmEdge* gt_contigs_graph_find_only_spm_edge(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming)
{
  GtContigsGraphSpmEdge *edge;
  gt_log_log("find_only_spm_edge(cnum="GT_WU",incoming=%u)",cnum,incoming);
  for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
  {
    if (!edge->deleted)
      return edge;
  }
  gt_assert(false);
  return NULL;
}

GtContigsGraphSpmEdge* gt_contigs_graph_find_spm_edge(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming, GtUword dest)
{
  GtContigsGraphSpmEdge *edge;
  gt_log_log("find_spm_edge(cnum="GT_WU",incoming=%u,dest="GT_WU")",
             cnum,incoming,dest);
  for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
  {
    if (edge->deleted)
      continue;
    if (edge->dest == (uint64_t)dest)
      return edge;
  }
  gt_assert(false);
  return NULL;
}

GtContigsGraphSpmEdge* gt_contigs_graph_find_deleted_spm_edge(
    GtContigsGraph *cg, GtUword cnum, unsigned int incoming)
{
  GtContigsGraphSpmEdge *edge;
  gt_log_log("find_deleted_spm_edge(cnum="GT_WU",incoming=%u)",cnum,incoming);
  gt_assert(cg->v_spm[incoming][cnum + 1UL].ptr >
      cg->v_spm[incoming][cnum].ptr + cg->v_spm[incoming][cnum].deg);
  for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
  {
    if (edge->deleted)
      return edge;
  }
  gt_assert(false);
  return NULL;
}

void gt_contigs_graph_rm_spm_edge(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming, GtContigsGraphSpmEdge *edge)
{
  GtContigsGraphSpmEdge *reverse_edge;
  unsigned int reverse_incoming;
  gt_log_log("rm spm edge "GT_WU" -- "GT_WU"", cnum, (GtUword)edge->dest);
  for (reverse_incoming = 0; reverse_incoming < 2U; reverse_incoming++)
  {
    for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, edge->dest,
          reverse_incoming, reverse_edge))
    {
      if (reverse_edge->deleted)
        continue;
      if ((GtUword)reverse_edge->dest == cnum &&
          (reverse_edge->reverse == (incoming == reverse_incoming)))
      {
        gt_assert(cg->v_spm[incoming][cnum].deg > 0);
        (cg->v_spm[incoming][cnum].deg)--;
        edge->deleted = true;

        gt_assert(cg->v_spm[reverse_incoming][edge->dest].deg > 0);
        (cg->v_spm[reverse_incoming][edge->dest].deg)--;
        reverse_edge->deleted = true;
        return;
      }
    }
  }
  gt_assert(false);
}

void gt_contigs_graph_rm_vertex(GtContigsGraph *cg, GtUword cnum)
{
  GtContigsGraphSpmEdge *edge;
  unsigned int incoming;
  gt_log_log("rm vertex "GT_WU"", cnum);
  for (incoming = 0; incoming < 2U; incoming++)
  {
    for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
    {
      if (edge->deleted)
        continue;
      gt_contigs_graph_rm_spm_edge(cg, cnum, incoming, edge);
    }
  }
  cg->v_m[cnum].deleted = true;
}

uint64_t gt_contigs_graph_nof_optional_neighbours(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming)
{
  uint64_t n = 0;
  GtContigsGraphSpmEdge *edge;
  for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
  {
    if (edge->deleted)
      continue;
    if (GT_CONTIGS_GRAPH_IS_OPTIONAL(cg, edge->dest))
      n++;
  }
  return n;
}

void gt_contigs_graph_rm_optional_neighbours(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming)
{
  GtContigsGraphSpmEdge *edge;
  for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, cnum, incoming, edge))
  {
    if (edge->deleted)
      continue;
    if (GT_CONTIGS_GRAPH_IS_OPTIONAL(cg, edge->dest))
    {
      /*gt_log_log("rm edge "GT_WU" -- %u", cnum, edge->dest);*/
      gt_contigs_graph_rm_spm_edge(cg, cnum, incoming, edge);
      /* if optional is now untraversable, rm it */
      if (cg->v_spm[0][edge->dest].deg == 0 ||
          cg->v_spm[1][edge->dest].deg == 0)
      {
        gt_contigs_graph_rm_vertex(cg, (GtUword)edge->dest);
      }
    }
  }
}

#define GT_CONTIGS_GRAPH_SET_DIR_MARK(CG, CNUM, INCOMING)\
  if ((INCOMING) == 0) \
    (CG)->v_m[CNUM].mark0 = true; \
  else \
    (CG)->v_m[CNUM].mark1 = true

#define GT_CONTIGS_GRAPH_DIR_MARK_IS_UNSET(CG, CNUM, INCOMING)\
  (((INCOMING) == 0  && !(CG)->v_m[CNUM].mark0) || \
   ((INCOMING) == 1U && !(CG)->v_m[CNUM].mark1))

void gt_contigs_graph_simplify_from_contig(GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming, bool restrict_rm_optionals)
{
  GtUword src = cnum;
  while (true)
  {
    GtContigsGraphSpmEdge *edge;
    if (cg->v_m[cnum].visited)
      break;
    cg->v_m[cnum].visited = true;
    if (cg->v_spm[incoming][cnum].deg == (uint64_t)1)
    {
      edge = gt_contigs_graph_find_only_spm_edge(cg, cnum, incoming);
      GT_CONTIGS_GRAPH_SET_DIR_MARK(cg, cnum, incoming);
      if (GT_CONTIGS_GRAPH_DIR_MARK_IS_UNSET(cg, edge->dest, incoming))
      {
        cnum = (GtUword)edge->dest;
        incoming = GT_CONTIGS_GRAPH_VALID_PATH_DIR(incoming, edge->reverse);
        /*gt_log_log("... continue from contig "GT_WU" (incoming=%u, deg=%u)",
                   cnum, incoming, cg->v_spm[incoming][cnum].deg);*/
      }
      else
      {
        break;
      }
    }
    else
    {
      if (gt_contigs_graph_nof_optional_neighbours(cg, cnum, incoming) ==
          cg->v_spm[incoming][cnum].deg - 1U)
      {
        if (restrict_rm_optionals && gt_contigs_graph_nof_optional_neighbours(
              cg, cnum, GT_CONTIGS_GRAPH_OTHER_DIR(incoming)) != cg->v_spm[
            incoming][cnum].deg - 1U)
        {
          break;
        }
        gt_contigs_graph_rm_optional_neighbours(cg, cnum, incoming);
        if (cg->v_spm[incoming][cnum].deg != (uint64_t)1)
        {
          GT_CONTIGS_GRAPH_SET_DIR_MARK(cg, cnum, incoming);
          break;
        }
        /*gt_log_log("... repeat contig "GT_WU" (incoming=%u, deg=%u)", cnum,
            incoming, cg->v_spm[incoming][cnum].deg);*/
      }
      else
      {
        break;
      }
    }
  }
  /* reset visited mark */
  cnum = src;
  while (true)
  {
    GtContigsGraphSpmEdge *edge;
    if (!cg->v_m[cnum].visited)
      break;
    cg->v_m[cnum].visited = false;
    if (cg->v_spm[incoming][cnum].deg != (uint64_t)1)
      break;
    edge = gt_contigs_graph_find_only_spm_edge(cg, cnum, incoming);
    cnum = (GtUword)edge->dest;
    incoming = GT_CONTIGS_GRAPH_VALID_PATH_DIR(incoming, edge->reverse);
  }
}

/*
 * From each (=1)-contig.
 * - traverse in both directions
 * - if a junction is found:
 *   - if only one way is not a (>=0)-contig, rm edge to (>=0)-contigs
 *     and the (>=0)-contig as well, if possible
 *   - otherwise break traversal
 *
 */
void gt_contigs_graph_simplify(GtContigsGraph *cg, bool restrict_rm_optionals)
{
  GtUword cnum;
  gt_assert(cg != NULL);
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    cg->v_m[cnum].mark0 = false;
    cg->v_m[cnum].mark1 = false;
    cg->v_m[cnum].visited = false;
  }
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    if (!cg->v_m[cnum].deleted &&
        GT_CONTIGS_GRAPH_IS_SINGLE_COPY(cg, cnum))
    {
      if (!cg->v_m[cnum].mark0)
      {
        /*gt_log_log("simplify from contig "GT_WU" (incoming=0, deg=%u)", cnum,
              cg->v_spm[0U][cnum].deg);*/
        gt_contigs_graph_simplify_from_contig(cg, cnum, 0U,
            restrict_rm_optionals);
      }
      if (!cg->v_m[cnum].mark1)
      {
        /*gt_log_log("simplify from contig "GT_WU" (incoming=1, deg=%u)", cnum,
              cg->v_spm[1U][cnum].deg);*/
        gt_contigs_graph_simplify_from_contig(cg, cnum, 1U,
            restrict_rm_optionals);
      }
    }
  }
}

#define GT_CONTIGS_GRAPH_PATHELEM_INC 256

void gt_contigs_graph_find_path_end(GtContigsGraphPathEndInfo *info,
    GtArrayGtContigsGraphPathElem *path, GtContigsGraph *cg,
    GtUword cnum, unsigned int incoming, bool use_only_internal)
{
  bool extended = false;
  gt_assert(cg->v_spm[incoming][cnum].deg == (uint64_t)1);
  info->cnum = cnum;
  info->dir = incoming;
  info->t = GT_CONTIGS_GRAPH_JUNCTION;
  info->depth = 0;
  gt_log_log("find_path_end(cnum="GT_WU", incoming=%u)", cnum, incoming);
  path->nextfreeGtContigsGraphPathElem = 0;
  do {
    GtContigsGraphSpmEdge *edge;
    GtContigsGraphPathElem *elem;
    GT_GETNEXTFREEINARRAY(elem, path, GtContigsGraphPathElem,
      GT_CONTIGS_GRAPH_PATHELEM_INC);
    if (incoming == 1U)
      elem->dest = info->cnum;
    edge = gt_contigs_graph_find_only_spm_edge(cg, info->cnum, info->dir);
    info->depth++;
    info->dest = info->cnum;
    info->cnum = (GtUword)edge->dest;
    elem->extended = extended;
    if (incoming == 0)
    {
      elem->dest = info->cnum;
      elem->dir = info->dir;
    }
    info->dir = GT_CONTIGS_GRAPH_VALID_PATH_DIR(info->dir, edge->reverse);
    if (incoming == 1U)
      elem->dir = GT_CONTIGS_GRAPH_OTHER_DIR(info->dir);
    if (info->cnum == cnum || cg->v_m[info->cnum].visited)
    {
      info->t = GT_CONTIGS_GRAPH_CIRCULAR;
      break;
    }
    if (cg->v_spm[GT_CONTIGS_GRAPH_OTHER_DIR(info->dir)][
        info->cnum].deg != (uint64_t)1)
    {
      extended = true;
      if (use_only_internal)
        break;
    }
    cg->v_m[info->cnum].visited = true;
  } while (cg->v_spm[info->dir][info->cnum].deg == (uint64_t)1);
  if (cg->v_spm[info->dir][info->cnum].deg == 0)
  {
    if (cg->v_spm[GT_CONTIGS_GRAPH_OTHER_DIR(info->dir)][
        info->cnum].deg == (uint64_t)1)
      info->t = GT_CONTIGS_GRAPH_SINGLE_END;
    else
    {
      gt_assert(cg->v_spm[GT_CONTIGS_GRAPH_OTHER_DIR(info->dir)][
          info->cnum].deg > (uint64_t)1);
      info->t = GT_CONTIGS_GRAPH_MULTIPLE_END;
    }
  }
  info->dir = GT_CONTIGS_GRAPH_OTHER_DIR(info->dir);
}

void gt_contigs_graph_create_composite_vertex(GtContigsGraph *cg,
    GtArrayGtContigsGraphPathElem *path, GtContigsGraphPathEndInfo *info)
{
  GtUword offset, cnum, nof_units, i, new_cnum;
  unsigned int dir;
  GT_UNUSED unsigned int enddir;
  gt_log_log("... create_composite_vertex");
  for (dir = 0; dir <= 1U; dir++)
    gt_log_log("info[%u]: cnum = "GT_WU", dest = "GT_WU", depth = "GT_WU", "
        "dir = %u, t = %s",
        dir,
        info[dir].cnum,
        info[dir].dest,
        info[dir].depth,
        info[dir].dir,
        info[dir].t == GT_CONTIGS_GRAPH_CIRCULAR ? "circular" :
          (info[dir].t == GT_CONTIGS_GRAPH_SINGLE_END ? "single_end" :
            (info[dir].t == GT_CONTIGS_GRAPH_MULTIPLE_END ? "multiple_end" :
              "junction")));
  if (info[0].t == GT_CONTIGS_GRAPH_CIRCULAR)
  {
    gt_assert(info[1].t == GT_CONTIGS_GRAPH_CIRCULAR);
    gt_assert(info[0].depth == info[1].depth);
    nof_units = 1UL + info[0].depth;
    /* it would be better to break circles at an optional,
     * if one is available */
  }
  else
  {
    nof_units = 1UL + info[0].depth + info[1].depth;
  }
  new_cnum = gt_contigs_graph_append_vertex(cg, 0, 0, nof_units);
#if 0
  new_cnum = gt_contigs_graph_append_vertex(cg,
      info[0].t == GT_CONTIGS_GRAPH_JUNCTION ? 1UL : 0,
      info[1].t == GT_CONTIGS_GRAPH_JUNCTION ? 1UL : 0,
      nof_units);
#endif
  gt_log_log("appended: vertex "GT_WU" (units: "GT_WU")", new_cnum, nof_units);
  cnum = info[1].cnum;
  gt_log_log("--- start traversal from "GT_WU"", cnum);
  dir = info[1].dir;
  offset = 0;
  for (i = 0; i < nof_units; i++)
  {
    GtContigsGraphSeqUnit *unit;
    unit = cg->units + cg->v_cmp[new_cnum - cg->nof_simple_v].ptr + i;
    unit->seqnum = cnum;
    if (i == 0)
      unit->revcomp = (info[1].dir == 1U);
    else if (i == nof_units - 1U)
      unit->revcomp = (info[0].dir == 0);
    else
      unit->revcomp = (dir == 1U);
    unit->offset = offset;
    if (i + 1 < nof_units)
    {
      GtContigsGraphSpmEdge *edge;
      GtUword dest;
      if (i < info[1].depth) {
        dest = path[1].spaceGtContigsGraphPathElem[info[1].depth - 1 - i].dest;
        dir = path[1].spaceGtContigsGraphPathElem[info[1].depth - 1 - i].dir;
      }
      else { /* if (i >= info[1].depth) */
        dest = path[0].spaceGtContigsGraphPathElem[i - info[1].depth].dest;
        dir = path[0].spaceGtContigsGraphPathElem[i - info[1].depth].dir;
      }
      edge = gt_contigs_graph_find_spm_edge(cg, cnum, dir, dest);
      cnum = (GtUword)edge->dest;
      dir = GT_CONTIGS_GRAPH_VALID_PATH_DIR(dir, edge->reverse);
      offset = (GtUword)edge->ovlen;
    }
  }
  gt_log_log("cnum = "GT_WU", info[0].cnum = "GT_WU"", cnum, info[0].cnum);
  gt_assert(cnum == info[0].cnum);
  cnum = info[1].cnum;
  gt_log_log("--- start traversal from "GT_WU"", cnum);
  gt_assert(nof_units >= 1UL);
  if (info[1].t == GT_CONTIGS_GRAPH_CIRCULAR)
  {
    nof_units--;
    gt_assert(nof_units >= 1UL);
  }
  for (i = 0; i < nof_units - 1UL; i++)
  {
    GtContigsGraphSpmEdge *edge;
    GtUword src = cnum;
    GtUword dest;
    GT_UNUSED bool extended;
    if (i < info[1].depth) {
      dest = path[1].spaceGtContigsGraphPathElem[info[1].depth - 1 - i].dest;
      dir = path[1].spaceGtContigsGraphPathElem[info[1].depth - 1 - i].dir;
      extended = path[1].spaceGtContigsGraphPathElem[info[1].depth - 1 - i].
        extended;
    }
    else { /* if (i >= info[1].depth) */
      dest = path[0].spaceGtContigsGraphPathElem[i - info[1].depth].dest;
      dir = path[0].spaceGtContigsGraphPathElem[i - info[1].depth].dir;
      extended = path[0].spaceGtContigsGraphPathElem[i - info[1].depth].
        extended;
    }
    edge = gt_contigs_graph_find_spm_edge(cg, cnum, dir, dest);
    cnum = (GtUword)edge->dest;
    gt_log_log("traverse edge "GT_WU" -- "GT_WU"", src, cnum);
    /*dir = GT_CONTIGS_GRAPH_VALID_PATH_DIR(dir, edge->reverse);*/
    if (cg->v_spm[0][src].deg == (uint64_t)1 &&
        cg->v_spm[1][src].deg == (uint64_t)1)
    {
      gt_log_log("mark vertex "GT_WU" (internal)", src);
      cg->v_m[src].selected = true;
     /* if (!extended && i > 1 && i < nof_units - 2)
        cg->v_m[src].processed = true;*/
    }
  }
  if (info[1].t != GT_CONTIGS_GRAPH_CIRCULAR)
  {
    gt_log_log("cnum = "GT_WU", info[0].cnum = "GT_WU"", cnum, info[0].cnum);
    gt_assert(cnum == info[0].cnum);
  }
  else
  {
    gt_log_log("cnum = "GT_WU", info[0].dest = "GT_WU"", cnum, info[0].dest);
    gt_assert(cnum == info[0].dest);
  }
  if (info[1].t == GT_CONTIGS_GRAPH_CIRCULAR)
  {
    gt_log_log("mark vertex "GT_WU" (circle end)", info[1].cnum);
    cg->v_m[info[1].cnum].selected = true;
  }
  if (info[0].t == GT_CONTIGS_GRAPH_SINGLE_END)
  {
    gt_log_log("mark vertex "GT_WU" (single end)", info[0].cnum);
    cg->v_m[info[0].cnum].selected = true;
  }
  if (info[1].t == GT_CONTIGS_GRAPH_SINGLE_END)
  {
    gt_log_log("mark vertex "GT_WU" (single end)", info[1].cnum);
    cg->v_m[info[1].cnum].selected = true;
  }
#if 0
  for (enddir = 0; enddir <= 1U; enddir++)
  {
    if (info[dir].t == GT_CONTIGS_GRAPH_JUNCTION)
    {
      GtContigsGraphSpmEdge *edge;
      cnum = info[enddir].cnum;
      dir = info[enddir].dir;
      /* set the edge on junction, recycling the previously deleted edge */
      edge = gt_contigs_graph_find_deleted_spm_edge(cg, cnum, dir);
      edge->deleted = false;
      cg->v_spm[dir][cnum].deg++;
      edge->dest = (uint32_t)new_cnum;
      edge->reverse = (dir == enddir);
      edge->ovlen = 0; /* this should be set to the junction contig length */
      /*cg->v_m[cnum].optional = true;*/
      /* set the edge on the composite contig */
      edge = cg->e_spm[enddir] + cg->v_spm[enddir][new_cnum].ptr;
      gt_assert(edge->deleted == false);
      edge->dest = (uint32_t)cnum;
      edge->reverse = (dir == enddir);
      edge->ovlen = 0; /* this should be set to the junction contig length */
    }
  }
#endif
}

void gt_contigs_graph_extend_contigs(GtContigsGraph *cg, bool use_only_internal)
{
  GtUword cnum, nof_v_before;
  GtArrayGtContigsGraphPathElem path[2];
  GT_INITARRAY(&(path[0]), GtContigsGraphPathElem);
  GT_INITARRAY(&(path[1]), GtContigsGraphPathElem);
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    cg->v_m[cnum].selected = false;
    cg->v_m[cnum].processed = false;
    cg->v_m[cnum].visited = false;
  }
  nof_v_before = cg->nof_v;
  gt_log_log("nof_v before extending contigs = "GT_WU"", cg->nof_v);
  for (cnum = 0; cnum < nof_v_before; cnum++)
  {
    GtContigsGraphPathEndInfo info[2];
    GtUword deg[2];
    unsigned int startdir;
    deg[0] = (GtUword)cg->v_spm[0][cnum].deg;
    deg[1] = (GtUword)cg->v_spm[1][cnum].deg;
    if ((cg->v_m[cnum].deleted) ||
        (cg->v_m[cnum].optional) ||
        (deg[0] > 1UL) || /* junction */
        (deg[1] > 1UL) || /* junction */
        (deg[0] == 0 && deg[1] == 0)) /* isolated */
      continue;
    gt_log_log("extend_contigs, cnum = "GT_WU"", cnum);
    for (startdir = 0; startdir < 2U; startdir++)
    {
      if (cg->v_spm[startdir][cnum].deg == 0)
      {
        info[startdir].depth = 0;
        info[startdir].t = GT_CONTIGS_GRAPH_SINGLE_END;
        info[startdir].dest = (GtUword)
          gt_contigs_graph_find_only_spm_edge(cg, cnum,
              GT_CONTIGS_GRAPH_OTHER_DIR(startdir))->dest;
        info[startdir].cnum = cnum;
        info[startdir].dir = GT_CONTIGS_GRAPH_OTHER_DIR(startdir);
        continue;
      }
      else
      {
        {
          GtUword cnum2;
          for (cnum2 = 0; cnum2 < nof_v_before; cnum2++)
          {
            cg->v_m[cnum2].visited = false;
          }
        }
        gt_contigs_graph_find_path_end(&(info[startdir]), &(path[startdir]),
            cg, cnum, startdir, use_only_internal);
      }
    }
    gt_contigs_graph_create_composite_vertex(cg, path, info);
    {
      GtUword cnum2;
      for (cnum2 = 0; cnum2 < nof_v_before; cnum2++)
      {
        if (cg->v_m[cnum2].processed && !cg->v_m[cnum2].deleted
        && info[0].t != GT_CONTIGS_GRAPH_CIRCULAR)
          gt_contigs_graph_rm_vertex(cg, cnum2);
      }
    }
  }
  for (cnum = 0; cnum < nof_v_before; cnum++)
  {
    if (cg->v_m[cnum].selected && !cg->v_m[cnum].deleted)
      gt_contigs_graph_rm_vertex(cg, cnum);
  }
  GT_FREEARRAY(&path[0], GtContigsGraphPathElem);
  GT_FREEARRAY(&path[1], GtContigsGraphPathElem);
}

int gt_contigs_graph_show_dot_subgraph(GtContigsGraph *cg,
    GtFile *outfp, GtUword *cnums, GtUword nofcnums,
    GtUword maxdepth, GtError *err)
{
  GtArray *stack;
  struct { GtUword v, d; } to_add, *current;
  GtUword next, cnum, i, v, d;
  GtContigsGraphSpmEdge *edge;
  unsigned int dir;

  stack = gt_array_new(sizeof (to_add));
  for (cnum = 0; cnum < cg->nof_v; cnum++)
  {
    cg->v_m[cnum].selected = false;
    cg->v_m[cnum].processed = false;
  }
  for (i = 0; i < nofcnums; i++)
  {
    cnum = cnums[i];
    if (cnum >= cg->nof_v)
    {
      gt_error_set(err, "Error by context output for contig "GT_WU": "
          "number of contigs is "GT_WU"", cnum, cg->nof_v);
      gt_array_delete(stack);
      return -1;
    }
    to_add.v = cnum;
    cg->v_m[cnum].selected = true;
    to_add.d = 1UL;
    gt_array_add(stack, to_add);
  }
  gt_file_xprintf(outfp, GT_CONTIGS_GRAPH_DOT_HEADER);
  for (next = 0; next < gt_array_size(stack); next++)
  {
    current = gt_array_get(stack, next);
    v = current->v;
    d = current->d;
    gt_contigs_graph_show_dot_for_contig(cg, v, outfp);
    cg->v_m[v].processed = true;
    to_add.d = d + 1UL;
    for (dir = 0; dir < 2U; dir++)
    {
      for (GT_CONTIGS_GRAPH_EACH_SPM_EDGE(cg, v, dir, edge))
      {
        if (cg->dot_show_deleted || !edge->deleted)
        {
          if (!cg->v_m[edge->dest].selected)
          {
            if (d < maxdepth)
            {
              to_add.v = (GtUword)edge->dest;
              gt_array_add(stack, to_add);
              cg->v_m[edge->dest].selected = true;
            }
            else
            {
              gt_file_xprintf(outfp, "  %u [style=invisible];\n", edge->dest);
            }
          }
        }
      }
    }
  }
  gt_array_delete(stack);
  gt_file_xprintf(outfp, GT_CONTIGS_GRAPH_DOT_FOOTER);
  return 0;
}

void gt_contigs_graph_delete(GtContigsGraph *cg)
{
  if (cg != NULL)
  {
    gt_free(cg->v_spm[0]);
    gt_free(cg->v_spm[1]);
    gt_free(cg->v_scf);
    gt_free(cg->v_d);
    gt_free(cg->v_m);
    gt_free(cg->e_spm[0]);
    gt_free(cg->e_spm[1]);
    gt_free(cg->e_scf);
    gt_free(cg->v_cmp);
    gt_free(cg->units);
    gt_reads_libraries_table_delete(cg->rlt);
    gt_free(cg);
  }
}
