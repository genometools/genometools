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

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "libgtcore/cstr.h"
#include "libgtcore/fasta.h"
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtltr/ltrfileout_stream.h"
#include "libgtltr/ltr_visitor.h"

#define FSWIDTH 60

struct LTRFileOutStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  Bioseq *bioseq;
  const char *fileprefix;
  GenFile *metadata_file,
          *tabout_file,
          *pbsout_file,
          *pptout_file,
          *ltr5out_file,
          *ltr3out_file;
  Hashtable *pdomout_files;
  LTRVisitor *lv;
  int tests_to_run;
  LTRElement element;
};

#define ltr_fileout_stream_cast(GS)\
        genome_stream_cast(ltr_fileout_stream_class(), GS)

int ltr_fileout_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                               Error *e)
{
  LTRFileOutStream *ls;
  Range lltr_rng, rltr_rng, rng, ppt_rng, pbs_rng;
  char *ppt_seq, *pbs_seq;
  Str *pdoms;
  int had_err;

  error_check(e);
  ls = ltr_fileout_stream_cast(gs);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (LTRElement));
  ls->element.pdoms = array_new(sizeof (GenomeFeature*));

  /* get annotations from parser */
  had_err = genome_stream_next_tree(ls->in_stream, gn, e);
  if (!had_err && *gn)
  {
    GenomeNodeIterator* gni;
    GenomeNode *mygn;
    /* fill LTRElement structure from GFF3 subgraph */
    gni = genome_node_iterator_new(*gn);
    for (mygn = *gn; mygn; mygn = genome_node_iterator_next(gni))
      genome_node_accept(mygn, (GenomeVisitor*) ls->lv, e);
    genome_node_iterator_delete(gni);
  }

  if (ls->element.mainnode)
  {
    unsigned long seqid, i;
    char *outseq,
         desc[BUFSIZ];
    Range ltr3_rng, ltr5_rng;

    /* find sequence in Bioseq multifasta */
    const char *sreg = str_get(genome_node_get_seqid((GenomeNode*)
                                                     ls->element.mainnode));
    sscanf(sreg,"seq%lu", &seqid);
    Seq *seq = bioseq_get_seq(ls->bioseq, seqid);
    ls->element.seqnr = seqid;

    ltrelement_format_description(&ls->element, desc, BUFSIZ-1);

    /* output basic retrotransposon data */
    lltr_rng = genome_node_get_range((GenomeNode*) ls->element.leftLTR);
    rltr_rng = genome_node_get_range((GenomeNode*) ls->element.rightLTR);
    rng = genome_node_get_range((GenomeNode*) ls->element.mainnode);
    genfile_xprintf(ls->tabout_file,
                    "%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
                    rng.start,
                    rng.end,
                    ltrelement_length(&ls->element),
                    lltr_rng.start,
                    lltr_rng.end,
                    ltrelement_leftltrlen(&ls->element),
                    rltr_rng.start,
                    rltr_rng.end,
                    ltrelement_rightltrlen(&ls->element));

    /* output PPT */
    if (ls->element.ppt)
    {
      Strand ppt_strand;
      ppt_strand = genome_feature_get_strand(ls->element.ppt);
      ppt_rng = genome_node_get_range((GenomeNode*) ls->element.ppt);
      ppt_seq = ltrelement_get_sequence(ppt_rng.start, ppt_rng.end,
                                        ppt_strand,
                                        seq, e);
      fasta_show_entry_generic(desc, ppt_seq, range_length(ppt_rng), FSWIDTH,
                               ls->pptout_file);
      genfile_xprintf(ls->tabout_file,
                      "%lu\t%lu\t%s\t%c\t%d\t",
                      ppt_rng.start,
                      ppt_rng.end,
                      ppt_seq,
                      STRANDCHARS[ppt_strand],
                      (ppt_strand == STRAND_FORWARD ?
                          abs(rltr_rng.start - ppt_rng.end) :
                          abs(lltr_rng.end - ppt_rng.start)));
      ma_free((char*) ppt_seq);
    } else genfile_xprintf(ls->tabout_file, "\t\t\t\t\t");

    /* output PBS */
    if (ls->element.pbs)
    {
      Strand pbs_strand;

      pbs_strand = genome_feature_get_strand(ls->element.pbs);
      pbs_rng = genome_node_get_range((GenomeNode*) ls->element.pbs);
      pbs_seq = ltrelement_get_sequence(pbs_rng.start, pbs_rng.end,
                                        pbs_strand,
                                        seq, e);
      fasta_show_entry_generic(desc, pbs_seq, range_length(pbs_rng), FSWIDTH,
                               ls->pbsout_file);
      genfile_xprintf(ls->tabout_file,
                      "%lu\t%lu\t%c\t%s\t%s\t%s\t%s\t%s\t",
                      pbs_rng.start,
                      pbs_rng.end,
                      STRANDCHARS[pbs_strand],
                      genome_feature_get_attribute((GenomeNode*)
                                              ls->element.pbs, "trna"),
                      pbs_seq,
                      genome_feature_get_attribute((GenomeNode*)
                                              ls->element.pbs, "trnaoffset"),
                      genome_feature_get_attribute((GenomeNode*)
                                              ls->element.pbs, "pbsoffset"),
                      genome_feature_get_attribute((GenomeNode*)
                                              ls->element.pbs, "edist"));
      ma_free((char*) pbs_seq);
    } else genfile_xprintf(ls->tabout_file, "\t\t\t\t\t\t\t\t");

    /* output protein domains */
    if (array_size(ls->element.pdoms) > 0)
    {

      pdoms = str_new();
      ltrelement_format_description(&ls->element, desc, BUFSIZ-1);
      if (STRAND_REVERSE == genome_feature_get_strand(ls->element.mainnode))
        array_reverse(ls->element.pdoms);
      for (i=0;i<array_size(ls->element.pdoms);i++)
      {
        GenFile *genfile = NULL;
        Range pdom_rng;
        GenomeFeature *gf = *(GenomeFeature**) array_get(ls->element.pdoms, i);
        char *pdom_seq;
        const char *pfamname = genome_feature_get_attribute((GenomeNode*)gf,
                                                            "pfamname");
        str_append_cstr(pdoms, pfamname);
        str_append_cstr(pdoms,"(");
        str_append_ulong(pdoms, (unsigned long) genome_feature_get_phase(gf));
        str_append_char(pdoms, STRANDCHARS[genome_feature_get_strand(gf)]);
        str_append_cstr(pdoms,")");
        if (i != array_size(ls->element.pdoms)-1)
          str_append_cstr(pdoms, "/");

        /* get protein domain output file name */
        genfile = (GenFile*) hashtable_get(ls->pdomout_files, pfamname);
        if (genfile == NULL)
        {
          /* no file opened for this domain yet, do it */
          char buffer[BUFSIZ];
          snprintf(buffer, BUFSIZ-1, "%s_pdom_%s.fas",
                   ls->fileprefix, pfamname);
          genfile = genfile_open(GFM_UNCOMPRESSED, buffer, "w+");
          hashtable_add(ls->pdomout_files, cstr_dup(pfamname), genfile);
        }
        pdom_rng = genome_node_get_range((GenomeNode*) gf);
        pdom_seq = ltrelement_get_sequence(pdom_rng.start,
                                           pdom_rng.end,
                                           genome_feature_get_strand(gf),
                                           seq, e);
        /* write out sequence to file */
        fasta_show_entry_generic(desc,
                                 pdom_seq,
                                 range_length(pdom_rng),
                                 FSWIDTH,
                                 genfile);
        ma_free(pdom_seq);
      }
      /* write domains to tabular outfile */
      genfile_xprintf(ls->tabout_file, "%s", str_get(pdoms));
      str_delete(pdoms);
    }

    /* output LTRs (we just expect them to exist) */
    switch (genome_feature_get_strand(ls->element.mainnode))
    {
      case STRAND_REVERSE:
        ltr5_rng = rltr_rng;
        ltr3_rng = lltr_rng;
        break;
      case STRAND_FORWARD:
      default:
        ltr5_rng = lltr_rng;
        ltr3_rng = rltr_rng;
        break;
    }
    outseq = ltrelement_get_sequence(ltr5_rng.start,
                              ltr5_rng.end,
                              genome_feature_get_strand(ls->element.mainnode),
                              seq, e);
    fasta_show_entry_generic(desc,
                             outseq,
                             range_length(ltr5_rng),
                             FSWIDTH,
                             ls->ltr5out_file);
    ma_free(outseq);
    outseq = ltrelement_get_sequence(ltr3_rng.start,
                              ltr3_rng.end,
                              genome_feature_get_strand(ls->element.mainnode),
                              seq, e);
    fasta_show_entry_generic(desc,
                             outseq,
                             range_length(ltr3_rng),
                             FSWIDTH,
                             ls->ltr3out_file);
    ma_free(outseq);

  }
  array_delete(ls->element.pdoms);
  genfile_xprintf(ls->tabout_file, "\n");
  return had_err;
}

static void write_metadata(GenFile *metadata_file,
                           int tests_to_run,
                           PPTOptions *ppt_opts,
                           PBSOptions *pbs_opts,
                           PdomOptions *pdom_opts,
                           const char *trnafilename,
                           const char *seqfilename,
                           const char *gfffilename)
{
  unsigned long i;
  char buffer[PATH_MAX+1];
  bool has_cwd;

  /* get working directory */
  memset(buffer, 0, PATH_MAX+1);
  has_cwd = (getcwd(buffer, PATH_MAX) != NULL);

  /* append working dir to relative paths if necessary */
  if (seqfilename[0] != '/' && has_cwd)
    genfile_xprintf(metadata_file,
                    "Sequence file used\t%s/%s\n", buffer, seqfilename);
  else
    genfile_xprintf(metadata_file,
                    "Sequence file used\t%s\n", seqfilename);
  if (gfffilename[0] != '/' && has_cwd)
    genfile_xprintf(metadata_file,
                    "GFF3 input used\t%s/%s\n", buffer, gfffilename);
  else
    genfile_xprintf(metadata_file,
                    "GFF3 input used\t%s\n", gfffilename);

  if (tests_to_run & LTRDIGEST_RUN_PPT)
  {
    genfile_xprintf(metadata_file,
                    "PPT minimum length\t%u\t6\n", ppt_opts->ppt_minlen);
    genfile_xprintf(metadata_file,
                    "U-box minimum length\t%u\t3\n", ppt_opts->ubox_minlen);
    genfile_xprintf(metadata_file,
                    "PPT search radius\t%u\t30\n", ppt_opts->radius);
  }

  if (tests_to_run & LTRDIGEST_RUN_PBS)
  {
    if (trnafilename[0] != '/' && has_cwd)
      genfile_xprintf(metadata_file,
                      "tRNA library for PBS detection\t%s/%s\n",
                      buffer, trnafilename);
    else
      genfile_xprintf(metadata_file,
                       "tRNA library for PBS detection\t%s\n",
                       trnafilename);
    genfile_xprintf(metadata_file,
                    "PBS/tRNA minimum alignment length\t%u\t11\n",
                    pbs_opts->ali_min_len);
    genfile_xprintf(metadata_file,
                    "PBS/tRNA maximum unit edit distance\t%u\t1\n",
                    pbs_opts->max_edist);
    genfile_xprintf(metadata_file,
                    "PBS max offset from 5' LTR end\t%u\t5\n",
                    pbs_opts->max_offset);
    genfile_xprintf(metadata_file,
                    "tRNA max offset from 3' end\t%u\t10\n",
                    pbs_opts->max_offset_trna);
    genfile_xprintf(metadata_file,
                    "PBS search radius\t%d\t30\n", pbs_opts->radius);
  }

  if (tests_to_run & LTRDIGEST_RUN_PDOM)
  {
    genfile_xprintf(metadata_file,
                    "Protein domains\t%lu (",
                    array_size(pdom_opts->plan7_ts));
    for (i=0;i<array_size(pdom_opts->plan7_ts);i++)
    {
      struct plan7_s *model = *(struct plan7_s **)
                                      array_get(pdom_opts->plan7_ts, i);
      genfile_xprintf(metadata_file, "%s", model->name);
      if (i != array_size(pdom_opts->plan7_ts)-1)
        genfile_xprintf(metadata_file, ", ");
    }
    genfile_xprintf(metadata_file, ")\n");
    genfile_xprintf(metadata_file,
                    "pHMM e-value cutoff \t%g\t%g\n",
                    pdom_opts->evalue_cutoff,
                    0.000001);
  }
  genfile_xprintf(metadata_file, "\n");
}

void ltr_fileout_stream_free(GenomeStream *gs)
{
  LTRFileOutStream *ls = ltr_fileout_stream_cast(gs);
  if (ls->tabout_file)
    genfile_close(ls->tabout_file);
  if (ls->metadata_file)
    genfile_close(ls->metadata_file);
  if (ls->pbsout_file)
    genfile_close(ls->pbsout_file);
  if (ls->pptout_file)
    genfile_close(ls->pptout_file);
  if (ls->ltr5out_file)
    genfile_close(ls->ltr5out_file);
  if (ls->ltr3out_file)
    genfile_close(ls->ltr3out_file);
  hashtable_delete(ls->pdomout_files);
  genome_visitor_delete((GenomeVisitor*) ls->lv);
  genome_stream_delete(ls->in_stream);
}

const GenomeStreamClass* ltr_fileout_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (LTRFileOutStream),
                                         ltr_fileout_stream_next_tree,
                                         ltr_fileout_stream_free };
  return &gsc;
}

GenomeStream* ltr_fileout_stream_new(GenomeStream *in_stream,
                                     int tests_to_run,
                                     Bioseq *bioseq,
                                     char *file_prefix,
                                     PPTOptions *ppt_opts,
                                     PBSOptions *pbs_opts,
                                     PdomOptions *pdom_opts,
                                     const char *trnafilename,
                                     const char *seqfilename,
                                     const char *gfffilename)
{
  GenomeStream *gs;
  LTRFileOutStream *ls;
  char fn[BUFSIZ];

  assert(file_prefix && in_stream && bioseq && ppt_opts
          && pbs_opts && pdom_opts);

  gs = genome_stream_create(ltr_fileout_stream_class(), true);
  ls = ltr_fileout_stream_cast(gs);

  /* ref GFF input stream and sequences*/
  ls->in_stream = genome_stream_ref(in_stream);
  ls->bioseq = bioseq;
  ls->tests_to_run = tests_to_run;

  /* open outfiles */
  ls->fileprefix = file_prefix;
  snprintf(fn, BUFSIZ-1, "%s_tabout.csv", file_prefix);
  ls->tabout_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");
  if (tests_to_run & LTRDIGEST_RUN_PPT)
  {
    snprintf(fn, BUFSIZ-1, "%s_ppt.fas", file_prefix);
    ls->pptout_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");
  }
  if (tests_to_run & LTRDIGEST_RUN_PBS)
  {
    snprintf(fn, BUFSIZ-1, "%s_pbs.fas", file_prefix);
    ls->pbsout_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");
  }
  snprintf(fn, BUFSIZ-1, "%s_5ltr.fas", file_prefix);
  ls->ltr5out_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");
  snprintf(fn, BUFSIZ-1, "%s_3ltr.fas", file_prefix);
  ls->ltr3out_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");
  snprintf(fn, BUFSIZ-1, "%s_conditions.csv", file_prefix);
  ls->metadata_file = genfile_open(GFM_UNCOMPRESSED, fn, "w+");

  /* create hashtable to hold protein domain output files */
  ls->pdomout_files = hashtable_new(HASH_STRING,
                                    ma_free_func,
                                    (FreeFunc) genfile_close);

  /* log run conditions in file */
  write_metadata(ls->metadata_file,
                 tests_to_run,
                 ppt_opts,
                 pbs_opts,
                 pdom_opts,
                 trnafilename,
                 seqfilename,
                 gfffilename);

  /* print tabular outfile headline */
  genfile_xprintf(ls->tabout_file,
              "element start\telement end\telement length\t"
              "lLTR start\tlLTR end\tlLTR length\t"
              "rLTR start\trLTR end\trLTR length\t"
              "PPT start\tPPT end\tPPT motif\tPPT strand\tPPT offset\t");
  if (tests_to_run & LTRDIGEST_RUN_PBS)
    genfile_xprintf(ls->tabout_file,
              "PBS start\tPBS end\tPBS strand\ttRNA\tRNA motif\tPBS offset\t"
              "tRNA offset\tPBS/tRNA edist\t");
  if (tests_to_run & LTRDIGEST_RUN_PDOM)
    genfile_xprintf(ls->tabout_file, "Protein domain hits\n");

  /* create visitor */
  ls->lv = (LTRVisitor*) ltr_visitor_new(&ls->element);
  return gs;
}
