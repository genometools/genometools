/*
  Copyright (c) 2004-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/translator_api.h"
#include "core/trans_table.h"
#include "core/unused_api.h"
#include "gth/gthstopcodon.h"
#include "gth/indent.h"
#include "gth/gthorf.h"
#include "gth/gthtrans.h"

#define TRANSLATIONLINEWIDTH    60
#define EXONSEPARATORSTRING     " : "
#define DOTLINECHAR             '.'
#define DOTSPACING              10

#define INSERT_EXON_BORDER(LINE)\
        if (gth_spliced_seq_pos_is_border(splicedseq, incounter))\
        {\
          for (i = 0; i < exonseparatorwidth; i++)\
          {\
            (LINE)[outcounter++] = (unsigned char) EXONSEPARATORSTRING[i];\
          }\
        }\
        incounter++;

#define INSERT_NEWLINE(LINE)\
        if (incounter % TRANSLATIONLINEWIDTH == 0)\
        {\
          (LINE)[outcounter++] = (unsigned char) '\n';\
        }

#define REPLACE_STOPCODON(LINE)\
        if (gs2out && (LINE)[outcounter] == GT_STOP_AMINO)\
        {\
          (LINE)[outcounter] = (unsigned char) GS2_STOPCODON;\
        }

#define ADD_CODONCHAR(FRAMEOUT, FRAMEIN, FRAMENUM)\
        if ((incounter > (FRAMENUM)) &&\
           (incounter < splicedseq->splicedseqlen - 1) &&\
           ((incounter - 1 - (FRAMENUM)) % GT_CODON_LENGTH == 0))\
        {\
          (FRAMEOUT)[outcounter] = (FRAMEIN)[(incounter - 1 - (FRAMENUM))\
                                             / GT_CODON_LENGTH];\
          REPLACE_STOPCODON(FRAMEOUT);\
        }\
        else\
        {\
          (FRAMEOUT)[outcounter] = (unsigned char) ' ';\
        }\
        outcounter++;

#define OUTCOUNTERCHECK\
        gt_assert(outcounter == outlen);

#define OUTPUT_EXONPOSSPACE\
        gt_file_xprintf(out->outfp, "%*s  ", out->widthforgenpos, "")

#define OUTPUT_LINE(LINE)\
        i = outcounter;\
        while (i < outlen)\
        {\
          if ((LINE)[i] == '\n')\
          {\
            break;\
          }\
          else\
          {\
            gt_file_xfputc((char) (LINE)[i++], out->outfp);\
          }\
        }\
        gt_file_xfputc('\n', out->outfp);

static void createoutputlines(char *dotline,
                              char *template_out,
                              char *frame0_out,
                              char *frame1_out,
                              char *frame2_out,
                              char *template_in,
                              char *frame0_in,
                              char *frame1_in,
                              char *frame2_in,
                              GthSplicedSeq *splicedseq,
                              unsigned long exonseparatorwidth,
                              GT_UNUSED unsigned long outlen, bool gs2out)
{
  unsigned long i,
       dotcounter = DOTSPACING,
       incounter  = 0,
       outcounter = 0;

  /* fill dot line */
  while (incounter < splicedseq->splicedseqlen) {
    if (dotcounter == DOTSPACING) {
      dotline[outcounter++] = (unsigned char) DOTLINECHAR;
      dotcounter = 1;
    }
    else {
      dotline[outcounter++] = (unsigned char) ' ';
      dotcounter++;
    }

    INSERT_EXON_BORDER(dotline);
    INSERT_NEWLINE(dotline);
  }
  OUTCOUNTERCHECK;

  /* fill template line */
  incounter  = 0;
  outcounter = 0;
  while (incounter < splicedseq->splicedseqlen) {
    template_out[outcounter++] = template_in[incounter];

    INSERT_EXON_BORDER(template_out);
    INSERT_NEWLINE(template_out);
  }
  OUTCOUNTERCHECK;

  /* fill frame0 line */
  incounter  = 0;
  outcounter = 0;
  while (incounter < splicedseq->splicedseqlen) {
    ADD_CODONCHAR(frame0_out, frame0_in, 0);

    INSERT_EXON_BORDER(frame0_out);
    INSERT_NEWLINE(frame0_out);
  }
  OUTCOUNTERCHECK;

  /* fill frame1 line */
  incounter  = 0;
  outcounter = 0;
  while (incounter < splicedseq->splicedseqlen) {
    ADD_CODONCHAR(frame1_out, frame1_in, 1);

    INSERT_EXON_BORDER(frame1_out);
    INSERT_NEWLINE(frame1_out);
  }

  /* fill frame2 line */
  incounter  = 0;
  outcounter = 0;
  while (incounter < splicedseq->splicedseqlen) {
    ADD_CODONCHAR(frame2_out, frame2_in, 2);

    INSERT_EXON_BORDER(frame2_out);
    INSERT_NEWLINE(frame2_out);
  }
}

static void showoutputlines(char *dotline,
                            char *template_out,
                            char *frame0_out,
                            char *frame1_out,
                            char *frame2_out,
                            unsigned long outlen,
                            bool gen_strand_forward,
                            unsigned long gen_total_length,
                            unsigned long gen_offset,
                            unsigned long *positionmapping,
                            GthOutput *out)
{
  unsigned long i, outcounter  = 0, origcounter = 0;

  while (outcounter < outlen) {
    OUTPUT_EXONPOSSPACE;
    OUTPUT_LINE(dotline);

    /* output exon position */
    gt_file_xprintf(out->outfp, "%*lu  ", out->widthforgenpos,
                    SHOWGENPOS(gen_strand_forward, gen_total_length,
                               gen_offset, positionmapping[origcounter]));
    OUTPUT_LINE(template_out);

    OUTPUT_EXONPOSSPACE;
    OUTPUT_LINE(frame0_out);

    OUTPUT_EXONPOSSPACE;
    OUTPUT_LINE(frame1_out);

    OUTPUT_EXONPOSSPACE;
    OUTPUT_LINE(frame2_out);

    (void) gt_file_xprintf(out->outfp, "\n\n");
    outcounter = i + 1;
    origcounter += TRANSLATIONLINEWIDTH;
  }
}

static void showtranslation(GthSplicedSeq *splicedseq,
                            char *frame0_in,
                            char *frame1_in,
                            char *frame2_in,
                            GtArray *exons,
                            bool gen_strand_forward,
                            unsigned long gen_total_length,
                            unsigned long gen_offset,
                            unsigned int indentlevel,
                            GthOutput *out)
{
  char *dotline, *template_out, *frame0_out, *frame1_out, *frame2_out;
  unsigned long i, exonseparatorwidth =  strlen(EXONSEPARATORSTRING),
                outlen = splicedseq->splicedseqlen +
                         ((gt_array_size(exons) - 1) * exonseparatorwidth) +
                         (splicedseq->splicedseqlen / TRANSLATIONLINEWIDTH);
  GtFile *outfp = out->outfp;

  dotline      = gt_malloc(sizeof (unsigned char) * outlen);
  template_out = gt_malloc(sizeof (unsigned char) * outlen);
  frame0_out   = gt_malloc(sizeof (unsigned char) * outlen);
  frame1_out   = gt_malloc(sizeof (unsigned char) * outlen);
  frame2_out   = gt_malloc(sizeof (unsigned char) * outlen);

  createoutputlines(dotline, template_out, frame0_out, frame1_out, frame2_out,
                    (char*) splicedseq->splicedseq, frame0_in, frame1_in,
                    frame2_in, splicedseq, exonseparatorwidth, outlen,
                    out->gs2out);

  if (out->xmlout) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<translation>\n");
    indentlevel++;

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<gDNA_template>");
    for (i = 0; i < outlen; i++) {
      if (template_out[i] != '\n') {
        gt_file_xfputc(template_out[i], outfp);
      }
    }
    gt_file_xprintf(outfp, "</gDNA_template>\n");

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<first_frame>");
    for (i = 0; i < outlen; i++) {
      if (frame0_out[i] != '\n') {
        gt_file_xfputc(frame0_out[i], outfp);
      }
    }
    gt_file_xprintf(outfp, "</first_frame>\n");

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<second_frame>");
    for (i = 0; i < outlen; i++) {
      if (frame1_out[i] != '\n') {
        gt_file_xfputc(frame1_out[i], outfp);
      }
    }
    gt_file_xprintf(outfp, "</second_frame>\n");

    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<third_frame>");
    for (i = 0; i < outlen; i++) {
      if (frame2_out[i] != '\n') {
        gt_file_xfputc(frame2_out[i], outfp);
      }
    }
    gt_file_xprintf(outfp, "</third_frame>\n");

    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</translation>\n");
  }
  else {
    showoutputlines(dotline, template_out, frame0_out, frame1_out, frame2_out,
                    outlen, gen_strand_forward, gen_total_length,
                    gen_offset, splicedseq->positionmapping, out);
  }

  gt_free(dotline);
  gt_free(template_out);
  gt_free(frame0_out);
  gt_free(frame1_out);
  gt_free(frame2_out);
}

void gt_outputtranslationandorf(unsigned long pglnum, const GthAGS *ags,
                                unsigned long agsnum,
                                unsigned long translationtable,
                                GthInput *input,
                                unsigned int indentlevel,
                                GthOutput *out)
{
  unsigned long i;
  unsigned int nframe;
  const unsigned char *gen_seq_orig;
  GtStr *frame[3];
  char translated;
  GtTranslatorStatus status;
  GtTranslator *translator;
  GtTransTable *transtable;
  GtCodonIterator *ci;
  GthSplicedSeq *spliced_seq;
  GtArray *ranges;
  GtFile *outfp = out->outfp;

  /* output header */
  if (out->xmlout) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<three_phase_translation "
                    "xmlns=\"http://www.genomethreader.org/GTH_output/"
                    "PGL_module/predicted_gene_location/AGS_information/"
                    "three_phase_translation/\">\n");
    indentlevel++;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<description PGL_serial=\"%lu\" "
                              "AGS_serial=\"%lu\" gDNA_strand=\"%c\"/>\n",
                       pglnum + OUTPUTOFFSET, agsnum + OUTPUTOFFSET,
                       SHOWSTRAND(gth_ags_is_forward(ags)));
  }
  else {
    gt_file_xprintf(outfp, "3-phase translation of AGS-%lu (%cstrand):\n\n",
                       agsnum + OUTPUTOFFSET,
                       SHOWSTRAND(gth_ags_is_forward(ags)));
  }

  ranges = gt_array_new(sizeof (GtRange));
  for (i = 0; i < gt_array_size(ags->exons); i++)
    gt_array_add(ranges, ((GthExonAGS*) gt_array_get(ags->exons, i))->range);

  /* get genomic sequence */
  gen_seq_orig = gth_input_original_genomic_sequence(input,
                                                     gth_ags_filenum(ags),
                                                     gth_ags_is_forward(ags));

  spliced_seq = gth_spliced_seq_new(gen_seq_orig, ranges);

  frame[0] = gt_str_new();
  frame[1] = gt_str_new();
  frame[2] = gt_str_new();

  /* prepare for translation */
  ci = gt_codon_iterator_simple_new((const char*) spliced_seq->splicedseq,
                                    spliced_seq->splicedseqlen, NULL);
  gt_assert(ci);
  transtable = gt_trans_table_new(translationtable, NULL);
  gt_assert(transtable);

  /* translate the template in all three frames */
  translator = gt_translator_new_with_table(transtable, ci);
  status = gt_translator_next(translator, &translated, &nframe, NULL);
  while (status == GT_TRANSLATOR_OK) {
    gt_str_append_char(frame[nframe], translated);
    status = gt_translator_next(translator, &translated, &nframe, NULL);
  }
  gt_assert(status != GT_TRANSLATOR_ERROR);
  gt_translator_delete(translator);
  gt_trans_table_delete(transtable);
  gt_codon_iterator_delete(ci);

  /* show the translation */
  showtranslation(spliced_seq, gt_str_get(frame[0]), gt_str_get(frame[1]),
                  gt_str_get(frame[2]), ags->exons, gth_ags_is_forward(ags),
                  gth_ags_total_length(ags), gth_ags_genomic_offset(ags),
                  indentlevel, out);

  /* show the (consolidated) ORFs */
  gthshowORFs(gt_str_get(frame[0]), gt_str_get(frame[1]), gt_str_get(frame[2]),
              gt_str_length(frame[0]), gt_str_length(frame[1]),
              gt_str_length(frame[2]), gth_ags_is_forward(ags),
              gth_ags_total_length(ags), gth_ags_genomic_offset(ags),
              gt_str_get(ags->gen_id), pglnum, agsnum, spliced_seq,
              indentlevel, out);

  if (out->xmlout) {
    indentlevel--;
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "</three_phase_translation>\n");
  }

  gth_spliced_seq_delete(spliced_seq);
  gt_array_delete(ranges);
  gt_str_delete(frame[0]);
  gt_str_delete(frame[1]);
  gt_str_delete(frame[2]);
}
