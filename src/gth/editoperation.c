/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include "gth/indent.h"
#include "gth/editoperation.h"

Eoptype gt_editoperation_type(Editoperation eop, bool proteineop)
{
  Editoperation maxlen = proteineop ? MAXIDENTICALLENGTH_PROTEIN
                                    : MAXIDENTICALLENGTH;

  if (eop & maxlen) {
    switch (eop & ~maxlen) {
      case 0:
        return EOP_TYPE_MATCH;
      case DELETIONEOP:
        return EOP_TYPE_INTRON;
      case DELETION_WITH_1_GAP_EOP:
        gt_assert(proteineop);
        return EOP_TYPE_INTRON_WITH_1_BASE_LEFT;
      case DELETION_WITH_2_GAPS_EOP:
        gt_assert(proteineop);
        return EOP_TYPE_INTRON_WITH_2_BASES_LEFT;
      default: gt_assert(0); return EOP_TYPE_MATCH;
    }
  }
  else {
    switch (eop) {
      case MISMATCHEOP:
        return EOP_TYPE_MISMATCH;
      case DELETIONEOP:
        return EOP_TYPE_DELETION;
      case INSERTIONEOP:
        return EOP_TYPE_INSERTION;
      case MISMATCH_WITH_1_GAP_EOP:
        return EOP_TYPE_MISMATCH_WITH_1_GAP;
      case MISMATCH_WITH_2_GAPS_EOP:
        gt_assert(proteineop);
        return EOP_TYPE_MISMATCH_WITH_2_GAPS;
      case DELETION_WITH_1_GAP_EOP:
        gt_assert(proteineop);
        return EOP_TYPE_DELETION_WITH_1_GAP;
      case DELETION_WITH_2_GAPS_EOP:
        gt_assert(proteineop);
        return EOP_TYPE_DELETION_WITH_2_GAPS;
      default: gt_assert(0); return EOP_TYPE_MISMATCH; /* illegal edit op. */
    }
  }
}

unsigned int gt_editoperation_length(Editoperation eop, bool proteineop)
{
  Editoperation maxlen = proteineop ? MAXIDENTICALLENGTH_PROTEIN
                                    : MAXIDENTICALLENGTH;
  if (eop & maxlen)
    return eop & maxlen;
  return 1;
}

void gt_editoperation_set_length(Editoperation *eop, unsigned int length,
                              bool proteineop)
{
  Editoperation maxlen = proteineop ? MAXIDENTICALLENGTH_PROTEIN
                                    : MAXIDENTICALLENGTH;
  gt_assert(eop && length);
  gt_assert(length <= maxlen);
  *eop = (*eop & ~maxlen) | length;
}

static void showmultieop(const char *eop_string, unsigned int eop_length,
                         bool proteineop, unsigned int indentlevel,
                         GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  if (proteineop) {
    gt_file_xprintf(outfp, "<Protein_eop_type>");
    gt_file_xprintf(outfp, "%s", eop_string);
    gt_file_xprintf(outfp, "</Protein_eop_type>\n");
  }
  else {
    gt_file_xprintf(outfp, "<DNA_eop_type>");
    gt_file_xprintf(outfp, "%s", eop_string);
    gt_file_xprintf(outfp, "</DNA_eop_type>\n");
  }

  gth_indent(outfp, indentlevel);
  if (proteineop)
  {
    gt_file_xprintf(outfp, "<Protein_eop_length>");
    gt_file_xprintf(outfp, "%u", eop_length);
    gt_file_xprintf(outfp, "</Protein_eop_length>\n");
  }
  else {
    gt_file_xprintf(outfp, "<DNA_eop_length>");
    gt_file_xprintf(outfp, "%u", eop_length);
    gt_file_xprintf(outfp, "</DNA_eop_length>\n");
  }
}

static void showoneeditopgeneric(GtFile *outfp, Editoperation eop,
                                 bool proteineop,
                                 bool xmlout, unsigned long indentlevel,
                                 bool nexteopisdefined, Editoperation nexteop,
                                 unsigned long *consecutive_eop_length)
{
  unsigned int eop_length;

  /* this code is necessary to show consecutive edit operations of the same
     type as one edit operation in the XML output */
  if (xmlout) {
    if (nexteopisdefined &&
        gt_editoperation_type(eop, proteineop) ==
        gt_editoperation_type(nexteop, proteineop)) {
      /* store length of this eop */
      *consecutive_eop_length += gt_editoperation_length(eop, proteineop);
      /* return, this consecutive eop is shown later */
      return;
    }
    else {
      /* store total length of this consecutive eop for output */
      eop_length = *consecutive_eop_length +
                   gt_editoperation_length(eop, proteineop);
      /* reset */
      *consecutive_eop_length = 0;
    }
  }
  else
    eop_length = gt_editoperation_length(eop, proteineop);

  switch (gt_editoperation_type(eop, proteineop)) {
    case EOP_TYPE_MATCH:
      if (xmlout)
        showmultieop("match", eop_length, proteineop, indentlevel, outfp);
      else
        gt_file_xprintf(outfp, "(M %u)", eop_length);
      break;
    case EOP_TYPE_INTRON:
      if (xmlout)
        showmultieop("intron", eop_length, proteineop, indentlevel, outfp);
      else {
        gt_file_xprintf(outfp, "(Intron%s %u)", proteineop ? "(0)" : "",
                        eop_length);
      }
      break;
    case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      if (xmlout) {
        showmultieop("intron_with_1_base_left", eop_length, proteineop,
                     indentlevel, outfp);
      }
      else
        gt_file_xprintf(outfp, "(Intron(1) %u)", eop_length);
      break;
    case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
      if (xmlout) {
        showmultieop("intron_with_2_bases_left", eop_length, proteineop,
                     indentlevel, outfp);
      }
      else
        gt_file_xprintf(outfp, "(Intron(2) %u)", eop_length);
      break;
    case EOP_TYPE_MISMATCH:
      if (xmlout)
        showmultieop("mismatch", eop_length, proteineop, indentlevel, outfp);
      else
        gt_file_xfputc('R',outfp);
      break;
    case EOP_TYPE_DELETION:
      if (xmlout)
        showmultieop("deletion", eop_length, proteineop, indentlevel, outfp);
      else
        gt_file_xfputc('D',outfp);
      break;
    case EOP_TYPE_INSERTION:
      if (xmlout)
        showmultieop("insertion", eop_length, proteineop, indentlevel, outfp);
      else
        gt_file_xfputc('I',outfp);
      break;
    case EOP_TYPE_MISMATCH_WITH_1_GAP:
      if (xmlout) {
        showmultieop("mismatch_with_1_gap", eop_length, proteineop, indentlevel,
                     outfp);
      }
      else
        gt_file_xprintf(outfp, "R1");
      break;
    case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      if (xmlout) {
        showmultieop("mismatch_with_2_gaps", eop_length, proteineop,
                     indentlevel, outfp);
      }
      else
        gt_file_xprintf(outfp, "R2");
      break;
    case EOP_TYPE_DELETION_WITH_1_GAP:
      if (xmlout) {
        showmultieop("deletion_with_1_gap", eop_length, proteineop, indentlevel,
                     outfp);
      }
      else
        gt_file_xprintf(outfp, "D1");
      break;
    case EOP_TYPE_DELETION_WITH_2_GAPS:
      if (xmlout) {
        showmultieop("deletion_with_2_gaps", eop_length, proteineop,
                     indentlevel, outfp);
      }
      else
        gt_file_xprintf(outfp, "D2");
      break;
    default: gt_assert(0);
  }
  if (!xmlout)
    gt_file_xfputc(' ', outfp);
}

void gt_editoperation_show(Editoperation *eops, unsigned long num_of_eops,
                        bool proteineops, bool xmlout, unsigned int indentlevel,
                        GtFile *outfp)
{
  long i;
  unsigned long consecutive_eop_length = 0;

  if (xmlout) {
    gth_indent(outfp, indentlevel);
    if (proteineops)
      gt_file_xprintf(outfp, "<Protein_eops>\n");
    else
      gt_file_xprintf(outfp, "<DNA_eops>\n");
  }

  for (i=(long) (num_of_eops - 1); i >= 0; i--) {
    if (i > 0) {
      showoneeditopgeneric(outfp, eops[i], proteineops, xmlout, indentlevel + 1,
                           true, eops[i-1], &consecutive_eop_length);
    }
    else {
      showoneeditopgeneric(outfp, eops[i], proteineops, xmlout, indentlevel + 1,
                           false, 0, &consecutive_eop_length);

    }
  }

  if (xmlout) {
    gth_indent(outfp, indentlevel);
    if (proteineops)
      gt_file_xprintf(outfp, "</Protein_eops>\n");
    else
      gt_file_xprintf(outfp, "</DNA_eops>\n");
  }
  else
    gt_file_xfputc('\n', outfp);
}
