/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/cstr_array.h"
#include "gth/indent.h"
#include "gth/showbool.h"
#include "gth/time.h"
#include "gth/call_info.h"
#include "gth/gthspeciestab.h"
#include "gth/run_header.h"
#include "gt_config.h"

/* The name of the splice site models. */
#define SPLICE_SITE_MODEL_NAME  "Bayesian"

/* The name of the generic species */
#define GENERIC_SPECIES_NAME    "generic"

static void show_overall_reference_type(GthAlphatype overallalphatype,
                                        unsigned int indentlevel,
                                        GtFile *outfp)
{
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<overall_reference_type>");
  switch (overallalphatype)
  {
    case DNA_ALPHA:
      gt_file_xprintf(outfp, "ESTcDNA");
      break;
    case PROTEIN_ALPHA:
      gt_file_xprintf(outfp, "Protein");
      break;
    case MIXED_ALPHA:
      gt_file_xprintf(outfp, "Mixed");
      break;
    default: gt_assert(0);
  }
  gt_file_xprintf(outfp, "</overall_reference_type>\n");
}

static void show_xml_run_header(GthCallInfo *call_info, GthInput *input,
                                const char *timestring, const char *gth_version,
                                unsigned int indentlevel, const char **args)
{
  GtFile *outfp = call_info->out->outfp;
  unsigned long i;

  gth_indent(outfp, indentlevel);
  if (call_info->intermediate) {
    gt_file_xprintf(outfp, "<header xmlns=\"http://www.GenomeThreader.org/"
                       "SplicedAlignment/header/\">\n");
  }
  else {
    gt_file_xprintf(outfp,
              "<header xmlns=\"http://www.genomethreader.org/GTH_output/"
              "header/\">\n");
  }

  /* at least one genomic file defined */
  gt_assert(gth_input_num_of_gen_files(input));
  /* at least one reference file defined */
  gt_assert(gth_input_num_of_ref_files(input));

  /* show a readable version of GthCallInfo. That is, it is shown with wich
     parameters the program was called */

  indentlevel++;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<source program=\"GenomeThreader\" version=\"%s\" "
                  "build_date=\"%s\" run_date=\"%s\"/>\n", gth_version,
                  GT_BUILT, timestring);

  /* show genomic file names */
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<gDNA_template_files>\n");
  indentlevel++;

  for (i = 0; i < gth_input_num_of_gen_files(input); i++) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<temp_name>%s</temp_name>\n",
                    gth_input_get_genomic_filename(input, i));
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</gDNA_template_files>\n");

  /* show reference file names */
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<reference_files>\n");
  indentlevel++;

  for (i = 0; i < gth_input_num_of_ref_files(input); i++) {
    gth_indent(outfp, indentlevel);
    gt_file_xprintf(outfp, "<file ref_name=\"%s\" type=\"%s\"/>\n",
                    gth_input_get_reference_filename(input, i),
                    gth_input_get_alphatype(input, i) == DNA_ALPHA
                    ? "ESTcDNA" : "Protein");
  }

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</reference_files>\n");

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<splice_site_parameters parameter_type=\"%s\" "
                  "species=\"%s\"/>\n", SPLICE_SITE_MODEL_NAME,
                  call_info->speciesnum ==  NUMOFSPECIES
                  ? GENERIC_SPECIES_NAME
                  : speciestab[call_info->speciesnum]);

  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<parameters>\n");
  indentlevel++;

  /* output name of BSSM file */
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<parameter name=\"bssmfile\" value=\"%s\"/>\n",
                  gth_input_bssmfilename(input));

  /* output name of scorematrix */
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<parameter name=\"scorematrixfile\" value=\"%s\"/>\n",
                  gt_str_get(call_info->scorematrixfile));

  /* output searchmode */
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "<parameter name=\"searchmode\" "
                  "value=\"forward=%s,reverse=%s)\"/>\n",
                  GTH_SHOWBOOL(gth_input_forward(input)),
                  GTH_SHOWBOOL(gth_input_reverse(input)));

  /* output arguments as comment */
  gt_file_xprintf(outfp, "<!--\n%c Arguments: ", COMMENTCHAR);
  gt_cstr_array_show_genfile(args, outfp);
  gt_file_xprintf(outfp, "-->\n");

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</parameters>\n");

  show_overall_reference_type(gth_input_overall_alphatype(input),
                              indentlevel, outfp);

  indentlevel--;
  gth_indent(outfp, indentlevel);
  gt_file_xprintf(outfp, "</header>\n");
}

void gth_run_header_show(GthCallInfo *call_info, GthInput *input,
                         const char *gth_version, unsigned int indentlevel,
                         const char **args)
{
  char *timestring;
  GtFile *outfp = call_info->out->outfp;

  /* determine time */
  timestring = gth_get_time();

  /* output XML header */
  if (call_info->out->xmlout) {
    show_xml_run_header(call_info, input, timestring, gth_version, indentlevel,
                        args);
  }
  else if (!call_info->out->gff3out) {
    gt_file_xprintf(outfp, "%c GenomeThreader %s (%s)\n", COMMENTCHAR,
                    gth_version, GT_BUILT);
    gt_file_xprintf(outfp, "%c Date run: %s\n", COMMENTCHAR, timestring);
    gt_file_xprintf(outfp, "%c Arguments: ", COMMENTCHAR);
    gt_cstr_array_show_genfile(args, outfp);
  }

  /* free */
  gt_free(timestring);
}
