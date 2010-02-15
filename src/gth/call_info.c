/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/cstr_api.h"
#include "core/ma.h"
#include "gth/default.h"
#include "gth/call_info.h"
#include "gth/gthspeciestab.h"

GthCallInfo* gth_call_info_new(const char *progname)
{
  GthCallInfo *call_info;
  call_info = gt_malloc(sizeof (GthCallInfo));
  call_info->out = gthoutput_new();

  call_info->progname                 = gt_cstr_dup(progname);
  call_info->scorematrixfile          = gt_str_new();
  call_info->dp_options_core          = gth_dp_options_core_new();
  call_info->dp_options_est           = gth_dp_options_est_new();
  call_info->dp_options_postpro       = gth_dp_options_postpro_new();
  call_info->translationtable         = DEFAULT_TRANSLATIONTABLE;
  call_info->out->skipalignmentout    = DEFAULT_SKIPALIGNMENTOUT;
  call_info->speciesnum               = NUMOFSPECIES;

  /* set the DP options to default values */
  call_info->out->showverbose         = NULL;
  call_info->out->showverboseVM       = NULL;

  call_info->out->comments            = DEFAULT_COMMENTS;
  call_info->out->verboseseqs         = DEFAULT_VERBOSESEQS;
  call_info->out->xmlout              = DEFAULT_XMLOUT;
  call_info->out->gff3out             = false;

  /* init the spliced alignment filter */
  call_info->sa_filter = gth_sa_filter_new();

  /* init the splice site model */
  call_info->splice_site_model = gth_splice_site_model_new();

  return call_info;
}

void gth_call_info_delete(GthCallInfo *call_info)
{
  if (!call_info) return;

  gt_free(call_info->progname);
  gt_str_delete(call_info->scorematrixfile);
  gth_dp_options_core_delete(call_info->dp_options_core);
  gth_dp_options_est_delete(call_info->dp_options_est);
  gth_dp_options_postpro_delete(call_info->dp_options_postpro);

  /* free the spliced alignment filter */
  gth_sa_filter_delete(call_info->sa_filter);

  /* free the splice site model */
  gth_splice_site_model_delete(call_info->splice_site_model);

  /* free the output data structure */
  gthoutput_delete(call_info->out);

  gt_free(call_info);
}
