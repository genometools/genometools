/*
  Copyright (c) 2007-2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef STYLE_H
#define STYLE_H

#include "lua.h"
#include "annotationsketch/style_api.h"
#include "extended/genome_node.h"

/* Creates a GtStyle object which reuses the given Lua state
   instead of creating a new one. */
GtStyle*       gt_style_new_with_state(lua_State*);

/* Identical to gt_style_get_color(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_color_with_track(const GtStyle *style,
                                                 const char *section,
                                                 const char *key,
                                                 GtColor *result,
                                                 GtFeatureNode *fn,
                                                 const GtStr *track_id,
                                                 GtError *err);
/* Retrieves a string value from <style> for key <key> in section <section>.
   The string is written to the <GtStr> object <result>, overwriting its prior
   contents. Optionally, a feature node pointer <fn> can be specified for
   handling in node-specific callbacks.
   Because color definitions can be functions, <gt_style_get_str()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the string was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing string, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_str(const GtStyle *style, const char *section,
                                    const char *key, GtStr *result,
                                    GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_str(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_str_with_track(const GtStyle *style,
                                               const char *section,
                                               const char *key,
                                               GtStr *result,
                                               GtFeatureNode *fn,
                                               const GtStr *track_id,
                                               GtError *err);

/* Retrieves a numeric value from <style> for key <key> in section <section>.
   The value is written to the location pointed to by <result>. Optionally, a
   feature node pointer <fn> can be specified for handling in node-specific
   callbacks.
   Because the definitions can be functions, <gt_style_get_num()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the number was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing number, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_num(const GtStyle *style, const char *section,
                                    const char *key, double *result,
                                    GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_num(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_num_with_track(const GtStyle *style,
                                               const char *section,
                                               const char *key,
                                               double *result,
                                               GtFeatureNode *fn,
                                               const GtStr *track_id,
                                               GtError *err);

/* Retrieves a boolean value from <style> for key <key> in section <section>.
   The value is written to the location pointed to by <result>. Optionally, a
   feature node pointer <fn> can be specified for handling in node-specific
   callbacks.
   Because the definitions can be functions, <gt_style_get_bool()> can fail at
   runtime. In this case, this function returns GT_STYLE_QUERY_ERROR and <err>
   is set accordingly.
   If the value was not specified in <style>, <result> is left untouched and
   GT_STYLE_QUERY_NOT_SET is returned so the caller can handle this case.
   In case of successful retrieval of an existing boolean, GT_STYLE_QUERY_OK
   is returned. */
GtStyleQueryStatus gt_style_get_bool(const GtStyle *style, const char *section,
                                     const char *key, bool *result,
                                     GtFeatureNode *fn, GtError *err);
/* Identical to gt_style_get_bool(), except that it also takes a <track_id>
   which is passed to a potential callback function in the style file. */
GtStyleQueryStatus gt_style_get_bool_with_track(const GtStyle *style,
                                                const char *section,
                                                const char *key,
                                                bool *result,
                                                GtFeatureNode *fn,
                                                const GtStr *track_id,
                                                GtError *err);

/* Set numeric value of key <key> in <section> to <number>. */
void     gt_style_set_num_p(GtStyle*, const char *section, const char *key,
                            double* number);

int            gt_style_unit_test(GtError*);

/* Deletes a GtStyle object but leaves the internal Lua state intact. */
void           gt_style_delete_without_state(GtStyle*);

#endif
