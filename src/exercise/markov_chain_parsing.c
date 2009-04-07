/*
  Copyright (c) 2004, 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004, 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/fa.h"
#include "core/parseutils.h"
#include "core/splitter.h"
#include "core/str.h"
#include "exercise/markov_chain_parsing.h"

#define FIRSTLINE  "# Number of states:"
#define THIRDLINE  "# State labels:"
#define FIFTHLINE  "# Transition matrix:"

static int read_next_line(GtStr *line, FILE *fp, const char *filename,
                          GtError *err)
{
  gt_error_check(err);
  gt_assert(line && fp);
  gt_str_reset(line);
  if (gt_str_read_next_line(line, fp) == EOF) {
    gt_error_set(err, "unexpected end of file \"%s\"", filename);
    return -1;
  }
  return 0;
}

static int scan_alphabet(GtStr *alphabet, const char *line, int num_of_states,
                         GtError *err)
{
  gt_error_check(err);
  gt_assert(alphabet && line);
  while (*line != '\0') {
    if (*line != ' ')
      gt_str_append_char(alphabet, *line);
    line++;
  }
  if (gt_str_length(alphabet) != num_of_states) {
    gt_error_set(err, "number of states = %d != %lu number of state labels",
              num_of_states, gt_str_length(alphabet));
    return -1;
  }
  return 0;
}

static GtMarkovChain* parse_gt_markov_chain_file(FILE *fp, const char *filename,
                                            GtError *err)
{
  unsigned long i, j;
  GtMarkovChain *mc = NULL;
  GtStr *line, *alphabet;
  GtSplitter *splitter;
  int num_of_states, had_err;
  gt_error_check(err);
  gt_assert(fp);
  line = gt_str_new();
  alphabet = gt_str_new();
  splitter = gt_splitter_new();

  /* read first line */
  had_err = read_next_line(line, fp, filename, err);

  /* check first line */
  if (!had_err && strcmp(gt_str_get(line), FIRSTLINE)) {
    gt_error_set(err, "first line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* read second line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* parse number of states */
  if (!had_err && gt_parse_int(&num_of_states, gt_str_get(line))) {
    gt_error_set(err, "could not parse number of states from file \"%s\"",
              filename);
    had_err = -1;
  }

  /* read third line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* check third line */
  if (!had_err && strcmp(gt_str_get(line), THIRDLINE)) {
    gt_error_set(err, "third line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* read fourth line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* determine alphabet from line */
  if (!had_err)
    had_err = scan_alphabet(alphabet, gt_str_get(line), num_of_states, err);

  /* read fifth line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* check fifth line */
  if (!had_err && strcmp(gt_str_get(line), FIFTHLINE)) {
    gt_error_set(err, "fifth line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* create markov chain object */
  if (!had_err)
    mc = gt_markov_chain_new(gt_str_get(alphabet));

  /* read in transition probabilities */
  for (i = 0; !had_err && i < num_of_states; i++) {
    double transition_prob;
    had_err = read_next_line(line, fp, filename, err);
    if (!had_err) {
      gt_splitter_split(splitter, gt_str_get(line), gt_str_length(line), ' ');
      if (gt_splitter_size(splitter) != num_of_states) {
        gt_error_set(err, "%lu line of file \"%s\" does not contain %d tokens",
                  5 + i + 1, filename, num_of_states);
        had_err = -1;
      }
    }
    for (j = 0; !had_err && j < num_of_states; j++) {
      if (gt_parse_double(&transition_prob,
                          gt_splitter_get_token(splitter, j))) {
        gt_error_set(err, "could not parse transition probability %lu on line "
                          "%lu from file \"%s\"", j + 1, 5 + i + 1, filename);
        had_err = -1;
      }
      if (!had_err) /* store transition probability */
        gt_markov_chain_set_transition_prob(mc, i, j, transition_prob);
    }
    gt_splitter_reset(splitter);
  }

  /* make sure markov chain is valid */
  if (!had_err && !gt_markov_chain_is_valid(mc)) {
    gt_error_set(err, "the markov chain defined in file \"%s\" is not valid",
              filename);
    had_err = -1;
  }

  gt_splitter_delete(splitter);
  gt_str_delete(alphabet);
  gt_str_delete(line);
  if (had_err) {
    gt_markov_chain_delete(mc);
    return NULL;
  }
  return mc;
}

GtMarkovChain* gt_markov_chain_parse(const char *filename, GtError *err)
{
  GtMarkovChain *mc;
  FILE *fp;
  gt_error_check(err);
  fp = gt_fa_xfopen(filename, "r");
  mc = parse_gt_markov_chain_file(fp, filename, err);
  gt_fa_xfclose(fp);
  return mc;
}
