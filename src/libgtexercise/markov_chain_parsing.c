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
#include "libgtcore/fa.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/splitter.h"
#include "libgtcore/str.h"
#include "libgtexercise/markov_chain_parsing.h"

#define FIRSTLINE  "# Number of states:"
#define THIRDLINE  "# State labels:"
#define FIFTHLINE  "# Transition matrix:"

static int read_next_line(Str *line, FILE *fp, const char *filename, Error *err)
{
  error_check(err);
  assert(line && fp);
  str_reset(line);
  if (str_read_next_line(line, fp) == EOF) {
    error_set(err, "unexpected end of file \"%s\"", filename);
    return -1;
  }
  return 0;
}

static int scan_alphabet(Str *alphabet, const char *line, int num_of_states,
                         Error *err)
{
  error_check(err);
  assert(alphabet && line);
  while (*line != '\0') {
    if (*line != ' ')
      str_append_char(alphabet, *line);
    line++;
  }
  if (str_length(alphabet) != num_of_states) {
    error_set(err, "number of states = %d != %lu number of state labels",
              num_of_states, str_length(alphabet));
    return -1;
  }
  return 0;
}

static MarkovChain* parse_markov_chain_file(FILE *fp, const char *filename,
                                            Error *err)
{
  unsigned long i, j;
  MarkovChain *mc = NULL;
  Str *line, *alphabet;
  Splitter *splitter;
  int num_of_states, had_err;
  error_check(err);
  assert(fp);
  line = str_new();
  alphabet = str_new();
  splitter = splitter_new();

  /* read first line */
  had_err = read_next_line(line, fp, filename, err);

  /* check first line */
  if (!had_err && strcmp(str_get(line), FIRSTLINE)) {
    error_set(err, "first line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* read second line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* parse number of states */
  if (!had_err && parse_int(&num_of_states, str_get(line))) {
    error_set(err, "could not parse number of states from file \"%s\"",
              filename);
    had_err = -1;
  }

  /* read third line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* check third line */
  if (!had_err && strcmp(str_get(line), THIRDLINE)) {
    error_set(err, "third line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* read fourth line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* determine alphabet from line */
  if (!had_err)
    had_err = scan_alphabet(alphabet, str_get(line), num_of_states, err);

  /* read fifth line */
  if (!had_err)
    had_err = read_next_line(line, fp, filename, err);

  /* check fifth line */
  if (!had_err && strcmp(str_get(line), FIFTHLINE)) {
    error_set(err, "fifth line of file \"%s\" does not match standard",
              filename);
    had_err = -1;
  }

  /* create markov chain object */
  if (!had_err)
    mc = markov_chain_new(str_get(alphabet));

  /* read in transition probabilities */
  for (i = 0; !had_err && i < num_of_states; i++) {
    double transition_prob;
    had_err = read_next_line(line, fp, filename, err);
    if (!had_err) {
      splitter_split(splitter, str_get(line), str_length(line), ' ');
      if (splitter_size(splitter) != num_of_states) {
        error_set(err, "%lu line of file \"%s\" does not contain %d tokens",
                  5 + i + 1, filename, num_of_states);
        had_err = -1;
      }
    }
    for (j = 0; !had_err && j < num_of_states; j++) {
      if (parse_double(&transition_prob, splitter_get_token(splitter, j))) {
        error_set(err, "could not parse transition probability %lu on line %lu "
                       "from file \"%s\"", j + 1, 5 + i + 1, filename);
        had_err = -1;
      }
      if (!had_err) /* store transition probability */
        markov_chain_set_transition_prob(mc, i, j, transition_prob);
    }
    splitter_reset(splitter);
  }

  /* make sure markov chain is valid */
  if (!had_err && !markov_chain_is_valid(mc)) {
    error_set(err, "the markov chain defined in file \"%s\" is not valid",
              filename);
    had_err = -1;
  }

  splitter_delete(splitter);
  str_delete(alphabet);
  str_delete(line);
  if (had_err) {
    markov_chain_delete(mc);
    return NULL;
  }
  return mc;
}

MarkovChain* markov_chain_parse(const char *filename, Error *err)
{
  MarkovChain *mc;
  FILE *fp;
  error_check(err);
  fp = fa_xfopen(filename, "r");
  mc = parse_markov_chain_file(fp, filename, err);
  fa_xfclose(fp);
  return mc;
}
