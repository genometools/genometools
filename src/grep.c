/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "grep.h"

static void grep_error(int errcode, regex_t *matcher)
{
  char sbuf[BUFSIZ], *buf;
  size_t bufsize;
  bufsize = regerror(errcode, matcher, NULL, 0);
  buf = malloc(bufsize);
  (void) regerror(errcode, matcher, buf ? buf : sbuf, buf ? bufsize : BUFSIZ);
  error("grep(): %s", buf ? buf : sbuf);
}

unsigned int grep(const char *pattern, const char *line)
{
  regex_t matcher;
  int rval;
  assert(pattern && line);
  if ((rval = regcomp(&matcher, pattern, REG_EXTENDED | REG_NOSUB)))
    grep_error(rval, &matcher);
  rval = regexec(&matcher, line, 0, NULL, 0);
  if (rval && rval != REG_NOMATCH)
    grep_error(rval, &matcher);
  regfree(&matcher);
  if (rval)
    return 0;
  return 1;
}

int grep_unit_test(void)
{
  assert( grep("a", "a"));
  assert(!grep("b", "a"));
  assert( grep("aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA"));
  assert(!grep("aba", "wenbapzbpqSayhzzaBaZZqyghaAAahhaA"));
  assert( grep("^aba", "abawenbapzbpqSayhzzZZqyghaAAahhaA"));
  assert(!grep("^aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA"));
  return EXIT_SUCCESS;
}
