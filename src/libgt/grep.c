/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "ensure.h"
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

bool grep(const char *pattern, const char *line)
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
    return false;
  return true;
}

int grep_unit_test(void)
{
  ensure( grep("a", "a"));
  ensure(!grep("b", "a"));
  ensure( grep("aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA"));
  ensure(!grep("aba", "wenbapzbpqSayhzzaBaZZqyghaAAahhaA"));
  ensure( grep("^aba", "abawenbapzbpqSayhzzZZqyghaAAahhaA"));
  ensure(!grep("^aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA"));
  return EXIT_SUCCESS;
}
