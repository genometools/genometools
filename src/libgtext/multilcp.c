/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include <stdlib.h>
#include <string.h>
#include "libgtcore/array2dim.h"
#include "libgtext/multilcp.h"

int** multilcp_compute(const char *u, int m, const char *v, int n)
{
  int i, j, **prefix;

  array2dim_malloc(prefix, m, n);
  for (i = 0; i < m; i++) {
    if (u[m-1-i] == v[n-1])
      prefix[m-1-i][n-1] = 1;
    else
      prefix[m-1-i][n-1] = 0;
  }
  for (j = 0; j < n; j++) {
    if (u[m-1] == v[n-1-j])
      prefix[m-1][n-1-j] = 1;
    else
      prefix[m-1][n-1-j] = 0;
  }
  for (i = m-2; i >= 0; i--) {
    for (j = n-2; j >= 0; j--) {
      if (u[i] == v[j])
        prefix[i][j] = 1 + prefix[i+1][j+1];
      else
        prefix[i][j] = 0;
    }
  }
  return prefix;
}

void multilcp_show(int **tab, int dim1, int dim2)
{
  int i, j;

  for (i = 0; i < dim1; i++) {
    for (j = 0; j < dim2; j++) {
      printf("%d",tab[i][j]);
      if (j < (dim2-1))
        printf(" ");
    }
    printf("\n");
  }
}
