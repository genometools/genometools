#ifndef GENSTREAM_H
#define GENSTREAM_H
#include <stdio.h>
#include <stdbool.h>
#include <zlib.h>

typedef struct
{
  bool isgzippedstream;
  union
  {
    FILE *fopenstream;
    gzFile gzippedstream;
  } stream;
} Genericstream;

#endif
