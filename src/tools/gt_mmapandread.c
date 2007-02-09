/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

int gt_mmapandread(int argc, char *argv[])
{
  int i, fd;
  void *map;
  struct stat sb;
  unsigned long j;
  char byte = 0;

  /* make sure at least one file is supplied */
  if (argc < 2) {
    fprintf(stderr, "Usage: %s file [...]\n", argv[0]);
    fprintf(stderr, "Map the supplied files into memory and read them once.\n");
    return EXIT_FAILURE;
  }

  /* iterate over all files */
  for (i = 1; i < argc; i++) {
    /* open file */
    fd = xopen(argv[i], O_RDONLY, 0);

    /* get file statistics */
    xfstat(fd, &sb);

    if (sb.st_size == 0)
      printf("file \"%s\" is empty\n", argv[i]);
    else if (!(sb.st_mode & S_IFREG))
      printf("\"%s\" is not a regular file\n", argv[i]);
    else {
      /* map file */
      map = xmmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

      /* read file */
      printf("reading file \"%s\"\n", argv[i]);
      j = 0;
      progressbar_start(&j, (unsigned long) sb.st_size);
      for (; j < (unsigned long) sb.st_size; j++)
        byte |= ((char*) map)[j];
      progressbar_stop();

      /* unmap file */
      xmunmap(map, sb.st_size);
    }

    /* close file */
    xclose(fd);
  }

  if (!byte)
    printf("all read files contained only null characters\n");

  return EXIT_SUCCESS;
}
