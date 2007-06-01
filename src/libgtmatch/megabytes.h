#ifndef MEGABYTES_H
#define MEGABYTES_H

/*
  The following macro transforms bytes into megabytes.
*/

#define MEGABYTES(V)  ((double) (V)/((((unsigned long) 1) << 20) - 1))

#endif
