#ifndef MINMAX_H
#define MINMAX_H

/*
  This file defines macros for maximum and minimum computation,
  if they are not already defined.
*/

#ifndef MAX
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif

#endif
