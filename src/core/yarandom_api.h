/* xscreensaver, Copyright (c) 1997, 1998, 2003 by Jamie Zawinski <jwz@jwz.org>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 */

#ifndef YARANDOM_API_H
#define YARANDOM_API_H

/* Yarandom module */

#undef random
#undef rand
#undef drand48
#undef srandom
#undef srand
#undef srand48
#undef frand
#undef RAND_MAX

/* Return a random number. */
unsigned int gt_ya_random(void);
/* Initialize random number generator using given seed. */
unsigned int gt_ya_rand_init(unsigned int);
/* Clean up static data for random number generator. */
void         gt_ya_rand_clean(void);

/* Maximum random number (2147483647) */
#define GT_RAND_MAX \
        0x7FFFFFFF
/* Return random number up to RAND_MAX. */
#define random() \
        ((GtWord) (gt_ya_random() & GT_RAND_MAX))

/*#define srandom(i) ya_rand_init(0)*/

/* Define these away to keep people from using the wrong APIs in GenomeTools.
 */
#define rand          __ERROR_use_random_not_rand_in_GenomeTools__
#define drand48       __ERROR_use_frand_not_drand48_in_GenomeTools__
#define srandom       __ERROR_do_not_call_srandom_in_GenomeTools__
#define srand         __ERROR_do_not_call_srand_in_GenomeTools__
#define srand48       __ERROR_do_not_call_srand48_in_GenomeTools__

#if defined (__GNUC__) && (__GNUC__ >= 2)
 /* Implement frand using GCC's statement-expression extension. */

# define gt_frand(f)                                                    \
  __extension__                                                         \
  ({ double tmp = ((((double) random()) * ((double) (f))) /             \
                   ((double) ((unsigned int)~0)));                      \
     tmp < 0 ? (-tmp) : tmp; })

#else /* not GCC2 - implement frand using a global variable.*/

/*@unused@*/ static double _frand_tmp_;
# define gt_frand(f)                                                    \
  (_frand_tmp_ = ((((double) random()) * ((double) (f))) /              \
                  ((double) ((unsigned int)~0))),                       \
   _frand_tmp_ < 0 ? (-_frand_tmp_) : _frand_tmp_)

#endif /* not GCC2 */

#endif
