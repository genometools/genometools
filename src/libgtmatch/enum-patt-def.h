#ifndef ENUM_PATT_DEF_H
#define ENUM_PATT_DEF_H

typedef struct Enumpatterniterator Enumpatterniterator;

Enumpatterniterator *newenumpatterniterator(unsigned long minpatternlen,
                                            unsigned long maxpatternlen,
                                            const Encodedsequence *encseq,
                                            Env *env);

const Uchar *nextEnumpatterniterator(unsigned long *patternlen,
                                     Enumpatterniterator *epi);

void freeEnumpatterniterator(Enumpatterniterator **epi,Env *env);

#endif
