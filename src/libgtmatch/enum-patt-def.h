#ifndef ENUM_PATT_DEF_H
#define ENUM_PATT_DEF_H

typedef struct Enumpatternstate Enumpatternstate;

Enumpatternstate *newenumpattern(unsigned long minpatternlen,
                                 unsigned long maxpatternlen,
                                 const Encodedsequence *encseq,
                                 Env *env);

const Uchar *nextsampledpattern(unsigned long *patternlen,
                                Enumpatternstate *eps);

void freeenumpattern(Enumpatternstate *eps,Env *env);

#endif
