#ifndef ESA_SEQREAD_H
#define ESA_SEQREAD_H
#include <stdbool.h>
#include "libgtcore/str.h"
#include "libgtcore/env.h"
#include "seqpos-def.h"
#include "encseq-def.h"

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

Sequentialsuffixarrayreader *newSequentialsuffixarrayreader(
                                        const Str *indexname,
                                        unsigned int demand,
                                        bool mapped,
                                        Env *env);

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar,
                                     Env *env);

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Env *env);

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Env *env);

Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

Readmode readmodeSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

const Alphabet *alphabetSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr);

#endif
