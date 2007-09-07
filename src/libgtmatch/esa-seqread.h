#ifndef ESA_SEQREAD_H
#define ESA_SEQREAD_H
#include <stdbool.h>
#include "libgtcore/str.h"
#include "libgtcore/env.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "sarr-def.h"

typedef enum
{
  SEQ_mappedboth,
  SEQ_scan,
  SEQ_suftabfrommemory
} Sequentialaccesstype;

typedef struct Sequentialsuffixarrayreader Sequentialsuffixarrayreader;

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                        const Str *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        Env *env);

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromRAM(
                                        const Encodedsequence *encseq,
                                        const Seqpos *suftab,
                                        Seqpos numberofsuffixes,
                                        Readmode readmode,
                                        Env *env);

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar,
                                     Env *env);

int nextSequentiallcpvalue(Seqpos *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           Env *env);

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              Env *env);

const Encodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *sarr);

Readmode readmodeSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

const Alphabet *alphabetSequentialsuffixarrayreader(
			  const Sequentialsuffixarrayreader *sarr);

unsigned long numofdbsequencesSequentialsuffixarrayreader(
                    const Sequentialsuffixarrayreader *sarr);

#endif
