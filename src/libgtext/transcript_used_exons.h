/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_USED_EXONS_H
#define TRANSCRIPT_USED_EXONS_H

#include <gtcore.h>

typedef struct TranscriptUsedExons TranscriptUsedExons;

TranscriptUsedExons* transcript_used_exons_new(Env*);
Dlist*               transcript_used_exons_get_all(TranscriptUsedExons*);
Dlist*               transcript_used_exons_get_single(TranscriptUsedExons*);
Dlist*               transcript_used_exons_get_initial(TranscriptUsedExons*);
Dlist*               transcript_used_exons_get_internal(TranscriptUsedExons*);
Dlist*               transcript_used_exons_get_terminal(TranscriptUsedExons*);
void                 transcript_used_exons_delete(TranscriptUsedExons*, Env*);

#endif
