/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_BITTABS_H
#define TRANSCRIPT_BITTABS_H

#include <libgt/bittab.h>
#include <libgt/env.h>

/* a container class for transcript bittabs */
typedef struct TranscriptBittabs TranscriptBittabs;

/* create an empy container */
TranscriptBittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal, Env*);

/* return the bittab for all exons */
Bittab*            transcript_bittabs_get_all(const TranscriptBittabs*);

/* return the bittab for single exons */
Bittab*            transcript_bittabs_get_single(const TranscriptBittabs*);

/* return the bittab for initial exons */
Bittab*            transcript_bittabs_get_initial(const TranscriptBittabs*);

/* return the bittab for internal exons */
Bittab*            transcript_bittabs_get_internal(const TranscriptBittabs*);

/* return the bittab for terminal exons */
Bittab*            transcript_bittabs_get_terminal(const TranscriptBittabs*);

void               transcript_bittabs_delete(TranscriptBittabs*, Env*);

#endif
