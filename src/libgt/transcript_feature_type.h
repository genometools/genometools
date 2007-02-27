/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSCRIPT_FEATURE_TYPE_H
#define TRANSCRIPT_FEATURE_TYPE_H

typedef enum {
  TRANSCRIPT_FEATURE_TYPE_SINGLE,
  TRANSCRIPT_FEATURE_TYPE_INITIAL,
  TRANSCRIPT_FEATURE_TYPE_INTERNAL,
  TRANSCRIPT_FEATURE_TYPE_TERMINAL,
  TRANSCRIPT_FEATURE_TYPE_UNDETERMINED
} TranscriptFeatureType;

#endif
