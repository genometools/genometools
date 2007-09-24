/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEARCHTSDANDMOTIF_H
#define SEARCHTSDANDMOTIF_H

int findcorrectboundaries(LTRharvestoptions *lo, LTRboundaries *boundaries,
    Sequentialsuffixarrayreader *ssar, const Seqpos *markpos, Env *env);

int searchformotifonlyinside(LTRharvestoptions *lo,
        LTRboundaries *boundaries,
        Sequentialsuffixarrayreader *ssar,
	const Seqpos *markpos,
	unsigned int *motifmismatchesleftLTR,
	unsigned int *motifmismatchesrightLTR,
	Env *env);

#endif
