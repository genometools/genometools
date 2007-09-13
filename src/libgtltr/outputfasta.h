/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef OUTPUTFASTA_H
#define OUTPUTFASTA_H

int showpredictionsmultiplefasta(const LTRharvestoptions *lo,
		       Seqpos *markpos,
		       bool innerregion,
		       unsigned int linewidth,
                       Sequentialsuffixarrayreader *ssar,
		       bool showseqnum,
		       Env *env);

#endif
