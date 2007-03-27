/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#include <libgtcore/env.h>
#include <libgtcore/option.h>

typedef struct OutputFileInfo OutputFileInfo;

OutputFileInfo* outputfileinfo_new(Env*);
void            outputfile_register_options(OptionParser*, GenFile **outfp,
                                            OutputFileInfo*, Env*);
void            outputfileinfo_delete(OutputFileInfo*, Env*);

#endif
