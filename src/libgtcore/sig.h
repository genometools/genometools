/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SIG_H
#define SIG_H

/* register a signal handler for all signals which can be caught */
void sig_register_all(void (*)(int));
void sig_unregister_all(void);

#endif
