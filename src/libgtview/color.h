/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file color.h
 * \author Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
 */
#ifndef COLOR_H
#define COLOR_H

#include <stdbool.h>

typedef struct {
  double red, green, blue;
} Color;

bool color_equals(Color, Color);

#endif
