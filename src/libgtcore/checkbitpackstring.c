/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**
** See LICENSE file or http://genometools.org/license.html for license details.
**
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/bitpackstring.h>
#include <libgtcore/env.h>
#include <libgtcore/ensure.h>

int
bitPackString_unit_test(Env *env)
{
  return bitPackStringInt_unit_test(env)
    || bitPackStringUInt_unit_test(env)
    || bitPackString8_unit_test(env)
    || bitPackString16_unit_test(env)
    || bitPackString32_unit_test(env)
    || bitPackString64_unit_test(env);
}

