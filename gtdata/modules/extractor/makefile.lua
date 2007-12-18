--[[
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

module(..., package.seeall)

require "extractor.file"

Makefile = extractor.File:new()

Makefile.filename = "Makefile"
Makefile.basename = "Makefile"
Makefile.filecontent = [[
CC:=gcc
CFLAGS:=-Wall -Werror -Os -pipe
LD:=$(CC)

SRC:=$(wildcard *.c)
OBJ:=$(SRC:%.c=%.o)
DEP:=$(SRC:%.c=%.d)

all: prog

prog: $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -lm -o $@

# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(DEP)

.PHONY: clean
clean:
	rm -f *.o
]]
