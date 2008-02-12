--[[
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

local makefile_template = [[
CC:=gcc
CFLAGS:=-Wall -Werror -Wunused-parameter -g -Os -pipe
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
	rm -f *.[od]
]]

Makefile = {}

function Makefile:new(progname)
  o = extractor.File:new("Makefile", true)
  if progname then
    o.progname = progname
    o.filecontent = makefile_template:gsub("prog", progname)
  else
    o.progname = "prog"
    o.filecontent = makefile_template
  end
  return o
end

