#
# Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg 
# See LICENSE file or http://genometools.org/license.html for license details
#

CCACHE:= $(shell type ccache >/dev/null && echo ccache)
CC:=$(CCACHE) gcc
LD:=gcc
INCLUDEOPT:= -I$(CURDIR)/src -I$(CURDIR)/obj
CFLAGS:=
GT_CFLAGS:= -Wall $(INCLUDEOPT)
LDFLAGS:=
LDLIBS:=-lm

VPATH:=$(VPATH):$(CURDIR)/src

LIBGT_SRC:=$(notdir $(wildcard src/libgt/*.c))
LIBGT_OBJ:=$(LIBGT_SRC:%.c=obj/%.o)

TOOLS_SRC:=$(notdir $(wildcard src/tools/*.c))
TOOLS_OBJ:=$(TOOLS_SRC:%.c=obj/%.o)

# process arguments
ifeq ($(opt),no)
  GT_CFLAGS += -g
else
  GT_CFLAGS += -Os
endif

ifeq ($(assert),no)
  GT_CFLAGS += -DNDEBUG
endif

# set prefix for install target
prefix ?= /usr/local

all: dirs lib/libgt.a bin/gt

dirs:
	@test -d include || mkdir -p include/libgt
	@test -d obj     || mkdir -p obj 
	@test -d lib     || mkdir -p lib 
	@test -d bin     || mkdir -p bin 

lib/libgt.a: obj/gt_build.h obj/gt_cc.h obj/gt_cflags.h obj/gt_version.h $(LIBGT_OBJ)
	ar ruv $@ $(LIBGT_OBJ)
ifdef RANLIB
	$(RANLIB) $@
endif
	ln src/gt.h include
	ln src/libgt/*.h include/libgt

bin/gt: obj/gt.o obj/gtr.o $(TOOLS_OBJ) lib/libgt.a
	$(LD) $(LDFLAGS) $^ $(LDLIBS) -o $@

obj/gt_build.h:
	@date +'#define GT_BUILT "%Y-%m-%d %H:%M:%S"' > $@

obj/gt_cc.h:
	@echo '#define GT_CC "'`$(CC) --version | head -n 1`\" > $@

obj/gt_cflags.h:
	@echo '#define GT_CFLAGS "$(CFLAGS) $(GT_CFLAGS)"' > $@

obj/gt_version.h: VERSION
	@echo '#define GT_VERSION "'`cat VERSION`\" > $@

# we create the dependency files on the fly
obj/%.o: src/%.c
	$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/libgt/%.c
	$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/tools/%.c
	$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read deps
-include obj/*.d

.PHONY: dist srcdist release gt libgt install splint test clean cleanup

dist: all
	tar cvzf gt-`cat VERSION`.tar.gz bin/gt_*

srcdist:
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
        gzip -9 > genometools-`cat VERSION`.tar.gz 

release:
	git tag "v`cat VERSION`"
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
        gzip -9 > genometools-`cat VERSION`.tar.gz 

gt: dirs bin/gt

libgt: dirs lib/libgt.a

install:
	test -d $(prefix)/bin || mkdir -p $(prefix)/bin
	cp bin/gt $(prefix)/bin
	test -d $(prefix)/include/libgt || mkdir -p $(prefix)/include/libgt
	cp src/gt.h $(prefix)/include	
	cp src/libgt/*.h $(prefix)/include/libgt
	test -d $(prefix)/lib || mkdir -p $(prefix)/lib
	cp lib/libgt.a $(prefix)/lib
	
splint:
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) $(CURDIR)/src/*.c

test: dirs bin/gt
	bin/gt -test
	cd testsuite && env -i ruby -I. testsuite.rb -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin

clean:
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_testsuite

cleanup: clean
	rm -rf include lib bin
