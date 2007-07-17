#
# Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg 
# See LICENSE file or http://genometools.org/license.html for license details
#

CC:=gcc
CXX:=g++
INCLUDEOPT:= -I$(CURDIR)/src -I$(CURDIR)/obj \
             -I$(CURDIR)/src/external/lua-5.1.2/src \
             -I$(CURDIR)/src/external/expat-2.0.0/lib\
             -I$(CURDIR)/src/external/bzip2-1.0.4\
             -I$(CURDIR)/src/external/agg-2.4/include\
             -I$(CURDIR)/src/external/libpng-1.2.18\
             -I/usr/include/cairo\
             -I/usr/local/include/cairo

CFLAGS:=
CXXFLAGS:=
GT_CFLAGS:= -g -Wall -Werror -pipe -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
            $(INCLUDEOPT)
GT_CXXFLAGS:= -g -pipe $(INCLUDEOPT)
STEST_FLAGS:=
LDFLAGS:=-L/usr/local/lib -L/usr/X11R6/lib
LDLIBS:=-lm -lz

# try to set RANLIB automatically
SYSTEM:=$(shell uname -s)
ifeq ($(SYSTEM),Darwin)
  RANLIB:=ranlib
  GT_CFLAGS+=-DHAVE_LLABS
endif
ifeq ($(SYSTEM),Linux)
  GT_CLAGS+=-D_ISOC99_SOURCE -D_XOPEN_SOURCE=600 -DHAVE_LLABS
endif

# the default GenomeTools libraries which are build
GTLIBS:=lib/libgtext.a\
        lib/libgtcore.a\
        lib/libbz2.a

# the core GenomeTools library (no other dependencies)
LIBGTCORE_SRC:=$(notdir $(wildcard src/libgtcore/*.c))
LIBGTCORE_OBJ:=$(LIBGTCORE_SRC:%.c=obj/%.o)

# the extended GenomeTools library (e.g., depends on Lua)
LIBGTEXT_C_SRC:=$(notdir $(wildcard src/libgtext/*.c))
LIBGTEXT_C_OBJ:=$(LIBGTEXT_C_SRC:%.c=obj/%.o)
LIBGTEXT_CXX_SRC:=$(notdir $(wildcard src/libgtext/*.cxx))
LIBGTEXT_CXX_OBJ:=$(LIBGTEXT_CXX_SRC:%.cxx=obj/%.o)

# the GenomeTools view library
LIBGTVIEW_C_SRC:=$(notdir $(wildcard src/libgtview/*.c))
LIBGTVIEW_C_OBJ:=$(LIBGTVIEW_C_SRC:%.c=obj/%.o)

TOOLS_SRC:=$(notdir $(wildcard src/tools/*.c))
TOOLS_OBJ:=$(TOOLS_SRC:%.c=obj/%.o)

LIBAGG_SRC:=$(notdir $(wildcard src/external/agg-2.4/src/*.cpp src/external/agg-2.4/src/ctrl/*.cpp))
LIBAGG_OBJ:=$(LIBAGG_SRC:%.cpp=obj/%.o)

LIBEXPAT_OBJ:=obj/xmlparse.o obj/xmlrole.o obj/xmltok.o

LIBLUA_OBJ=obj/lapi.o obj/lcode.o obj/ldebug.o obj/ldo.o obj/ldump.o \
           obj/lfunc.o obj/lgc.o obj/llex.o obj/lmem.o obj/lobject.o \
           obj/lopcodes.o obj/lparser.o obj/lstate.o obj/lstring.o   \
           obj/ltable.o obj/ltm.o obj/lundump.o obj/lvm.o obj/lzio.o \
           obj/lauxlib.o obj/lbaselib.o obj/ldblib.o obj/liolib.o    \
           obj/lmathlib.o obj/loslib.o obj/ltablib.o obj/lstrlib.o   \
           obj/loadlib.o obj/linit.o

LIBPNG_OBJ=obj/png.o obj/pngset.o obj/pngget.o obj/pngrutil.o obj/pngtrans.o \
           obj/pngwutil.o obj/pngread.o obj/pngrio.o obj/pngwio.o            \
           obj/pngwrite.o obj/pngrtran.o obj/pngwtran.o obj/pngmem.o         \
           obj/pngerror.o obj/pngpread.o

LIBRNV_OBJ := obj/rn.o obj/rnc.o obj/rnd.o obj/rnl.o obj/rnv.o obj/rnx.o obj/drv.o  \
              obj/ary.o obj/xsd.o obj/xsd_tm.o obj/dxl.o obj/dsl.o obj/sc.o obj/u.o \
              obj/ht.o obj/er.o obj/xmlc.o obj/s.o obj/m.o obj/rx.o

LIBBZ2_OBJ := obj/blocksort.o obj/huffman.o obj/crctable.o obj/randtable.o \
              obj/compress.o obj/decompress.o obj/bzlib.o

SERVER=gordon@genometools.org
WWWBASEDIR=/var/www/servers

# process arguments
ifneq ($(opt),no)
  GT_CFLAGS += -Os
  GT_CXXFLAGS += -Os
endif

ifeq ($(assert),no)
  GT_CFLAGS += -DNDEBUG
  GT_CXXFLAGS += -DNDEBUG
endif

ifeq ($(static),yes)
  LDFLAGS += -static
endif

ifeq ($(libgtview),yes)
  GTLIBS := $(GTLIBS) lib/libgtview.a
  GT_CFLAGS += -DLIBGTVIEW
  LDLIBS += -lcairo
  STEST_FLAGS += -libgtview
endif

# set prefix for install target
prefix ?= /usr/local

all: dirs $(GTLIBS) bin/skproto bin/gt bin/rnv

dirs:
	@test -d obj     || mkdir -p obj
	@test -d lib     || mkdir -p lib
	@test -d bin     || mkdir -p bin

lib/libexpat.a: $(LIBEXPAT_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBEXPAT_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libagg.a: $(LIBAGG_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBAGG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libbz2.a: $(LIBBZ2_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBBZ2_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtcore.a: obj/gt_build.h obj/gt_cc.h obj/gt_cflags.h obj/gt_version.h \
                 $(LIBGTCORE_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBGTCORE_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtext.a: $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtview.a: $(LIBGTVIEW_C_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBGTVIEW_C_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libpng.a: $(LIBPNG_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBPNG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/librnv.a: $(LIBRNV_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBRNV_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

bin/skproto: obj/skproto.o obj/gt_skproto.o lib/libgtcore.a lib/libbz2.a
	@echo "[link $@]"
	@$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

bin/gt: obj/gt.o obj/gtr.o $(TOOLS_OBJ) $(GTLIBS)
	@echo "[link $@]"
	@$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

bin/rnv: obj/xcl.o lib/librnv.a lib/libexpat.a
	@$(CC) $(LDFLAGS) $^ -o $@
	@echo "[link $@]"

obj/gt_build.h:
	@date +'#define GT_BUILT "%Y-%m-%d %H:%M:%S"' > $@

obj/gt_cc.h:
	@echo '#define GT_CC "'`$(CC) --version | head -n 1`\" > $@

obj/gt_cflags.h:
	@echo '#define GT_CFLAGS "$(CFLAGS) $(GT_CFLAGS)"' > $@

obj/gt_version.h: VERSION
	@echo '#define GT_VERSION "'`cat VERSION`\" > $@

Doxyfile: Doxyfile.in VERSION
	@echo '[rebuild $@]'
	@sed -e "s/@VERSION@/`cat VERSION`/g" <$< >$@

src/libgtcore/bitpackstringop8.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/checkbitpackstring8.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/bitpackstringop16.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/checkbitpackstring16.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/bitpackstringop32.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/checkbitpackstring32.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/bitpackstringop64.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

src/libgtcore/checkbitpackstring64.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl $^

# we create the dependency files on the fly
obj/%.o: src/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/libgtcore/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/libgtext/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/libgtext/%.cxx
	@echo "[compile $@]"
	@$(CXX) -c $< -o $@  $(CXXFLAGS) $(GT_CXXFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/libgtview/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/tools/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/external/agg-2.4/src/%.cpp
	@echo "[compile $@]"
	@$(CXX) -c $< -o $@  $(CXXFLAGS) $(GT_CXXFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/external/agg-2.4/src/ctrl/%.cpp
	@echo "[compile $@]"
	@$(CXX) -c $< -o $@  $(CXXFLAGS) $(GT_CXXFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/external/bzip2-1.0.4/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/external/expat-2.0.0/lib/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -DHAVE_MEMMOVE -MT $@ -MMD \
        -MP -MF $(@:.o=.d)

obj/%.o: src/external/lua-5.1.2/src/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -DLUA_USE_POSIX -MT $@ -MMD \
        -MP -MF $(@:.o=.d)

obj/%.o: src/external/libpng-1.2.18/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: src/external/rnv-1.7.8/%.c
	@echo "[compile $@]"
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) -DUNISTD_H="<unistd.h>" \
        -DEXPAT_H="<expat.h>" -DRNV_VERSION="\"1.7.8\"" -MT $@ -MMD -MP   \
        -MF $(@:.o=.d)

# read deps
-include obj/*.d

.SUFFIXES:
.PHONY: dist srcdist release gt install splint test clean cleanup apidoc

dist: all
	tar cvzf gt-`cat VERSION`.tar.gz bin/gt

srcdist:
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
        gzip -9 > genometools-`cat VERSION`.tar.gz

apidoc: Doxyfile
	doxygen

release: apidoc
	git tag "v`cat VERSION`"
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
	gzip -9 > genometools-`cat VERSION`.tar.gz
	scp genometools-`cat VERSION`.tar.gz $(SERVER):$(WWWBASEDIR)/genometools.org/htdocs/pub
	rsync -rv doc/api/html/ $(SERVER):$(WWWBASEDIR)/genometools.org/htdocs/apidoc
	git push --tags

installwww:
# install genometools.org website
	rsync -rv www/genometools.org/ $(SERVER):$(WWWBASEDIR)/genometools.org

gt: dirs bin/gt

install:
	test -d $(prefix)/bin || mkdir -p $(prefix)/bin
	cp bin/gt $(prefix)/bin
	cp -r gtdata $(prefix)/bin
	test -d $(prefix)/include/libgtcore || mkdir -p $(prefix)/include/libgtcore
	cp src/gtcore.h $(prefix)/include
	cp src/libgtcore/*.h $(prefix)/include/libgtcore
	test -d $(prefix)/include/libgtext || mkdir -p $(prefix)/include/libgtext
	cp src/gtext.h $(prefix)/include
	cp src/libgtext/*.h $(prefix)/include/libgtext
	cp src/gt.h $(prefix)/include
	test -d $(prefix)/lib || mkdir -p $(prefix)/lib
	cp lib/libgtcore.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtcore.a
endif
	cp lib/libgtext.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtext.a
endif

splint:
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) $(CURDIR)/src/*.c \
        $(CURDIR)/src/libgtcore/*.c \
        $(CURDIR)/src/libgtext/*.c \
        $(CURDIR)/src/tools/*.c

test: all
	bin/gt -test
	cd testsuite && env -i ruby -I. testsuite.rb -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin $(STEST_FLAGS)

clean:
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests

cleanup: clean
	rm -rf lib bin Doxyfile doc/api
