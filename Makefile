#
# Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

INCLUDEOPT:=-I$(CURDIR)/src -I$(CURDIR)/obj \
            -I$(CURDIR)/src/external/zlib-1.2.3 \
            -I$(CURDIR)/src/external/md5-1.1.2/src \
            -I$(CURDIR)/src/external/lua-5.1.3/src \
            -I$(CURDIR)/src/external/luafilesystem-1.4.1/src \
            -I$(CURDIR)/src/external/lpeg-0.7 \
            -I$(CURDIR)/src/external/expat-2.0.1/lib \
            -I$(CURDIR)/src/external/bzip2-1.0.5 \
            -I$(CURDIR)/src/external/agg-2.4/include \
            -I$(CURDIR)/src/external/libpng-1.2.18 \
            -I$(CURDIR)/src/external/libtecla-1.6.1
# these variables are exported by the configuration script
CC:=gcc
CXX:=g++
EXP_CFLAGS:=$(CFLAGS)
EXP_LDFLAGS:=$(LDFLAGS)
EXP_CXXFLAGS:=$(CXXFLAGS)
EXP_CPPFLAGS:=$(CPPFLAGS)
EXP_LDLIBS:=$(LIBS) -lm
# ...while those starting with GT_ are for internal purposes only
GT_CFLAGS:=-g -Wall -Werror -Wunused-parameter -pipe -fPIC -Wpointer-arith
# expat needs -DHAVE_MEMMOVE
# lua needs -DLUA_USE_POSIX
# rnv needs -DUNISTD_H="<unistd.h>" -DEXPAT_H="<expat.h>" -DRNV_VERSION="\"1.7.8\""
# tecla needs -DHAVE_CURSES_H -DHAVE_TERM_H -DUSE_TERMINFO
EXT_FLAGS:= -DHAVE_MEMMOVE -DLUA_USE_POSIX -DUNISTD_H="<unistd.h>" \
            -DEXPAT_H="<expat.h>" -DRNV_VERSION=\"1.7.8\" \
            -DHAVE_CURSES_H -DHAVE_TERM_H -DUSE_TERMINFO
EXP_CPPFLAGS+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 $(EXT_FLAGS)
GT_CPPFLAGS:=$(INCLUDEOPT)
GT_CXXFLAGS:=-g -pipe
GT_LDFLAGS:=-Llib
STEST_FLAGS:=
EXP_LDFLAGS+=$(foreach dir, \
	$(shell test -d /usr/local/lib && echo /usr/local/lib ; \
	test -d /usr/X11R6/lib && echo /usr/X11R6/lib),-L$(dir))
BUILDSTAMP:=$(shell date +'"%Y-%m-%d %H:%M:%S"')

# try to set RANLIB automatically
SYSTEM:=$(shell uname -s)
  MACHINE:=$(shell uname -m)
ifeq ($(SYSTEM),Darwin)
  RANLIB:=ranlib
  SHARED:=-dynamiclib -undefined dynamic_lookup
  SHARED_OBJ_NAME_EXT:=.dylib
  ifeq ($(universal),yes)
    MACHINE:="Universal_Binary"
    GT_CFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
    GT_LDFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
  endif
else
  SHARED_OBJ_NAME_EXT:=.so
  SHARED:=-shared
endif

# the default GenomeTools libraries which are build
GTLIBS:=lib/libgtexercise.a\
        lib/libgtext.a\
        lib/libgtmgth.a\
        lib/libgtmatch.a\
        lib/libgtltr.a\
        lib/libgtcore.a\
        lib/libgtlua.a\
        lib/libexpat.a

# libraries for which we build replacements (that also appear in dependencies)
EXP_LDLIBS+=-lz -lbz2
OVERRIDELIBS:=lib/libbz2.a

# compiled executables
GTMAIN_SRC:=src/gt.c src/gtr.c src/gtt.c
GTMAIN_OBJ:=$(GTMAIN_SRC:%.c=obj/%.o)
GTMAIN_DEP:=$(GTMAIN_SRC:%.c=obj/%.d)

EXAMPLE_SRC:=src/example.c
EXAMPLE_OBJ:=$(EXAMPLE_SRC:%.c=obj/%.o)
EXAMPLE_DEP:=$(EXAMPLE_SRC:%.c=obj/%.d)

SKPROTO_SRC:=src/skproto.c src/tools/gt_skproto.c
SKPROTO_OBJ:=$(SKPROTO_SRC:%.c=obj/%.o)
SKPROTO_DEP:=$(SKPROTO_SRC:%.c=obj/%.d)

# the core GenomeTools library (no other dependencies)
AUTOGEN_LIBGTCORE_SRC:= src/libgtcore/bitpackstringop8.c \
	src/libgtcore/checkbitpackstring8.c \
	src/libgtcore/bitpackstringop16.c src/libgtcore/checkbitpackstring16.c \
	src/libgtcore/bitpackstringop32.c src/libgtcore/checkbitpackstring32.c \
	src/libgtcore/bitpackstringop64.c src/libgtcore/checkbitpackstring64.c
LIBGTCORE_SRC:=$(wildcard src/libgtcore/*.c)
LIBGTCORE_SRC:=$(filter-out $(AUTOGEN_LIBGTCORE_SRC), $(LIBGTCORE_SRC)) \
	 $(AUTOGEN_LIBGTCORE_SRC)
LIBGTCORE_OBJ:=$(LIBGTCORE_SRC:%.c=obj/%.o)
LIBGTCORE_DEP:=$(LIBGTCORE_SRC:%.c=obj/%.d)
LIBGTCORE_LIBDEP=-lbz2 -lz

# the extended GenomeTools library (e.g., depends on Lua)
LIBGTEXT_C_SRC:=$(wildcard src/libgtext/*.c)
LIBGTEXT_C_OBJ:=$(LIBGTEXT_C_SRC:%.c=obj/%.o)
LIBGTEXT_C_DEP:=$(LIBGTEXT_C_SRC:%.c=obj/%.d)
LIBGTEXT_CXX_SRC:=$(wildcard src/libgtext/*.cxx)
LIBGTEXT_CXX_OBJ:=$(LIBGTEXT_CXX_SRC:%.cxx=obj/%.o)
LIBGTEXT_CXX_DEP:=$(LIBGTEXT_CXX_SRC:%.cxx=obj/%.d)
LIBGTEXT_LIBDEP=-lgtcore -lbz2 -lz

# the exercise GenomeTools library
LIBGTEXERCISE_SRC:=$(wildcard src/libgtexercise/*.c)
LIBGTEXERCISE_OBJ:=$(LIBGTEXERCISE_SRC:%.c=obj/%.o)
LIBGTEXERCISE_DEP:=$(LIBGTEXERCISE_SRC:%.c=obj/%.d)

# the MetaGenomeThreader library
LIBGTMGTH_SRC:=$(wildcard src/libgtmgth/*.c)
LIBGTMGTH_OBJ:=$(LIBGTMGTH_SRC:%.c=obj/%.o)
LIBGTMGTH_DEP:=$(LIBGTMGTH_SRC:%.c=obj/%.d)

# the GenomeTools matching library
LIBGTMATCH_SRC:=$(wildcard src/libgtmatch/*.c)
LIBGTMATCH_OBJ:=$(LIBGTMATCH_SRC:%.c=obj/%.o)
LIBGTMATCH_DEP:=$(LIBGTMATCH_SRC:%.c=obj/%.d)

# the GenomeTools LTRharvest library
LIBGTLTR_SRC:=$(wildcard src/libgtltr/*.c)
LIBGTLTR_OBJ:=$(LIBGTLTR_SRC:%.c=obj/%.o)
LIBGTLTR_DEP:=$(LIBGTLTR_SRC:%.c=obj/%.d)

# the GenomeTools view library
LIBGTVIEW_C_SRC:=$(wildcard src/libgtview/*.c)
LIBGTVIEW_C_OBJ:=$(LIBGTVIEW_C_SRC:%.c=obj/%.o)
LIBGTVIEW_C_DEP:=$(LIBGTVIEW_C_SRC:%.c=obj/%.d)
LIBGTVIEW_LIBDEP=-lcairo -lgtext -lgtcore -lbz2 -lz

# the GenomeTools Lua library
LIBGTLUA_C_SRC:=$(wildcard src/libgtlua/*.c)
LIBGTLUA_C_OBJ:=$(LIBGTLUA_C_SRC:%.c=obj/%.o)
LIBGTLUA_C_DEP:=$(LIBGTLUA_C_SRC:%.c=obj/%.d)

TOOLS_SRC:=$(wildcard src/tools/*.c)
TOOLS_OBJ:=$(TOOLS_SRC:%.c=obj/%.o)
TOOLS_DEP:=$(TOOLS_SRC:%.c=obj/%.d)

LIBAGG_SRC:=$(wildcard src/external/agg-2.4/src/*.cpp src/external/agg-2.4/src/ctrl/*.cpp)
LIBAGG_OBJ:=$(LIBAGG_SRC:%.cpp=obj/%.o)
LIBAGG_DEP:=$(LIBAGG_SRC:%.cpp=obj/%.d)

EXPAT_DIR:=src/external/expat-2.0.1/lib
LIBEXPAT_SRC:=$(EXPAT_DIR)/xmlparse.c $(EXPAT_DIR)/xmlrole.c \
              $(EXPAT_DIR)/xmltok.c
LIBEXPAT_OBJ:=$(LIBEXPAT_SRC:%.c=obj/%.o)
LIBEXPAT_DEP:=$(LIBEXPAT_SRC:%.c=obj/%.d)

LUA_DIR:=src/external/lua-5.1.3/src
LIBLUA_SRC=$(LUA_DIR)/lapi.c $(LUA_DIR)/lcode.c $(LUA_DIR)/ldebug.c \
           $(LUA_DIR)/ldo.c $(LUA_DIR)/ldump.c $(LUA_DIR)/lfunc.c \
           $(LUA_DIR)/lgc.c $(LUA_DIR)/llex.c $(LUA_DIR)/lmem.c \
           $(LUA_DIR)/lobject.c $(LUA_DIR)/lopcodes.c $(LUA_DIR)/lparser.c \
           $(LUA_DIR)/lstate.c $(LUA_DIR)/lstring.c $(LUA_DIR)/ltable.c \
           $(LUA_DIR)/ltm.c $(LUA_DIR)/lundump.c $(LUA_DIR)/lvm.c \
           $(LUA_DIR)/lzio.c $(LUA_DIR)/lauxlib.c $(LUA_DIR)/lbaselib.c \
           $(LUA_DIR)/ldblib.c $(LUA_DIR)/liolib.c $(LUA_DIR)/lmathlib.c \
           $(LUA_DIR)/loslib.c $(LUA_DIR)/ltablib.c $(LUA_DIR)/lstrlib.c \
           $(LUA_DIR)/loadlib.c $(LUA_DIR)/linit.c \
           src/external/md5-1.1.2/src/md5.c\
           src/external/md5-1.1.2/src/md5lib.c\
           src/external/md5-1.1.2/src/des56.c\
           src/external/md5-1.1.2/src/ldes56.c\
           src/external/luafilesystem-1.4.1/src/lfs.c\
           src/external/lpeg-0.7/lpeg.c
LIBLUA_OBJ:=$(LIBLUA_SRC:%.c=obj/%.o)
LIBLUA_DEP:=$(LIBLUA_SRC:%.c=obj/%.d)

LUAMAIN_SRC:=$(LUA_DIR)/lua.c
LUAMAIN_OBJ:=$(LUAMAIN_SRC:%.c=obj/%.o)
LUAMAIN_DEP:=$(LUAMAIN_SRC:%.c=obj/%.d)

PNG_DIR:=src/external/libpng-1.2.18
LIBPNG_SRC:=$(PNG_DIR)/png.c $(PNG_DIR)/pngset.c $(PNG_DIR)/pngget.c \
            $(PNG_DIR)/pngrutil.c $(PNG_DIR)/pngtrans.c $(PNG_DIR)/pngwutil.c \
            $(PNG_DIR)/pngread.c $(PNG_DIR)/pngrio.c $(PNG_DIR)/pngwio.c \
            $(PNG_DIR)/pngwrite.c $(PNG_DIR)/pngrtran.c $(PNG_DIR)/pngwtran.c \
            $(PNG_DIR)/pngmem.c $(PNG_DIR)/pngerror.c $(PNG_DIR)/pngpread.c
LIBPNG_OBJ:=$(LIBPNG_SRC:%.c=obj/%.o)
LIBPNG_DEP:=$(LIBPNG_SRC:%.c=obj/%.d)

TECLA_DIR:=src/external/libtecla-1.6.1
LIBTECLA_SRC:=$(TECLA_DIR)/chrqueue.c $(TECLA_DIR)/cplfile.c \
              $(TECLA_DIR)/cplmatch.c $(TECLA_DIR)/direader.c \
              $(TECLA_DIR)/errmsg.c $(TECLA_DIR)/expand.c \
              $(TECLA_DIR)/freelist.c $(TECLA_DIR)/getline.c \
              $(TECLA_DIR)/hash.c $(TECLA_DIR)/history.c \
              $(TECLA_DIR)/homedir.c $(TECLA_DIR)/ioutil.c \
              $(TECLA_DIR)/keytab.c $(TECLA_DIR)/pathutil.c \
              $(TECLA_DIR)/pcache.c $(TECLA_DIR)/stringrp.c \
              $(TECLA_DIR)/strngmem.c $(TECLA_DIR)/version.c
LIBTECLA_OBJ:=$(LIBTECLA_SRC:%.c=obj/%.o)
LIBTECLA_DEP:=$(LIBTECLA_SRC:%.c=obj/%.d)

RNV_DIR:=src/external/rnv-1.7.8
LIBRNV_SRC:=$(RNV_DIR)/rn.c $(RNV_DIR)/rnc.c $(RNV_DIR)/rnd.c $(RNV_DIR)/rnl.c \
            $(RNV_DIR)/rnv.c $(RNV_DIR)/rnx.c $(RNV_DIR)/drv.c \
            $(RNV_DIR)/ary.c $(RNV_DIR)/xsd.c $(RNV_DIR)/xsd_tm.c \
            $(RNV_DIR)/dxl.c $(RNV_DIR)/dsl.c $(RNV_DIR)/sc.c $(RNV_DIR)/u.c \
            $(RNV_DIR)/ht.c $(RNV_DIR)/er.c $(RNV_DIR)/xmlc.c $(RNV_DIR)/s.c \
            $(RNV_DIR)/m.c $(RNV_DIR)/rx.c
LIBRNV_OBJ:=$(LIBRNV_SRC:%.c=obj/%.o)
LIBRNV_DEP:=$(LIBRNV_SRC:%.c=obj/%.d)

RNVMAIN_SRC:=$(RNV_DIR)/xcl.c
RNVMAIN_OBJ:=$(RNVMAIN_SRC:%.c=obj/%.o)
RNVMAIN_DEP:=$(RNVMAIN_SRC:%.c=obj/%.d)

BZ2_DIR:=src/external/bzip2-1.0.5
LIBBZ2_SRC:=$(BZ2_DIR)/blocksort.c $(BZ2_DIR)/huffman.c $(BZ2_DIR)/crctable.c \
            $(BZ2_DIR)/randtable.c $(BZ2_DIR)/compress.c \
            $(BZ2_DIR)/decompress.c $(BZ2_DIR)/bzlib.c
LIBBZ2_OBJ:=$(LIBBZ2_SRC:%.c=obj/%.o)
LIBBZ2_DEP:=$(LIBBZ2_SRC:%.c=obj/%.d)

ZLIB_DIR:=src/external/zlib-1.2.3
ZLIB_SRC:=$(ZLIB_DIR)/adler32.c $(ZLIB_DIR)/compress.c $(ZLIB_DIR)/crc32.c \
          $(ZLIB_DIR)/gzio.c $(ZLIB_DIR)/uncompr.c $(ZLIB_DIR)/deflate.c \
          $(ZLIB_DIR)/trees.c $(ZLIB_DIR)/zutil.c $(ZLIB_DIR)/inflate.c \
          $(ZLIB_DIR)/infback.c $(ZLIB_DIR)/inftrees.c $(ZLIB_DIR)/inffast.c
ZLIB_OBJ:=$(ZLIB_SRC:%.c=obj/%.o)
ZLIB_DEP:=$(ZLIB_SRC:%.c=obj/%.d)

# the objects which are included into the single GenomeThreader shared library
GTSHAREDLIB_OBJ:=$(LIBGTCORE_OBJ) $(LIBGTEXT_C_OBJ) $(LIBLUA_OBJ)
GTSHAREDLIB_LIBDEP:=$(LIBGTCORE_LIBDEP)

SERVER=gordon@genometools.org
WWWBASEDIR=/var/www/servers

# process arguments
ifeq ($(assert),no)
  EXP_CPPFLAGS += -DNDEBUG
endif

ifeq ($(cov),yes)
  export CCACHE_DISABLE # ccache cannot handle coverage objects
  GT_CFLAGS += -fprofile-arcs -ftest-coverage
  GT_LDFLAGS += -fprofile-arcs -ftest-coverage
  STEST_FLAGS += -gcov
  opt=no
  # hacks to link shared libs with cov=yes
  ifeq ($(SYSTEM),Linux)
    GT_LDFLAGS += -fprofile-arcs -ftest-coverage
  endif
  ifeq ($(MACHINE),amd64)
    GT_LDFLAGS += -fPIC
  endif
endif

ifneq ($(opt),no)
  GT_CFLAGS += -Os
  GT_CXXFLAGS += -Os
endif

ifeq ($(prof),yes)
  GT_CFLAGS += -pg
  GT_LDFLAGS += -pg
endif

ifeq ($(static),yes)
  GT_LDFLAGS += -static
endif

ifeq ($(curl),yes)
  EXP_CPPFLAGS += -DCURLDEF
  GT_CPPFLAGS += -I/usr/include/curl -I/usr/local/include/curl
  EXP_LDLIBS += -lcurl
endif

ifneq ($(curses),no)
  GTLIBS := $(GTLIBS) lib/libtecla.a
  EXP_CPPFLAGS += -DCURSES
  EXP_LDLIBS += -lncurses
endif

ifdef gttestdata
  STEST_FLAGS += -gttestdata $(gttestdata)
endif

ifeq ($(memcheck),yes)
  STEST_FLAGS += -memcheck
endif

# system specific stuff (concerning 64bit compilation)
ifeq ($(64bit),yes)
  ifneq ($(MACHINE),x86_64)
    m64=yes
  endif
else
  ifeq ($(MACHINE),x86_64)
    m32=yes
  endif
endif

ifeq ($(64bit),yes)
  BIT=64bit
else
  BIT=32bit
endif

ifeq ($(m32),yes)
  GT_CFLAGS += -m32
  GT_LDFLAGS += -m32
endif

ifeq ($(m64),yes)
  GT_CFLAGS += -m64
  GT_LDFLAGS += -m64
endif

ifeq ($(libgtview),yes)
  GTLIBS := $(GTLIBS) lib/libgtview.a
  GTSHAREDLIB_OBJ := $(GTSHAREDLIB_OBJ) $(LIBGTVIEW_C_OBJ)
  GTSHAREDLIB_LIBDEP:= $(GTSHAREDLIB_LIBDEP) -lcairo
  EXP_CPPFLAGS += -DLIBGTVIEW
  GT_CPPFLAGS += -I/usr/include/cairo -I/usr/local/include/cairo
  EXP_LDLIBS:=-lcairo $(EXP_LDLIBS)
  STEST_FLAGS += -libgtview
else
  OVERRIDELIBS += lib/libz.a # using own zlib together with cairo doesn't work
endif

# set prefix for install target
prefix ?= /usr/local

all: $(GTLIBS) lib/libgt$(SHARED_OBJ_NAME_EXT) bin/skproto bin/gt \
     bin/example bin/lua bin/rnv

lib/libexpat.a: $(LIBEXPAT_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBEXPAT_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libagg.a: $(LIBAGG_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBAGG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libbz2.a: $(LIBBZ2_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBBZ2_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libz.a: $(ZLIB_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(ZLIB_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtcore.a: obj/gt_config.h  $(LIBGTCORE_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTCORE_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgt$(SHARED_OBJ_NAME_EXT): obj/gt_config.h $(GTSHAREDLIB_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $(SHARED) $(GTSHAREDLIB_OBJ) \
	-o $@ $(GTSHAREDLIB_LIBDEP)

lib/libgtext.a: $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtexercise.a: $(LIBGTEXERCISE_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTEXERCISE_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtmgth.a: $(LIBGTMGTH_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTMGTH_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtmatch.a: $(LIBGTMATCH_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTMATCH_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtltr.a: $(LIBGTLTR_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTLTR_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif


lib/libgtview.a: $(LIBGTVIEW_C_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTVIEW_C_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtlua.a: $(LIBGTLUA_C_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTLUA_C_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libpng.a: $(LIBPNG_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBPNG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libtecla.a: $(LIBTECLA_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBTECLA_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/librnv.a: $(LIBRNV_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBRNV_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

define PROGRAM_template
$(1): $(2)
	@echo "[link $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) $$(EXP_LDFLAGS) $$(GT_LDFLAGS) $$(filter-out $$(OVERRIDELIBS),$$^) \
	  $$(filter-out $$(patsubst lib%.a,-l%,$$(notdir $$(OVERRIDELIBS))),\
	  $$(EXP_LDLIBS)) $$(OVERRIDELIBS) -o $$@

$(1)_static: $(2)
	@echo "[link $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) $$(EXP_LDFLAGS) $$(GT_LDFLAGS) $$(filter-out $$(OVERRIDELIBS),$$^) \
	  $$(filter-out $$(patsubst lib%.a,-l%,$$(notdir $$(OVERRIDELIBS))),\
	  $$(EXP_LDLIBS)) $$(OVERRIDELIBS) -static -o $$@
endef

$(eval $(call PROGRAM_template, bin/skproto, $(SKPROTO_OBJ) lib/libgtcore.a \
                                             $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/gt, $(GTMAIN_OBJ) $(TOOLS_OBJ) $(GTLIBS) \
                                        $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/example, $(EXAMPLE_OBJ) $(GTLIBS) \
                                             $(OVERRIDELIBS)))

bin/lua: $(LUAMAIN_OBJ) $(LIBLUA_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $^ -lm -o $@

bin/rnv: $(RNVMAIN_OBJ) lib/librnv.a lib/libexpat.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $^ -o $@

obj/gt_config.h:
	@echo '[create $@]'
	@test -d $(@D) || mkdir -p $(@D)
	@(echo '#define GT_BUILT $(BUILDSTAMP)' ;\
	echo '#define GT_CC "'`$(CC) --version | head -n 1`\" ;\
	echo '#define GT_CFLAGS "$(EXP_CFLAGS) $(GT_CFLAGS)"' ;\
	echo '$(EXP_CPPFLAGS) $(GT_CPPFLAGS)' | \
	sed -e 's/\([^\]\)"/\1\\"/g' -e 's/^"/\\"/g' -e 's/$$/"/' \
	    -e 's/^/#define GT_CPPFLAGS "/'; \
	  echo '#define GT_VERSION "'`cat VERSION`\" ) > $@

bitpackstringop_Dependencies=src/libgtcore/bitpackstringop.template \
	 src/libgtcore/bitpackstringvectorreadop.gen \
	 src/libgtcore/bitpackstringvectorwriteop.gen \
	 scripts/template2c.pl

src/libgtcore/bitpackstringop8.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $<

src/libgtcore/checkbitpackstring8.c: \
 src/libgtcore/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $<

src/libgtcore/bitpackstringop16.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $<

src/libgtcore/checkbitpackstring16.c: \
 src/libgtcore/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $<

src/libgtcore/bitpackstringop32.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $<

src/libgtcore/checkbitpackstring32.c: \
 src/libgtcore/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $<

src/libgtcore/bitpackstringop64.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $<

src/libgtcore/checkbitpackstring64.c: \
 src/libgtcore/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $<

src/libgtcore/checkbitpackstring-int.c: \
 src/libgtcore/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl '-int' $<

define COMPILE_template
$(1): $(2)
	@echo "[compile $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) -c $$< -o $$@ $$(EXP_CPPFLAGS) $$(GT_CPPFLAGS) $$(EXP_CFLAGS) \
	  $$(GT_CFLAGS) $(3)
	@$$(CC) -c $$< -o $$(@:.o=.d) $$(EXP_CPPFLAGS) $$(GT_CPPFLAGS) \
        $(3) -MM -MP -MT $$@
endef

$(eval $(call COMPILE_template, obj/%.o, %.c))

obj/%.o: %.cxx
	@echo "[compile $@]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CXXFLAGS) $(GT_CXXFLAGS)
	@$(CXX) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GT_CPPFLAGS) -MM -MP \
	  -MT $@

obj/%.o: %.cpp
	@echo "[compile $@]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CXXFLAGS) $(GT_CXXFLAGS)
	@$(CXX) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GT_CPPFLAGS) -MM -MP \
	  -MT $@

obj/src/libgtcore/versionfunc.o: obj/gt_config.h

# read dependencies
-include $(GTMAIN_DEP) \
         $(EXAMPLE_DEP) \
         $(SKPROTO_DEP) \
	 $(LIBGTCORE_DEP) \
	 $(LIBGTEXT_C_DEP) \
	 $(LIBGTEXT_CXX_DEP) \
	 $(LIBGTEXERCISE_DEP) \
	 $(LIBGTMGTH_DEP) \
	 $(LIBGTMATCH_DEP) \
	 $(LIBGTLTR_DEP) \
	 $(LIBGTLUA_C_DEP) \
	 $(TOOLS_DEP) \
	 $(LIBAGG_DEP) \
	 $(LIBEXPAT_DEP) \
	 $(LIBLUA_DEP) \
         $(LUAMAIN_DEP) \
	 $(LIBPNG_DEP) \
	 $(LIBTECLA_DEP) \
	 $(LIBRNV_DEP) \
	 $(RNVMAIN_DEP) \
	 $(LIBBZ2_DEP) \
	 $(ZLIB_DEP)

ifeq ($(libgtview),yes)
-include $(LIBGTVIEW_C_DEP) $(LIBGTVIEW_CXX_DEP)
endif

.SUFFIXES:
.PHONY: dist srcdist release gt install docs installwww splint test clean cleanup

VERSION:="`cat $(CURDIR)/VERSION`"
SYSTEMNAME:="$(SYSTEM)_$(MACHINE)"
GTDISTBASENAME:="gt-$(VERSION)-$(SYSTEMNAME)-${BIT}"
DISTDIR:="$(CURDIR)/dist/$(SYSTEMNAME)"
GTDISTDIR:="$(DISTDIR)/$(GTDISTBASENAME)"

dist: all
	@echo "[build distribution]"
	@rm -rf $(GTDISTDIR)
	@rm -rf $(DISTDIR)/$(GTDISTBASENAME).tar.gz
	@mkdir -p $(GTDISTDIR)/bin
	@cp $(CURDIR)/doc/dist_readme.txt $(GTDISTDIR)/README
	@cp $(CURDIR)/LICENSE $(GTDISTDIR)
	@cp $(CURDIR)/CONTRIBUTORS $(GTDISTDIR)
	@cp $(CURDIR)/CHANGELOG $(GTDISTDIR)
	@cp $(CURDIR)/bin/gt $(GTDISTDIR)/bin
	@strip $(GTDISTDIR)/bin/gt
	@cp -r $(CURDIR)/gtdata $(GTDISTDIR)
	@cd $(DISTDIR) && tar cf $(GTDISTBASENAME).tar $(GTDISTBASENAME)
	@cd $(DISTDIR) && gzip -f -9 $(GTDISTBASENAME).tar
	@echo "$(DISTDIR)/$(GTDISTBASENAME).tar.gz"

srcdist:
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
        gzip -9 > genometools-`cat VERSION`.tar.gz

release:
	git tag "v`cat VERSION`"
	git archive --format=tar --prefix="genometools-`cat VERSION`"/ HEAD | \
	gzip -9 > genometools-`cat VERSION`.tar.gz
	scp "genometools-`cat VERSION`.tar.gz" $(SERVER):$(WWWBASEDIR)/genometools.org/htdocs/pub
	git push --tags origin master

docs: bin/gt
	bin/gt gtscripts/gtdoc.lua -html $(CURDIR) \
        > www/genometools.org/htdocs/docs.html

installwww:
# install genometools.org website
	rsync -rv www/genometools.org/ $(SERVER):$(WWWBASEDIR)/genometools.org

gt: bin/gt

install: all
	test -d $(prefix)/bin || mkdir -p $(prefix)/bin
	cp bin/gt $(prefix)/bin
	cp -r gtdata $(prefix)/bin
	test -d $(prefix)/include/libgtcore \
	  || mkdir -p $(prefix)/include/libgtcore
	cp src/gtcore.h obj/gt_config.h $(prefix)/include
	cp src/libgtcore/*.h $(prefix)/include/libgtcore
	test -d $(prefix)/include/libgtext \
          || mkdir -p $(prefix)/include/libgtext
	cp src/gtext.h $(prefix)/include
	cp src/libgtext/*.h $(prefix)/include/libgtext
	cp src/genometools.h $(prefix)/include
	test -d $(prefix)/include/libgtmatch \
	  || mkdir -p $(prefix)/include/libgtmatch
	cp src/gtmatch.h $(prefix)/include
	cp src/libgtmatch/*.h src/libgtmatch/*.pr $(prefix)/include/libgtmatch
	test -d $(prefix)/lib || mkdir -p $(prefix)/lib
	cp lib/libgtcore.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtcore.a
endif
	cp lib/libgtext.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtext.a
endif
	cp lib/libgtmatch.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtmatch.a
endif
	@echo '[build config script $(@F)]'
	sed -e 's!@CC@!$(CC)!' -e 's!@CFLAGS@!$(EXP_CFLAGS)!' \
	  -e 's!@CPPFLAGS@!$(subst ",\\",-I"$(prefix)/include" $(EXP_CPPFLAGS))!' \
	  -e 's!@CXX@!$(CXX)!' -e 's!@CXXFLAGS@!$(EXP_CXXFLAGS)!' \
	  -e 's!@LDFLAGS@!-L$(prefix)/lib $(EXP_LDFLAGS)!' \
	  -e 's!@LIBS@!$(EXP_LDLIBS)!' -e "s!@VERSION@!`cat VERSION`!" \
	  -e 's!@BUILDSTAMP@!$(BUILDSTAMP)!' \
	  -e 's!@SYSTEM@!$(SYSTEM)!' <src/genometools-config.in \
	  >$(prefix)/bin/genometools-config
	chmod 755 $(prefix)/bin/genometools-config

cflags:
	@echo ${GT_CFLAGS}

splint: obj/gt_config.h
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) \
	$(CURDIR)/src/*.c \
        $(CURDIR)/src/libgtcore/*.c \
        $(CURDIR)/src/libgtext/*.c \
        $(CURDIR)/src/tools/*.c

EISFILES=${shell ls ${CURDIR}/src/libgtmatch/*.c | grep eis-}\
         ${CURDIR}/src/libgtmatch/sfx-opt.c\
         ${CURDIR}/src/libgtmatch/sfx-run.c\
         ${CURDIR}/src/libgtmatch/encseq-specialsrank.c

SKTOOLS=${shell grep -l Kurtz src/tools/*.c}

spgt:${addprefix obj/,${notdir ${subst .c,.splint,\
	             ${filter-out ${EISFILES},${wildcard ${CURDIR}/src/libgtmatch/*.c}}\
                     ${wildcard ${CURDIR}/src/libgtltr/*.c}\
                                ${SKTOOLS}}}}

scgt:
	src_check src/libgtmatch/*
	src_check src/libgtltr/*
	src_check src/tools/*

splintclean:
	find obj -name '*.splint' | xargs rm -f

obj/%.splint: ${CURDIR}/src/libgtmatch/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/tools/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/libgtltr/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.prepro: ${CURDIR}/src/libgtmatch/%.c
	@echo "[generate $@]"
	$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS) -E -g3
	indent $@

RUBY:=ruby

test: all
	GT_MEM_BOOKKEEPING=on bin/gt -test
	cd testsuite && env -i GT_MEM_BOOKKEEPING=on MAKE=$(MAKE) PATH=$(PATH) \
          CCACHE_DISABLE=yes HOME=$(HOME) \
          $(RUBY) -I. testsuite.rb \
          -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin -cur $(CURDIR) \
          -gtruby $(CURDIR)/gtruby $(STEST_FLAGS)

clean:
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests

cleanup: clean
	rm -rf lib bin
