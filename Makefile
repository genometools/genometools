#
# Copyright (c) 2006-2013 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2006-2013 Center for Bioinformatics, University of Hamburg
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
            -I$(CURDIR)/src/external/zlib-1.2.7 \
            -I$(CURDIR)/src/external/md5-1.1.2/src \
            -I$(CURDIR)/src/external/lua-5.1.5/src \
            -I$(CURDIR)/src/external/luafilesystem-1.5.0/src \
            -I$(CURDIR)/src/external/lpeg-0.10.2 \
            -I$(CURDIR)/src/external/expat-2.0.1/lib \
            -I$(CURDIR)/src/external/bzip2-1.0.6 \
            -I$(CURDIR)/src/external/libtecla-1.6.1 \
            -I$(CURDIR)/src/external/samtools-0.1.18 \
            -I$(CURDIR)/src/external/sqlite-3.8.0.1

ifeq ($(shell pkg-config --version > /dev/null 2> /dev/null; echo $$?),0)
  HAS_PKGCONFIG:=yes
else
  HAS_PKGCONFIG:=no
endif

ifneq ($(cairo),no)
  ifeq ($(HAS_PKGCONFIG),yes)
	INCLUDEOPT+=$(shell pkg-config --cflags-only-I pango) \
	            $(shell pkg-config --cflags-only-I cairo) \
	            $(shell pkg-config --cflags-only-I pangocairo) \
	            $(shell pkg-config --cflags-only-I glib-2.0)
  endif
endif

ifneq ($(fpic),no)
  FPIC:=-fPIC
endif

# these variables are exported by the configuration script
ifndef CC
  CC:=gcc
endif
ifndef CXX
  CXX:=g++
endif
EXP_CFLAGS:=$(CFLAGS)
EXP_LDFLAGS:=$(LDFLAGS)
EXP_CXXFLAGS:=$(CXXFLAGS)
EXP_CPPFLAGS:=$(CPPFLAGS)
EXP_LDLIBS:=$(LIBS) -lm
# ...while those starting with GT_ are for internal purposes only
GT_CFLAGS:=-g -Wall -Wunused-parameter -pipe $(FPIC) -Wpointer-arith
# expat needs -DHAVE_MEMMOVE
# tecla needs -DHAVE_CURSES_H -DHAVE_TERM_H -DUSE_TERMINFO
# zlib needs -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN
EXT_FLAGS:= -DHAVE_MEMMOVE -DHAVE_CURSES_H -DHAVE_TERM_H -DUSE_TERMINFO \
            -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN
EXP_CPPFLAGS+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 $(EXT_FLAGS)
GT_CPPFLAGS:=$(INCLUDEOPT)
GT_CXXFLAGS:=-g -pipe
GT_LDFLAGS:=-Llib
STEST_FLAGS:=
EXP_LDFLAGS+=$(foreach dir, \
	$(shell test -d /usr/local/lib && echo /usr/local/lib ; \
	        test -d /usr/X11R6/lib && echo /usr/X11R6/lib ; \
                test -d /opt/local/lib && echo /opt/local/lib),-L$(dir))
BUILDSTAMP:=$(shell date +'"%Y-%m-%d %H:%M:%S"')

# try to set RANLIB automatically
SYSTEM:=$(shell uname -s)
MACHINE:=$(shell uname -m)
ifeq ($(SYSTEM),Darwin)
  RANLIB:=ranlib
  NO_STATIC_LINKING:=defined
  SHARED:=-dynamiclib -undefined dynamic_lookup
  SHARED_OBJ_NAME_EXT:=.dylib
  ifeq ($(universal),yes)
    MACHINE:="Universal_Binary"
    GT_CFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
    GT_LDFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
  endif
  ifeq ($(ppc),yes)
    MACHINE:="Power_Macintosh"
    GT_CFLAGS+=-arch ppc -arch_errors_fatal
    GT_LDFLAGS+=-arch ppc -arch_errors_fatal
  endif
else
  SHARED_OBJ_NAME_EXT:=.so
  SHARED:=-shared
endif

# compiled executables
GTMAIN_SRC:=src/gt.c src/gtr.c src/gtt.c src/interactive.c
GTMAIN_OBJ:=$(GTMAIN_SRC:%.c=obj/%.o)
GTMAIN_DEP:=$(GTMAIN_SRC:%.c=obj/%.d)

EXAMPLES_SRC:=src/example.c
EXAMPLES_DEP:=$(EXAMPLES_SRC:%.c=obj/%.d)

TOOLS_SRC:=$(wildcard src/tools/*.c)
TOOLS_OBJ:=$(TOOLS_SRC:%.c=obj/%.o)
TOOLS_DEP:=$(TOOLS_SRC:%.c=obj/%.d)

EXPAT_DIR:=src/external/expat-2.0.1/lib
LIBEXPAT_SRC:=$(EXPAT_DIR)/xmlparse.c $(EXPAT_DIR)/xmlrole.c \
              $(EXPAT_DIR)/xmltok.c
LIBEXPAT_OBJ:=$(LIBEXPAT_SRC:%.c=obj/%.o)
LIBEXPAT_DEP:=$(LIBEXPAT_SRC:%.c=obj/%.d)

LIBLUA_SRC=src/lualib.c\
           src/external/md5-1.1.2/src/md5.c\
           src/external/md5-1.1.2/src/md5lib.c\
           src/external/md5-1.1.2/src/des56.c\
           src/external/md5-1.1.2/src/ldes56.c\
           src/external/luafilesystem-1.5.0/src/lfs.c\
           src/external/lpeg-0.10.2/lpeg.c
LIBLUA_OBJ:=$(LIBLUA_SRC:%.c=obj/%.o)
LIBLUA_DEP:=$(LIBLUA_SRC:%.c=obj/%.d)

LUAMAIN_SRC:=src/external/lua-5.1.5/etc/all.c
LUAMAIN_OBJ:=$(LUAMAIN_SRC:%.c=obj/%.o)
LUAMAIN_DEP:=$(LUAMAIN_SRC:%.c=obj/%.d)

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

BZ2_DIR:=src/external/bzip2-1.0.6
LIBBZ2_SRC:=$(BZ2_DIR)/blocksort.c $(BZ2_DIR)/huffman.c $(BZ2_DIR)/crctable.c \
            $(BZ2_DIR)/randtable.c $(BZ2_DIR)/compress.c \
            $(BZ2_DIR)/decompress.c $(BZ2_DIR)/bzlib.c
LIBBZ2_OBJ:=$(LIBBZ2_SRC:%.c=obj/%.o)
LIBBZ2_DEP:=$(LIBBZ2_SRC:%.c=obj/%.d)

SQLITE3_DIR:=src/external/sqlite-3.8.0.1
SQLITE3_SRC:=$(SQLITE3_DIR)/sqlite3.c
SQLITE3_OBJ:=$(SQLITE3_SRC:%.c=obj/%.o)
SQLITE3_DEP:=$(SQLITE3_SRC:%.c=obj/%.d)

ZLIB_DIR:=src/external/zlib-1.2.7
ZLIB_SRC:=$(ZLIB_DIR)/adler32.c $(ZLIB_DIR)/compress.c $(ZLIB_DIR)/crc32.c \
          $(ZLIB_DIR)/gzclose.c $(ZLIB_DIR)/gzlib.c $(ZLIB_DIR)/gzread.c \
          $(ZLIB_DIR)/gzwrite.c $(ZLIB_DIR)/uncompr.c $(ZLIB_DIR)/deflate.c \
          $(ZLIB_DIR)/trees.c $(ZLIB_DIR)/zutil.c $(ZLIB_DIR)/inflate.c \
          $(ZLIB_DIR)/infback.c $(ZLIB_DIR)/inftrees.c $(ZLIB_DIR)/inffast.c
ZLIB_OBJ:=$(ZLIB_SRC:%.c=obj/%.o)
ZLIB_DEP:=$(ZLIB_SRC:%.c=obj/%.d)

SAMTOOLS_DIR:=src/external/samtools-0.1.18
SAMTOOLS_SRC:=$(SAMTOOLS_DIR)/bgzf.c \
              $(SAMTOOLS_DIR)/kstring.c \
              $(SAMTOOLS_DIR)/bam_aux.c \
              $(SAMTOOLS_DIR)/bam.c \
              $(SAMTOOLS_DIR)/bam_import.c \
              $(SAMTOOLS_DIR)/sam.c \
              $(SAMTOOLS_DIR)/bam_index.c \
              $(SAMTOOLS_DIR)/bam_pileup.c \
              $(SAMTOOLS_DIR)/bam_lpileup.c \
              $(SAMTOOLS_DIR)/bam_md.c \
              $(SAMTOOLS_DIR)/razf.c \
              $(SAMTOOLS_DIR)/faidx.c \
              $(SAMTOOLS_DIR)/bedidx.c \
              $(SAMTOOLS_DIR)/knetfile.c \
              $(SAMTOOLS_DIR)/bam_sort.c \
              $(SAMTOOLS_DIR)/sam_header.c \
              $(SAMTOOLS_DIR)/bam_reheader.c \
              $(SAMTOOLS_DIR)/kprobaln.c \
              $(SAMTOOLS_DIR)/bam_cat.c
SAMTOOLS_OBJ:=$(SAMTOOLS_SRC:%.c=obj/%.o)
SAMTOOLS_DEP:=$(SAMTOOLS_SRC:%.c=obj/%.d)

# add necessary shared lib dependencies then not building them ourselves
ifeq ($(useshared),yes)
  DEPLIBS:=-lbz2 -lz -lexpat -llua5.1-lpeg -llua5.1 -llua5.1-md5 \
           -llua5.1-filesystem -llua5.1-des56 -lbam
else
  DEPLIBS:=
endif
EXP_LDLIBS += $(DEPLIBS)
GTSHAREDLIB_LIBDEP += $(DEPLIBS)

SERVER=satta@genometools.org
WWWBASEDIR=/var/www/servers

# process arguments
ifeq ($(assert),no)
  EXP_CPPFLAGS += -DNDEBUG
endif

ifeq ($(shell hmmscan -h > /dev/null 2> /dev/null; echo $$?),0)
  STEST_FLAGS += -hmmer
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
    GT_LDFLAGS += $(FPIC)
  endif
endif

ifneq ($(opt),no)
  ifeq ($(findstring clang,$(CC)),clang)
    ifeq ($(lto),yes)
      # clang/LLVM supports link-time optimization with -O4
      # Note that on Linux this usually requires additional configuration of
      # LLVM! Consult http://llvm.org/docs/GoldPlugin.html if problems arise.
      GT_CFLAGS += -O4
      GT_CXXFLAGS += -O4
    else
      GT_CFLAGS += -O3
      GT_CXXFLAGS += -O3
    endif
  else
    GT_CFLAGS += -O3
    GT_CXXFLAGS += -O3
  endif
endif

ifeq ($(prof),yes)
  GT_CFLAGS += -pg
  GT_LDFLAGS += -pg
endif

ifeq ($(curl),yes)
  EXP_CPPFLAGS += -DCURLDEF
  GT_CPPFLAGS += -I/usr/include/curl -I/usr/local/include/curl
  EXP_LDLIBS += -lcurl
endif

ifneq ($(curses),no)
  EXP_CPPFLAGS += -DCURSES
  ifeq ($(findstring CYGWIN,$(SYSTEM)),CYGWIN)
    EXP_CPPFLAGS += $(shell pkg-config --cflags-only-I ncurses)
  endif
  EXP_LDLIBS += -lncurses
  GTLIBS := lib/libtecla.a
endif

ifdef testthreads
  STEST_FLAGS += -threads $(testthreads)
endif

ifdef gttestdata
  STEST_FLAGS += -gttestdata $(gttestdata)
endif

ifeq ($(memcheck),yes)
  STEST_FLAGS += -memcheck
endif

# system specific stuff (concerning 64bit compilation)
ifeq ($(SYSTEM),Darwin)
  ifeq ($(64bit),yes)
    m64=yes
  else
    m32=yes
  endif
endif

ifeq ($(SYSTEM),Darwin)
  EXT_FLAGS += -DLUA_DL_DYLD
else
  EXT_FLAGS += -DLUA_DL_DLOPEN
  ifneq ($(SYSTEM),FreeBSD)
  ifneq ($(SYSTEM),Windows)
    LUA_LDLIB:=-ldl
    EXP_LDLIBS += -ldl
  endif
  endif
endif

ifeq ($(64bit),yes)
  ifneq ($(MACHINE),x86_64)
    m64=yes
  endif
  BIT=64bit
  SPLINTD:=-D_LP64
else
  ifeq ($(MACHINE),x86_64)
    m32=yes
  endif
  BIT=32bit
endif

ifeq ($(m32),yes)
  GT_CFLAGS += -m32
  GT_LDFLAGS += -m32
  SQLITE_CFLAGS += -m32
endif

ifeq ($(m64),yes)
  ifneq ($(MACHINE),ia64)
    ifneq ($(MACHINE),alpha)
      GT_CFLAGS += -m64
      GT_LDFLAGS += -m64
      SQLITE_CFLAGS += -m64
    endif
  endif
endif

ifeq ($(findstring clang,$(CC)),clang)
  # do not complain about unnecessary options
  GT_CFLAGS += -Qunused-arguments -Wno-parentheses
  GT_CPPFLAGS += -Qunused-arguments -Wno-parentheses
endif

ifneq ($(sharedlib),no)
  SHARED_LIBGENOMETOOLS := lib/libgenometools$(SHARED_OBJ_NAME_EXT)
endif

LIBGENOMETOOLS_DIRS:= src/core \
                      src/extended \
                      src/gtlua \
                      src/match \
                      src/gth \
                      src/ltr \
                      src/mgth

ifneq ($(cairo),no)
  ifeq ($(SYSTEM),Darwin)
    EXP_LDFLAGS:=
  endif
  ifeq ($(HAS_PKGCONFIG),yes)
    GTSHAREDLIB_LIBDEP:= $(GTSHAREDLIB_LIBDEP) \
                         $(shell pkg-config --libs pango) \
                         $(shell pkg-config --libs cairo) \
                         $(shell pkg-config --libs pangocairo)
    EXP_LDLIBS:=$(EXP_LDLIBS) \
                $(shell pkg-config --libs pango) \
                $(shell pkg-config --libs cairo) \
                $(shell pkg-config --libs pangocairo)
  endif
  ANNOTATIONSKETCH_EXAMPLES := bin/examples/sketch_constructed \
                               bin/examples/sketch_parsed_with_ctrack \
                               bin/examples/sketch_parsed_with_ordering \
                               bin/examples/sketch_parsed
  ANNOTATIONSKETCH_MANUAL := doc/manuals/annotationsketch.pdf
  LIBGENOMETOOLS_DIRS:=$(LIBGENOMETOOLS_DIRS) src/annotationsketch
  SKETCH_EXTRA_BINARIES := bin/examples/sketch_parsed \
                           bin/examples/sketch_constructed
else
  EXP_CPPFLAGS += -DWITHOUT_CAIRO
  STEST_FLAGS += -nocairo
  CAIRO_FILTER_OUT:=src/gtlua/annotationsketch_lua.c \
                    src/gtlua/canvas_lua.c \
                    src/gtlua/diagram_lua.c \
                    src/gtlua/image_info_lua.c \
                    src/gtlua/layout_lua.c
endif

ifeq ($(threads),yes)
  EXP_CPPFLAGS += -DGT_THREADS_ENABLED
  EXP_LDLIBS += -lpthread
  GTSHAREDLIB_LIBDEP += -lpthread
endif

ifneq ($(with-sqlite),no)
  ifeq ($(useshared),yes)
    EXP_LDLIBS += -lsqlite3
    GTSHAREDLIB_LIBDEP += -lsqlite3
  else
    LIBGENOMETOOLS_DIRS := $(SQLITE3_DIR) $(LIBGENOMETOOLS_DIRS)
    GT_CPPFLAGS +=  -I$(CURDIR)/$(SQLITE3_DIR)
    OVERRIDELIBS += lib/libsqlite.a
  endif
  EXP_CPPFLAGS += -DHAVE_SQLITE
  ifeq ($(threads),yes)
    EXP_CPPFLAGS += -DSQLITE_THREADSAFE=1
  else
    EXP_CPPFLAGS += -DSQLITE_THREADSAFE=0
  endif
  EXP_LDLIBS += -lpthread
  ifneq ($(SYSTEM),FreeBSD)
    EXP_LDLIBS += -ldl
  endif
else
  SQLITE_FILTER_OUT:=src/extended/rdb_sqlite.c
endif

ifeq ($(with-mysql),yes)
  GTSHAREDLIB_LIBDEP:= $(GTSHAREDLIB_LIBDEP) -lmysqlclient
  EXP_CPPFLAGS += -DHAVE_MYSQL
  EXP_LDLIBS += -lmysqlclient
else
  MYSQL_FILTER_OUT:=src/extended/rdb_mysql.c
  ifeq ($(with-sqlite),no)
    STEST_FLAGS += -nordb
  endif
endif

# the GenomeTools library
LIBGENOMETOOLS_PRESRC:=$(foreach DIR,$(LIBGENOMETOOLS_DIRS),$(wildcard $(DIR)/*.c))
# remove AnnotationSketch-only bindings
LIBGENOMETOOLS_PRESRC:=$(filter-out $(CAIRO_FILTER_OUT),\
                         $(LIBGENOMETOOLS_PRESRC))
# remove MySQL-only files
LIBGENOMETOOLS_PRESRC:=$(filter-out $(MYSQL_FILTER_OUT),\
                         $(LIBGENOMETOOLS_PRESRC))
# remove SQLite-only files
LIBGENOMETOOLS_PRESRC:=$(filter-out $(SQLITE_FILTER_OUT),\
                         $(LIBGENOMETOOLS_PRESRC))

ifeq ($(amalgamation),yes)
  LIBGENOMETOOLS_SRC:=obj/amalgamation.c
else
  LIBGENOMETOOLS_SRC:=$(LIBGENOMETOOLS_PRESRC)
endif
LIBGENOMETOOLS_OBJ:=$(LIBGENOMETOOLS_SRC:%.c=obj/%.o)
LIBGENOMETOOLS_DEP:=$(LIBGENOMETOOLS_SRC:%.c=obj/%.d)

ifneq ($(useshared),yes)
  LIBGENOMETOOLS_OBJ += $(LIBLUA_OBJ) \
                        $(SAMTOOLS_OBJ) \
                        $(LIBEXPAT_OBJ) \
                        $(LIBBZ2_OBJ) \
                        $(ZLIB_OBJ)
  LIBGENOMETOOLS_DEP += $(LIBLUA_DEP) \
                        $(SAMTOOLS_DEP) \
                        $(LIBEXPAT_DEP) \
                        $(LIBBZ2_DEP) \
                        $(ZLIB_DEP)
endif

ifneq ($(with-sqlite),no)
  ifneq ($(useshared),yes)
    LIBGENOMETOOLS_OBJ := $(LIBGENOMETOOLS_OBJ) lib/libsqlite.a
  endif
endif

GT_CFLAGS_NO_WERROR:=$(GT_CFLAGS) -w

ifneq ($(errorcheck),no)
  GT_CFLAGS += -Werror
endif

# set prefix for install target
prefix ?= /usr/local

ifdef DESTDIR
  prefix:=$(DESTDIR)$(prefix)
endif

# allow to set patch program
patch ?= patch

ifneq ($(useshared),yes)
  ADDITIONAL_BINARIES:=bin/lua
else
  ADDITIONAL_BINARIES:=
endif

all: lib/libgenometools.a $(SHARED_LIBGENOMETOOLS) \
     bin/gt bin/examples/custom_stream \
     bin/examples/gff3sort bin/examples/gff3validator bin/examples/noop \
     $(ANNOTATIONSKETCH_EXAMPLES) \
     $(ADDITIONAL_BINARIES)

lib/libexpat.a: $(LIBEXPAT_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(LIBEXPAT_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libsqlite.a: $(SQLITE3_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(SQLITE3_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libbz2.a: $(LIBBZ2_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(LIBBZ2_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libz.a: $(ZLIB_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(ZLIB_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgenometools.a: obj/gt_config.h  $(LIBGENOMETOOLS_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) r $@ $(LIBGENOMETOOLS_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

ifneq ($(useshared),yes)
  ADDITIONAL_ZLIBS:=lib/libbz2.a lib/libz.a
  ADDITIONAL_SO_DEPS:=$(ADDITIONAL_ZLIBS) \
                      lib/libexpat.a
endif

lib/libgenometools$(SHARED_OBJ_NAME_EXT): obj/gt_config.h \
                                          $(LIBGENOMETOOLS_OBJ) \
                                          $(ADDITIONAL_SO_DEPS) \
                                          $(ADDITIONAL_ZLIBS)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $(ADDITIONAL_SO_DEPS) $(SHARED) $(LIBGENOMETOOLS_OBJ) \
	-o $@ $(GTSHAREDLIB_LIBDEP)

lib/libtecla.a: $(LIBTECLA_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(LIBTECLA_OBJ)
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
endef

$(eval $(call PROGRAM_template, bin/gt, $(GTMAIN_OBJ) $(TOOLS_OBJ) \
                                        lib/libgenometools.a \
                                        $(GTLIBS) \
                                        $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/custom_stream, \
                                obj/src/examples/custom_stream.o \
                                lib/libgenometools.a\
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/gff3sort, \
                                obj/src/examples/gff3sort.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/gff3validator, \
                                obj/src/examples/gff3validator.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/noop, \
                                obj/src/examples/noop.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_constructed, \
                                obj/src/examples/sketch_constructed.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed, \
                                obj/src/examples/sketch_parsed.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed_with_ctrack, \
                                obj/src/examples/sketch_parsed_with_ctrack.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed_with_ordering, \
                                obj/src/examples/sketch_parsed_with_ordering.o \
                                lib/libgenometools.a \
                                $(ADDITIONAL_ZLIBS)))

bin/lua: $(LUAMAIN_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $^ -lm $(LUA_LDLIB) -o $@

obj/gt_config.h: VERSION
	@echo '[create $@]'
	@test -d $(@D) || mkdir -p $(@D)
	@(echo '#ifndef GT_CONFIG_H' ;\
	echo '#define GT_CONFIG_H' ;\
        echo '#define GT_BUILT $(BUILDSTAMP)' ;\
	echo '#define GT_CC "'`$(CC) --version | head -n 1`\" ;\
	echo '#define GT_CFLAGS "$(EXP_CFLAGS) $(GT_CFLAGS)"' ;\
	echo '$(EXP_CPPFLAGS) $(GT_CPPFLAGS)' | \
	sed -e 's/\([^\]\)"/\1\\"/g' -e 's/^"/\\"/g' -e 's/$$/"/' \
	    -e 's/^/#define GT_CPPFLAGS "/'; \
	  echo '#define GT_VERSION "'`cat VERSION`\" ) > $@
	@cat VERSION | \
          sed 's/\([0-9]*\)\.[0-9]*\.[0-9]*/#define GT_MAJOR_VERSION \1/' >> $@
	@cat VERSION | \
          sed 's/[0-9]*\.\([0-9]*\)\.[0-9]*/#define GT_MINOR_VERSION \1/' >> $@
	@cat VERSION | \
          sed 's/[0-9]*\.[0-9]*\.\([0-9]*\)/#define GT_MICRO_VERSION \1/' >> $@
	@echo '#endif' >> $@

obj/amalgamation.c: $(LIBGENOMETOOLS_PRESRC)
	@echo '[create $@]'
	@test -d $(@D) || mkdir -p $(@D)
	@scripts/create_amalgamation $(LIBGENOMETOOLS_PRESRC) > $@

bitpackstringop_Dependencies=src/core/bitpackstringop.template \
	 src/core/bitpackstringvectorreadop.gen \
	 src/core/bitpackstringvectorwriteop.gen \
	 scripts/template2c.pl

src/core/bitpackstringop8.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $<

src/core/checkbitpackstring8.c: \
 src/core/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $<

src/core/bitpackstringop16.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $<

src/core/checkbitpackstring16.c: \
 src/core/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $<

src/core/bitpackstringop32.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $<

src/core/checkbitpackstring32.c: \
 src/core/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $<

src/core/bitpackstringop64.c: $(bitpackstringop_Dependencies)
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $<

src/core/checkbitpackstring64.c: \
 src/core/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $<

src/core/checkbitpackstring-int.c: \
 src/core/checkbitpackstring.template scripts/template2c.pl
	@echo '[rebuild $@]'
	@scripts/template2c.pl '-int' $<

# SQLite needs special attention
obj/$(SQLITE3_DIR)/%.o: $(SQLITE3_DIR)/%.c
	@echo "[compile $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) -c $< -o $@ -DHAVE_MALLOC_USABLE_SIZE $(EXP_CPPFLAGS) \
	  $(GT_CPPFLAGS) $(EXP_CFLAGS) $(SQLITE_CFLAGS) -DSQLITE_ENABLE_UNLOCK_NOTIFY  $(3) $(FPIC)
	@$(CC) -c $< -o $(@:.o=.d) -DHAVE_MALLOC_USABLE_SIZE $(EXP_CPPFLAGS) \
	  $(GT_CPPFLAGS) $(3) -MM -MP -MT $@ $(FPIC)

define COMPILE_template
$(1): $(2)
	@echo "[compile $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) -c $$< -o $$@ $$(EXP_CPPFLAGS) $$(GT_CPPFLAGS) $$(EXP_CFLAGS) \
	  $$(GT_CFLAGS) $(3)
	@$$(CC) -c $$< -o $$(@:.o=.d) $$(EXP_CPPFLAGS) $$(GT_CPPFLAGS) \
        $(3) -MM -MP -MT $$@
endef

obj/$(SAMTOOLS_DIR)/%.o: $(SAMTOOLS_DIR)/%.c
	@echo "[compile $(notdir $@)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS_NO_WERROR)
	@$(CC) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GT_CPPFLAGS) -MM -MP \
	  -MT $@

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

obj/src/core/versionfunc.o: obj/gt_config.h

# read dependencies
-include $(GTMAIN_DEP) \
         $(EXAMPLES_DEP) \
         $(TOOLS_DEP) \
         $(LUAMAIN_DEP) \
         $(LIBTECLA_DEP) \
         $(LIBBZ2_DEP) \
         $(ZLIB_DEP) \
         $(SAMTOOLS_DEP) \
         $(LIBGENOMETOOLS_DEP) \
         obj/src/examples/custom_stream.d \
         obj/src/examples/gff3sort.d \
         obj/src/examples/gff3validator.d \
         obj/src/examples/noop.d \
         obj/src/examples/sketch_constructed.d \
         obj/src/examples/sketch_parsed.d

.PRECIOUS:
.SUFFIXES:
.PHONY: dist srcdist release gt install docs manuals installwww push \
        splint test clean cleanup

VERSION:="`cat $(CURDIR)/VERSION`"
SYSTEMNAME:="$(SYSTEM)_$(MACHINE)"
GTDISTBASENAME:="gt-$(VERSION)-$(SYSTEMNAME)-${BIT}"
DISTDIR:="$(CURDIR)/dist/$(SYSTEMNAME)"
SCRIPTSDIR:="$(CURDIR)/scripts"
GTDISTDIR:="$(DISTDIR)/$(GTDISTBASENAME)"

dist: all manuals
	@echo "[build distribution]"
	@rm -rf $(GTDISTDIR)
	@rm -rf $(DISTDIR)/$(GTDISTBASENAME).tar.gz
	@mkdir -p $(GTDISTDIR)/bin $(GTDISTDIR)/doc
	@cp $(CURDIR)/doc/dist_readme.txt $(GTDISTDIR)/README
	@cp $(CURDIR)/LICENSE $(GTDISTDIR)
	@cp $(CURDIR)/CONTRIBUTORS $(GTDISTDIR)
	@cp $(CURDIR)/CHANGELOG $(GTDISTDIR)
	@cp $(CURDIR)/bin/gt $(GTDISTDIR)/bin
	@strip $(GTDISTDIR)/bin/gt
	@cp $(CURDIR)/doc/manuals/*.pdf $(GTDISTDIR)/doc
	@cp -r $(CURDIR)/gtdata $(GTDISTDIR)
	@cp -r $(CURDIR)/gtpython $(GTDISTDIR)
	@cp -r $(CURDIR)/gtruby $(GTDISTDIR)
	@$(MAKE) prefix=$(GTDISTDIR) install
	@cd $(DISTDIR) && $(SCRIPTSDIR)/tar_root.sh $(GTDISTBASENAME)
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
	git push --tags github master

docs: bin/gt $(SKETCH_EXTRA_BINARIES)
	bin/gt gtscripts/gtdoc.lua -html $(CURDIR) \
        > www/genometools.org/htdocs/libgenometools.html
	bin/gt gtscripts/gtdoc.lua -lua -html $(CURDIR) \
        > www/genometools.org/htdocs/docs.html
ifdef SKETCH_EXTRA_BINARIES
	bin/examples/sketch_parsed gtdata/sketch/default.style \
          www/genometools.org/htdocs/images/parsed.png \
          testdata/eden.gff3
	bin/examples/sketch_parsed \
	  www/genometools.org/htdocs/annotationsketch/callbacks.style \
	  www/genometools.org/htdocs/images/callbacks.png \
          www/genometools.org/htdocs/annotationsketch/callback_examples_with_score.gff3
	bin/examples/sketch_constructed gtdata/sketch/default.style \
	  www/genometools.org/htdocs/images/constructed.png
endif
	sed -nf scripts/incl.sed \
	  www/genometools.org/htdocs/examples_tmpl.html | \
          sed 'N;N;s/\n//' > /tmp/tmp.sed.$$$$ && \
	sed -f /tmp/tmp.sed.$$$$ \
	  www/genometools.org/htdocs/examples_tmpl.html > \
	  www/genometools.org/htdocs/examples.html; rm -f /tmp/tmp.sed.$$$$
	bin/gt gtscripts/gtdoc.lua -tex $(CURDIR) \
	> doc/manuals/api_reference.tex
	bin/gt gtscripts/gtdoc.lua -lua -tex $(CURDIR) \
	> doc/manuals/gtscript_reference.tex
	$(MAKE) -C $(CURDIR)/doc/devguide
	bin/gt -list > /tmp/list.txt
	bin/gt -createman /tmp/gtmanpages
	asciidoc --backend=xhtml11 -a linkcss -a stylesdir=. \
	  -a stylesheet=style.css -a badges -a icons \
	  -f www/genometools.org/htdocs/tool_list.conf \
	  --out-file www/genometools.org/htdocs/tools.html /tmp/list.txt
	test -d www/genometools.org/htdocs/tools || mkdir -p www/genometools.org/htdocs/tools
	rm -f www/genometools.org/htdocs/tools/*
	scripts/create_html /tmp/gtmanpages www/genometools.org/htdocs/
	rm -rf /tmp/list.txt /tmp/gtmanpages

doc/manuals/annotationsketch.pdf: docs
	$(MAKE) -C $(CURDIR)/doc/manuals annotationsketch

manuals: $(ANNOTATIONSKETCH_MANUAL)
	$(MAKE) -C $(CURDIR)/doc/manuals

manpages: bin/gt
	bin/gt -createman /tmp/gtmanpages
	test -d $(CURDIR)/doc/manpages || mkdir -p $(CURDIR)/doc/manpages
	rm -f $(CURDIR)/doc/manpages/*
	scripts/create_manpages /tmp/gtmanpages $(CURDIR)/doc/manpages

installwww:
# install genometools.org website
	rsync -rv www/genometools.org/ $(SERVER):$(WWWBASEDIR)/genometools.org

push:
	git push origin master
	git push github master

gt: bin/gt

install: all
	test -d $(prefix)/bin || mkdir -p $(prefix)/bin
	cp bin/gt $(prefix)/bin
	cp -r gtdata $(prefix)/bin
	test -d $(prefix)/include/genometools/core \
	  || mkdir -p $(prefix)/include/genometools/core
	cp src/core/*_api.h $(prefix)/include/genometools/core
	test -d $(prefix)/include/genometools/extended \
          || mkdir -p $(prefix)/include/genometools/extended
	cp src/extended/*_api.h $(prefix)/include/genometools/extended
	test -d $(prefix)/include/genometools/annotationsketch \
          || mkdir -p $(prefix)/include/genometools/annotationsketch
	cp src/annotationsketch/*_api.h \
          $(prefix)/include/genometools/annotationsketch
	test -d $(prefix)/include/genometools/ltr \
          || mkdir -p $(prefix)/include/genometools/ltr
	cp src/ltr/*_api.h $(prefix)/include/genometools/ltr
	cp obj/gt_config.h $(prefix)/include/genometools
	cp src/genometools.h $(prefix)/include/genometools
	test -d $(prefix)/lib || mkdir -p $(prefix)/lib
	cp lib/libgenometools.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgenometools.a
endif
	cp lib/libgenometools$(SHARED_OBJ_NAME_EXT) $(prefix)/lib
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

installmanpages: manpages
	test ! -d $(CURDIR)/doc/manpages || \
	  ((test -d $(prefix)/share/man/man1 || mkdir -p $(prefix)/share/man/man1) \
	    && cp $(CURDIR)/doc/manpages/* $(prefix)/share/man/man1)

cflags:
	@echo ${GT_CFLAGS}

splint: obj/gt_config.h
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) \
	$(CURDIR)/src/*.c \
        $(CURDIR)/src/core/*.c \
        $(CURDIR)/src/extended/*.c \
        $(CURDIR)/src/tools/*.c

EISFILES=${shell ls ${CURDIR}/src/match/*.c | grep eis-}

SKTOOLS=${shell grep -l Kurtz src/tools/*.c}
SKCORE=${shell grep -l 'Stefan Kurtz' src/core/*.c}
DWCORE=${shell grep -l Willrodt src/core/*.c}
DWTOOLS=${shell grep -l Willrodt src/tools/*.c}
GGTOOLS=${shell grep -l Gonnella src/tools/*.c}
SKEXT=${shell grep -l 'Stefan Kurtz' src/extended/*.c}
DWEXT=${shell grep -l Willrodt src/extended/*.c}
OEEXT=${shell grep -l Eigenbrod src/extended/*.c}
GGEXT=${shell grep -l Gonnella src/extended/*.c}

ALLSPLINT=${addprefix obj/,${notdir ${subst .c,.splint,\
             ${filter-out ${EISFILES},${wildcard ${CURDIR}/src/match/*.c}}\
             ${wildcard ${CURDIR}/src/ltr/*.c}\
             ${SKTOOLS} ${SKCORE} ${SKEXT} \
						 ${DWTOOLS} ${DWCORE} ${DWEXT} \
             ${GGTOOLS} ${GGEXT} ${OEEXT} }}}

ALLSCANBUILD=${subst .splint,.sb, ${ALLSPLINT}}

spgt:${ALLSPLINT}

sbgt:${ALLSCANBUILD}

scgt:
	src_check src/core/*
	src_check src/match/*
	src_check src/ltr/*
	src_check src/extended/*
	src_check src/tools/*

splintclean:
	find obj -name '*.splint' | xargs rm -f

sbclean:
	find obj -name '*.sb' | xargs rm -f

obj/%.sb: src/match/%.c
	@echo "scan-build $<"
	@scan-build -analyze-headers --status-bugs --use-cc $(CC) $(CC) -c $< $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS) -o obj/${subst .c,.o,$<} > /dev/null
	@touch $@

obj/%.sb: src/ltr/%.c
	@echo "scan-build $<"
	@scan-build -analyze-headers --status-bugs --use-cc $(CC) $(CC) -c $< $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS) -o obj/${subst .c,.o,$<} > /dev/null
	@touch $@

obj/%.sb: src/tools/%.c
	@echo "scan-build $<"
	@scan-build -analyze-headers --status-bugs --use-cc $(CC) $(CC) -c $< $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS) -o obj/${subst .c,.o,$<} > /dev/null
	@touch $@

obj/%.sb: src/core/%.c
	@echo "scan-build $<"
	@scan-build -analyze-headers --status-bugs --use-cc $(CC) $(CC) -c $< $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS) -o obj/${subst .c,.o,$<} > /dev/null
	@touch $@

obj/%.sb: src/extended/%.c
	@echo "scan-build $<"
	@scan-build -analyze-headers --status-bugs --use-cc $(CC) $(CC) -c $< $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS) -o obj/${subst .c,.o,$<} > /dev/null
	@touch $@

obj/%.splint: ${CURDIR}/src/match/%.c
	@echo "splint $<"
	@splint $(SPLINTD) $(EXP_CPPFLAGS) $(INCLUDEOPT) -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/tools/%.c
	@echo "splint $<"
	@splint $(SPLINTD) $(EXP_CPPFLAGS) $(INCLUDEOPT) -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/ltr/%.c
	@echo "splint $<"
	@splint $(SPLINTD) $(EXP_CPPFLAGS) $(INCLUDEOPT) -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/core/%.c
	@echo "splint $<"
	@splint $(SPLINTD) $(EXP_CPPFLAGS) $(INCLUDEOPT) -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/extended/%.c
	@echo "splint $<"
	@splint $(SPLINTD) $(EXP_CPPFLAGS) $(INCLUDEOPT) -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@


DWHEADER=${filter-out %impl.h, ${shell find ${CURDIR} -name '*.h' | \
				   xargs grep -l Willrodt}}
ALLHEADER=${addprefix obj/,${notdir ${subst .h,.check,\
					${DWHEADER}}}}

headercheck:${ALLHEADER}

headerclean:
	find obj -name '*.check' | xargs rm -f

obj/%.check: ${CURDIR}/src/match/%.h
	@echo "check include $<"
	@src_check_header.rb $<
	@touch $@

obj/%.check: ${CURDIR}/src/tools/%.h
	@echo "check include $<"
	@src_check_header.rb $<
	@touch $@

obj/%.check: ${CURDIR}/src/ltr/%.h
	@echo "check include $<"
	@src_check_header.rb $<
	@touch $@

obj/%.check: ${CURDIR}/src/core/%.h
	@echo "check include $<"
	@src_check_header.rb $<
	@touch $@

obj/%.check: ${CURDIR}/src/extended/%.h
	@echo "check include $<"
	@src_check_header.rb $<
	@touch $@

obj/%.prepro: ${CURDIR}/src/match/%.c
	@echo "[generate $@]"
	$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS) -E -g3
	indent $@

obj/%.prepro: ${CURDIR}/src/tools/%.c
	@echo "[generate $@]"
	$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS) -E -g3
	indent $@

obj/%.prepro: ${CURDIR}/src/core/%.c
	@echo "[generate $@]"
	$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS) -E -g3
	indent $@

obj/%.prepro: ${CURDIR}/src/extended/%.c
	@echo "[generate $@]"
	$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) \
	  $(EXP_CFLAGS) $(GT_CFLAGS) -E -g3
	indent $@

RUBY:=ruby

ifndef testrange
  testrange:=0..5000
endif

test: all
	GT_MEM_BOOKKEEPING=on bin/gt -test
	cd testsuite && env -i GT_MEM_BOOKKEEPING=on MAKE=$(MAKE) PATH=$(PATH) \
          CCACHE_DISABLE=yes HOME=$(HOME) \
          $(RUBY) -I. testsuite.rb \
          -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin -cur $(CURDIR) \
          -gtruby $(CURDIR)/gtruby $(STEST_FLAGS) \
          -select $(testrange)

clean:
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests
	$(MAKE) -s -C $(CURDIR)/doc/devguide clean
	$(MAKE) -s -C $(CURDIR)/doc/manuals cleanup

cleangenerated:
	rm -f doc/manuals/api_reference.tex \
          doc/manuals/gtscript_reference.tex
	find doc . -name "*.toc" -delete
	rm -f www/genometools.org/htdocs/images/callbacks.png \
          www/genometools.org/htdocs/images/parsed.png    \
          www/genometools.org/htdocs/images/constructed.png \
          doc/manuals/annotationsketch.pdf
	rm -f www/genometools.org/htdocs/examples.html
	rm -rf doc/manpages

gtkviewer:
	@echo "[compile $(notdir $@)]"
	@$(CC) -o bin/examples/gtkviewer $(GT_CPPFLAGS) $(GT_LDFLAGS) \
  src/examples/gtkviewer.c  -lcairo `pkg-config --cflags --libs gtk+-2.0` \
  -lgenometools

cleanup: clean cleangenerated
	rm -rf lib bin
	rm -rf gtpython/build
