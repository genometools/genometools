#
# Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
            -I$(CURDIR)/src/external/lua-5.1.4/src \
            -I$(CURDIR)/src/external/luafilesystem-1.4.1/src \
            -I$(CURDIR)/src/external/lpeg-0.9 \
            -I$(CURDIR)/src/external/expat-2.0.1/lib \
            -I$(CURDIR)/src/external/bzip2-1.0.5 \
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
GT_CFLAGS:=-g -Wall -Wunused-parameter -pipe -fPIC -Wpointer-arith
GT_CFLAGS_NO_WERROR:=$(GT_CFLAGS) -w
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
	        test -d /usr/X11R6/lib && echo /usr/X11R6/lib ; \
                test -d /opt/local/lib && echo /opt/local/lib),-L$(dir))
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

# libraries for which we build replacements (that also appear in dependencies)
EXP_LDLIBS+=-lz -lbz2
OVERRIDELIBS:=lib/libbz2.a

# compiled executables
GTMAIN_SRC:=src/gt.c src/gtr.c src/gtt.c src/interactive.c
GTMAIN_OBJ:=$(GTMAIN_SRC:%.c=obj/%.o)
GTMAIN_DEP:=$(GTMAIN_SRC:%.c=obj/%.d)

EXAMPLES_SRC:=src/example.c
EXAMPLES_DEP:=$(EXAMPLES_SRC:%.c=obj/%.d)

SKPROTO_SRC:=src/skproto.c src/tools/gt_skproto.c
SKPROTO_OBJ:=$(SKPROTO_SRC:%.c=obj/%.o)
SKPROTO_DEP:=$(SKPROTO_SRC:%.c=obj/%.d)

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
           src/external/luafilesystem-1.4.1/src/lfs.c\
           src/external/lpeg-0.9/lpeg.c
LIBLUA_OBJ:=$(LIBLUA_SRC:%.c=obj/%.o)
LIBLUA_DEP:=$(LIBLUA_SRC:%.c=obj/%.d)

LUAMAIN_SRC:=src/external/lua-5.1.4/etc/all.c
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

HMMER_BASE:=src/external/hmmer-2.3.2
HMMER_DIR:=$(HMMER_BASE)/src
SQUID_DIR:=$(HMMER_BASE)/squid
HMMER_SRC:=$(HMMER_DIR)/alphabet.c \
           $(HMMER_DIR)/core_algorithms.c $(HMMER_DIR)/debug.c \
           $(HMMER_DIR)/display.c $(HMMER_DIR)/emit.c \
           $(HMMER_DIR)/fast_algorithms.c $(HMMER_DIR)/histogram.c \
           $(HMMER_DIR)/hmmio.c $(HMMER_DIR)/mathsupport.c $(HMMER_DIR)/masks.c \
           $(HMMER_DIR)/misc.c $(HMMER_DIR)/plan7.c $(HMMER_DIR)/plan9.c \
           $(HMMER_DIR)/postprob.c $(HMMER_DIR)/prior.c $(HMMER_DIR)/pvm.c \
           $(HMMER_DIR)/threads.c $(HMMER_DIR)/tophits.c $(HMMER_DIR)/trace.c \
           $(SQUID_DIR)/a2m.c $(SQUID_DIR)/aligneval.c $(SQUID_DIR)/alignio.c \
           $(SQUID_DIR)/clustal.c $(SQUID_DIR)/cluster.c $(SQUID_DIR)/dayhoff.c \
           $(SQUID_DIR)/eps.c $(SQUID_DIR)/file.c $(SQUID_DIR)/getopt.c \
           $(SQUID_DIR)/gki.c $(SQUID_DIR)/hsregex.c \
           $(SQUID_DIR)/iupac.c $(SQUID_DIR)/msa.c $(SQUID_DIR)/msf.c \
           $(SQUID_DIR)/phylip.c $(SQUID_DIR)/revcomp.c $(SQUID_DIR)/rk.c \
           $(SQUID_DIR)/selex.c $(SQUID_DIR)/seqencode.c $(SQUID_DIR)/shuffle.c \
           $(SQUID_DIR)/sqerror.c $(SQUID_DIR)/sqio.c $(SQUID_DIR)/squidcore.c \
           $(SQUID_DIR)/sre_ctype.c $(SQUID_DIR)/sre_math.c $(SQUID_DIR)/sre_random.c \
           $(SQUID_DIR)/sre_string.c $(SQUID_DIR)/ssi.c $(SQUID_DIR)/stack.c \
           $(SQUID_DIR)/stockholm.c $(SQUID_DIR)/stopwatch.c \
           $(SQUID_DIR)/translate.c $(SQUID_DIR)/types.c \
           $(SQUID_DIR)/vectorops.c $(SQUID_DIR)/weight.c
HMMER_OBJ:=$(HMMER_SRC:%.c=obj/%.o)
HMMER_DEP:=$(HMMER_SRC:%.c=obj/%.d)

ZLIB_DIR:=src/external/zlib-1.2.3
ZLIB_SRC:=$(ZLIB_DIR)/adler32.c $(ZLIB_DIR)/compress.c $(ZLIB_DIR)/crc32.c \
          $(ZLIB_DIR)/gzio.c $(ZLIB_DIR)/uncompr.c $(ZLIB_DIR)/deflate.c \
          $(ZLIB_DIR)/trees.c $(ZLIB_DIR)/zutil.c $(ZLIB_DIR)/inflate.c \
          $(ZLIB_DIR)/infback.c $(ZLIB_DIR)/inftrees.c $(ZLIB_DIR)/inffast.c
ZLIB_OBJ:=$(ZLIB_SRC:%.c=obj/%.o)
ZLIB_DEP:=$(ZLIB_SRC:%.c=obj/%.d)

# the objects which are included into the single GenomeTools shared library
GTSHAREDLIB_LIBDEP:=-lbz2 -lz

SERVER=gordon@genometools.org
WWWBASEDIR=/var/www/servers

# process arguments
ifeq ($(assert),no)
  EXP_CPPFLAGS += -DNDEBUG
endif

ifneq ($(errorcheck),no)
  GT_CFLAGS += -Werror
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
  GT_CFLAGS += -O3
  GT_CXXFLAGS += -O3
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
  EXP_CPPFLAGS += -DCURSES
  EXP_LDLIBS += -lncurses
  GTLIBS := lib/libtecla.a
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

ifneq ($(sharedlib),no)
  SHARED_LIBGENOMETOOLS := lib/libgenometools$(SHARED_OBJ_NAME_EXT)
endif

LIBGENOMETOOLS_DIRS:= src/core \
                      src/extended \
                      src/gtlua \
                      src/match \
                      src/gth \
                      src/mgth \
                      src/ltr

ifeq ($(with-hmmer),yes)
  LIBGENOMETOOLS_DIRS := src/external/hmmer-2.3.2  $(LIBGENOMETOOLS_DIRS)
  EXP_CPPFLAGS += -DHAVE_HMMER
  GT_CPPFLAGS +=  -I$(CURDIR)/$(HMMER_DIR) -I$(CURDIR)/$(SQUID_DIR)
  OVERRIDELIBS += lib/libhmmer.a
  STEST_FLAGS += -hmmer
endif

ifneq ($(cairo),no)
  GTSHAREDLIB_LIBDEP:= $(GTSHAREDLIB_LIBDEP) -lcairo -lfontconfig
  GT_CPPFLAGS += -I/usr/include/cairo -I/usr/local/include/cairo \
                 -I/opt/local/include/cairo -I/usr/include/fontconfig \
                 -I/usr/local/include/fontconfig -I/opt/local/include/fontconfig
  EXP_LDLIBS:=-lcairo -lfontconfig $(EXP_LDLIBS)
  ANNOTATIONSKETCH_EXAMPLES := bin/examples/sketch_constructed \
                               bin/examples/sketch_parsed_with_ctrack \
                               bin/examples/sketch_parsed_with_ordering \
                               bin/examples/sketch_parsed
  ANNOTATIONSKETCH_MANUAL := doc/manuals/annotationsketch.pdf
  LIBGENOMETOOLS_DIRS:=$(LIBGENOMETOOLS_DIRS) src/annotationsketch
else
  OVERRIDELIBS += lib/libz.a # using own zlib together with cairo doesn't work
  EXP_CPPFLAGS += -DWITHOUT_CAIRO
  STEST_FLAGS += -nocairo
endif

ifeq ($(threads),yes)
  EXP_CPPFLAGS += -DGT_THREADS_ENABLED
  EXP_LDLIBS += -lpthread
  GTSHAREDLIB_LIBDEP += -lpthread
  HMMERTHREADS_DEF := #define HMMER_THREADS
else
  HMMERTHREADS_DEF := #undef HMMER_THREADS
endif

# the GenomeTools library
LIBGENOMETOOLS_PRESRC:=$(foreach DIR,$(LIBGENOMETOOLS_DIRS),$(wildcard $(DIR)/*.c))
ifeq ($(amalgamation),yes)
  LIBGENOMETOOLS_SRC:=obj/amalgamation.c
else
  LIBGENOMETOOLS_SRC:=$(LIBGENOMETOOLS_PRESRC)
endif
LIBGENOMETOOLS_OBJ:=$(LIBGENOMETOOLS_SRC:%.c=obj/%.o) \
                    $(LIBLUA_OBJ) \
                    $(LIBEXPAT_OBJ)
LIBGENOMETOOLS_DEP:=$(LIBGENOMETOOLS_SRC:%.c=obj/%.d) \
                    $(LIBLUA_DEP) \
                    $(LIBEXPAT_DEP)

ifeq ($(with-hmmer),yes)
  LIBGENOMETOOLS_OBJ := lib/libhmmer.a $(LIBGENOMETOOLS_OBJ) lib/libhmmer.a
endif

# set prefix for install target
prefix ?= /usr/local

# allow to set patch program
patch ?= patch

all: lib/libgenometools.a $(SHARED_LIBGENOMETOOLS) \
     bin/skproto bin/gt bin/lua bin/rnv bin/examples/custom_stream \
     bin/examples/gff3validator bin/examples/noop $(ANNOTATIONSKETCH_EXAMPLES)

lib/libexpat.a: $(LIBEXPAT_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBEXPAT_OBJ)
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

lib/libhmmer.a: hmmer_get $(HMMER_DIR)/config.h.in $(HMMER_DIR)/config.h $(SQUID_DIR)/squidconf.h $(SQUID_DIR)/squid.h $(HMMER_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(HMMER_OBJ)
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

lib/libgenometools.a: obj/gt_config.h  $(LIBGENOMETOOLS_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGENOMETOOLS_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgenometools$(SHARED_OBJ_NAME_EXT): obj/gt_config.h \
                                          $(LIBGENOMETOOLS_OBJ) lib/libbz2.a \
                                          lib/libz.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $(SHARED) $(LIBGENOMETOOLS_OBJ) \
	-o $@ $(GTSHAREDLIB_LIBDEP)

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

$(eval $(call PROGRAM_template, bin/skproto, $(SKPROTO_OBJ) \
                                             lib/libgenometools.a\
                                             $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/gt, $(GTMAIN_OBJ) $(TOOLS_OBJ) \
                                        lib/libgenometools.a \
                                        $(GTLIBS) \
                                        $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/custom_stream, \
                                obj/src/examples/custom_stream.o \
                                lib/libgenometools.a \
                                $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/gff3validator, \
                                obj/src/examples/gff3validator.o \
                                lib/libgenometools.a \
                                $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/noop, \
                                obj/src/examples/noop.o \
                                lib/libgenometools.a \
                                $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_constructed, \
                                obj/src/examples/sketch_constructed.o \
                                lib/libgenometools.a $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed, \
                                obj/src/examples/sketch_parsed.o \
                                lib/libgenometools.a $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed_with_ctrack, \
                                obj/src/examples/sketch_parsed_with_ctrack.o \
                                lib/libgenometools.a $(OVERRIDELIBS)))

$(eval $(call PROGRAM_template, bin/examples/sketch_parsed_with_ordering, \
                                obj/src/examples/sketch_parsed_with_ordering.o \
                                lib/libgenometools.a $(OVERRIDELIBS)))


bin/lua: $(LUAMAIN_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $^ -lm -o $@

bin/rnv: $(RNVMAIN_OBJ) lib/librnv.a lib/libexpat.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GT_LDFLAGS) $^ -o $@

$(SQUID_DIR)/squidconf.h:
	@echo '[create $(@F)]'
	@scripts/generate_hmmer_squidconf_h $(SQUID_DIR)/squidconf.h.in > $@

$(SQUID_DIR)/squid.h:
	@echo '[create $(@F)]'
	@scripts/generate_hmmer_squid_h  $(SQUID_DIR)/squid.h.in > $@

$(HMMER_DIR)/config.h: $(HMMER_DIR)/config.h.in
	@echo '[create $(@F)]'
	@sed  -e 's/#undef PACKAGE_VERSION/#define PACKAGE_VERSION "2.3.2"/'\
	      -e 's/#undef HMMER_THREADS/$(HMMERTHREADS_DEF)/'\
	      -e 's/#undef PACKAGE_NAME/#define PACKAGE_NAME "HMMER"/'\
	      -e 's/#undef PACKAGE_TARNAME/#define PACKAGE_TARNAME "hmmer"/'\
	      -e 's/#undef PACKAGE_STRING/#define PACKAGE_STRING "HMMER 2.3.2"/'\
	      -e 's/#undef PACKAGE_DATE/#define PACKAGE_DATE "Oct 2003"/'\
	      -e 's/#undef PACKAGE_COPYRIGHT/#define PACKAGE_COPYRIGHT "Copyright (C) 1992-2003 HHMI\/Washington University School of Medicine"/'\
	      -e 's/#undef PACKAGE_LICENSE/#define PACKAGE_LICENSE "Freely distributed under the GNU General Public License (GPL)"/' $(HMMER_DIR)/config.h.in  > $@

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

$(SQUID_DIR)/%.c: $(HMMER_DIR)/config.h
	@true

$(HMMER_DIR)/%.c: $(HMMER_DIR)/config.h
	@true

obj/$(SQUID_DIR)/%.o: $(SQUID_DIR)/%.c
	@echo "[compile $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS_NO_WERROR) $(3)
	@$(CC) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(3) -MM -MP \
	  -MT $@

obj/$(HMMER_DIR)/%.o: $(HMMER_DIR)/%.c
	@echo "[compile $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) -c $< -o $@ $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(EXP_CFLAGS) \
	  $(GT_CFLAGS_NO_WERROR) $(3)
	@$(CC) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GT_CPPFLAGS) $(3) -MM -MP \
	  -MT $@

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

obj/src/core/versionfunc.o: obj/gt_config.h

# read dependencies
-include $(GTMAIN_DEP) \
         $(EXAMPLES_DEP) \
         $(SKPROTO_DEP) \
	 $(TOOLS_DEP) \
         $(LUAMAIN_DEP) \
	 $(LIBTECLA_DEP) \
	 $(LIBRNV_DEP) \
	 $(RNVMAIN_DEP) \
	 $(LIBBZ2_DEP) \
         $(HMMER_DEP) \
	 $(ZLIB_DEP) \
         $(LIBGENOMETOOLS_DEP) \
         obj/src/examples/custom_stream.d \
         obj/src/examples/gff3validator.d \
         obj/src/examples/noop.d \
         obj/src/examples/sketch_constructed.d \
         obj/src/examples/sketch_parsed.d

ifeq ($(libannotationsketch),yes)
-include $(LIBANNOTATIONSKETCH_C_DEP) $(LIBANNOTATIONSKETCH_CXX_DEP)
endif

.PRECIOUS: $(HMMER_DIR)/%.c $(SQUID_DIR)/%.c $(HMMER_DIR)/%.h.in $(SQUID_DIR)/%.h.in
.SUFFIXES:
.PHONY: dist srcdist release gt install docs manuals installwww push \
        splint test clean cleanup hmmer_get

VERSION:="`cat $(CURDIR)/VERSION`"
SYSTEMNAME:="$(SYSTEM)_$(MACHINE)"
GTDISTBASENAME:="gt-$(VERSION)-$(SYSTEMNAME)-${BIT}"
DISTDIR:="$(CURDIR)/dist/$(SYSTEMNAME)"
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
	git push --tags github master

docs: bin/gt bin/examples/sketch_parsed bin/examples/sketch_constructed
	bin/gt gtscripts/gtdoc.lua -html $(CURDIR) \
        > www/genometools.org/htdocs/libgenometools.html
	bin/gt gtscripts/gtdoc.lua -lua -html $(CURDIR) \
        > www/genometools.org/htdocs/docs.html
	bin/examples/sketch_parsed gtdata/sketch/default.style \
          www/genometools.org/htdocs/images/parsed.png \
          testdata/eden.gff3
	bin/examples/sketch_parsed \
	  www/genometools.org/htdocs/annotationsketch/callbacks.style \
	  www/genometools.org/htdocs/images/callbacks.png \
          www/genometools.org/htdocs/annotationsketch/callback_examples_with_score.gff3
	bin/examples/sketch_constructed gtdata/sketch/default.style \
	  www/genometools.org/htdocs/images/constructed.png
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

doc/manuals/annotationsketch.pdf: docs
	$(MAKE) -C $(CURDIR)/doc/manuals annotationsketch

manuals: $(ANNOTATIONSKETCH_MANUAL)
	$(MAKE) -C $(CURDIR)/doc/manuals

installwww:
# install genometools.org website
	rsync -rv www/genometools.org/ $(SERVER):$(WWWBASEDIR)/genometools.org

patch:
	if test -f $(CURDIR)/src/patches/$(SYSTEMNAME).patch; then $(patch) -p1 < $(CURDIR)/src/patches/$(SYSTEMNAME).patch; fi

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

cflags:
	@echo ${GT_CFLAGS}

splint: obj/gt_config.h
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) \
	$(CURDIR)/src/*.c \
        $(CURDIR)/src/core/*.c \
        $(CURDIR)/src/extended/*.c \
        $(CURDIR)/src/tools/*.c

EISFILES=${shell ls ${CURDIR}/src/match/*.c | grep eis-}\
         ${CURDIR}/src/match/sfx-opt.c\
         ${CURDIR}/src/match/sfx-run.c\
         ${CURDIR}/src/match/encseq-specialsrank.c

# Ask Sascha to splint the following files which are currently included.

LTRdigestfiles=${CURDIR}/src/ltr/gt_ltrdigest.c\
               ${CURDIR}/src/ltr/gt_ltrdigest.h\
               ${CURDIR}/src/ltr/ltr_visitor.c\
               ${CURDIR}/src/ltr/ltr_visitor.h\
               ${CURDIR}/src/ltr/ltrdigest_def.h\
               ${CURDIR}/src/ltr/ltrdigest_stream.c\
               ${CURDIR}/src/ltr/ltrdigest_stream.h\
               ${CURDIR}/src/ltr/ltrelement.c\
               ${CURDIR}/src/ltr/ltrelement.h\
               ${CURDIR}/src/ltr/ltrfileout_stream.c\
               ${CURDIR}/src/ltr/ltrfileout_stream.h\
               ${CURDIR}/src/ltr/pbs.c\
               ${CURDIR}/src/ltr/pbs.h\
               ${CURDIR}/src/ltr/pdom.c\
               ${CURDIR}/src/ltr/pdom.h\
               ${CURDIR}/src/ltr/ppt.c\
               ${CURDIR}/src/ltr/ppt.h

SKTOOLS=${shell grep -l Kurtz src/tools/*.c}

MOVEDENCSEQFILES = ${CURDIR}/src/core/mapspec-gen.c\
                   ${CURDIR}/src/core/readmode.c\
                   ${CURDIR}/src/core/encodedsequence.c

ALLSPLINT=${addprefix obj/,${notdir ${subst .c,.splint,\
	             ${filter-out ${EISFILES},${wildcard ${CURDIR}/src/match/*.c}}\
                     ${filter-out ${LTRdigestfiles},${wildcard ${CURDIR}/src/ltr/*.c}}\
                                ${SKTOOLS} ${MOVEDENCSEQFILES}}}}\
     obj/redblack.splint


spgt:${ALLSPLINT}

scgt:
	src_check src/match/*
	src_check src/ltr/*
	src_check src/tools/*

splintclean:
	find obj -name '*.splint' | xargs rm -f

obj/%.splint: ${CURDIR}/src/match/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/tools/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/ltr/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/core/%.c
	@echo "splint $<"
	@splint -DBIGSEQPOS -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.splint: ${CURDIR}/src/extended/%.c
	@echo "splint $<"
	@splint -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
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

RUBY:=ruby

test: all
	GT_MEM_BOOKKEEPING=on bin/gt -test
	cd testsuite && env -i GT_MEM_BOOKKEEPING=on MAKE=$(MAKE) PATH=$(PATH) \
          CCACHE_DISABLE=yes HOME=$(HOME) \
          $(RUBY) -I. testsuite.rb \
          -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin -cur $(CURDIR) \
          -gtruby $(CURDIR)/gtruby $(STEST_FLAGS)

clean:
	rm -rf $(HMMER_DIR)/config.h $(SQUID_DIR)/squidconf.h $(SQUID_DIR)/squid.h
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests

gtkviewer:
	$(CC) -o bin/examples/gtkviewer $(GT_CPPFLAGS) $(GT_LDFLAGS) `pkg-config --cflags --libs gtk+-2.0` -lgenometools src/examples/gtkviewer.c

cleanup: clean
	rm -rf lib bin

hmmer_get:
	@echo "[check for HMMER source distribution]"
	@test ! -f "$(HMMER_DIR)/config.h.in" && \
	  cd src/external && \
	  echo "[retrieve HMMER source distribution]" && \
	  wget -q ftp://selab.janelia.org/pub/software/hmmer/CURRENT/hmmer-2.3.2.tar.gz && \
	  tar -xzf hmmer-2.3.2.tar.gz && \
	  rm hmmer-2.3.2.tar.gz || true
