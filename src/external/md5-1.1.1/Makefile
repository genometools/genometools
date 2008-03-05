# $Id: Makefile,v 1.7 2007/10/11 00:02:56 carregal Exp $

CONFIG= ./config

include $(CONFIG)

ifeq "$(LUA_VERSION_NUM)" "500"
COMPAT_O= $(COMPAT_DIR)/compat-5.1.o
endif

MD5_OBJS= src/md5.o src/md5lib.o $(COMPAT_O)
MD5_LUAS= src/md5.lua
MD5_LIBNAME = core.so

DES56_OBJS= src/des56.o src/ldes56.o
DES56_LIBNAME= des56.so

all: src/$(MD5_LIBNAME) src/$(DES56_LIBNAME)

src/$(MD5_LIBNAME) : $(MD5_OBJS)
	export MACOSX_DEPLOYMENT_TARGET="10.3"; $(CC) $(CFLAGS) $(LIB_OPTION) -o src/$(MD5_LIBNAME) $(MD5_OBJS)

src/$(DES56_LIBNAME) : $(DES56_OBJS)
	export MACOSX_DEPLOYMENT_TARGET="10.3"; $(CC) $(CFLAGS) $(LIB_OPTION) -o src/$(DES56_LIBNAME) $(DES56_OBJS)

$(COMPAT_DIR)/compat-5.1.o: $(COMPAT_DIR)/compat-5.1.c
	$(CC) -c $(CFLAGS) -o $@ $(COMPAT_DIR)/compat-5.1.c

install: src/$(MD5_LIBNAME) src/$(DES56_LIBNAME)
	mkdir -p $(LUA_LIBDIR)/md5
	cp src/$(MD5_LIBNAME) $(LUA_LIBDIR)/md5/core.so
	mkdir -p $(LUA_DIR)
	cp $(MD5_LUAS) $(LUA_DIR)
	cp src/$(DES56_LIBNAME) $(LUA_LIBDIR)

clean:
	rm -f $(MD5_OBJS) src/$(MD5_LIBNAME) $(DES56_OBJS) src/$(DES56_LIBNAME)
