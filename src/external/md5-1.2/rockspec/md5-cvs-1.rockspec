package = "MD5"
version = "cvs-1"
source = {
   cvs_tag = "HEAD",
   url = "cvs://:pserver:anonymous@cvs.luaforge.net:/cvsroot/md5"
}
description = {
   summary = "Basic cryptographic library",
   detailed = [[
      MD5 offers basic cryptographic facilities for Lua 5.1: a hash
      (digest) function, and a pair crypt/decrypt.
   ]],
   license = "MIT/X11",
   homepage = "http://www.keplerproject.org/md5/"
}
dependencies = {
   "lua >= 5.1"
}
build = {
   type = "make",
   variables = {
      LUA_VERSION_NUM="501",
   },
   build_variables = {
     LIB_OPTION = "$(LIBFLAG)",
     CFLAGS = "$(CFLAGS) -I$(LUA_INCDIR)",
   },
   install_variables = {
      LUA_LIBDIR = "$(LIBDIR)",
      LUA_DIR = "$(LUADIR)"
   },
   platforms = {
     win32 = {
       build_variables = {
	 LUA_LIB = "$(LUA_LIBDIR)\\lua5.1.lib"
       }
     }
   }
}
