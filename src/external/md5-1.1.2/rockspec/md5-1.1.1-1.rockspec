package = "MD5"
version = "1.1.1-1"
source = {
   url = "http://luaforge.net/frs/download.php/2746/md5-1.1.1.tar.gz"
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
   }
}
