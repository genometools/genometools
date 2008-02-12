/*
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
*/

/*
  This interactive Lua interpreter was heavily inspired by the original Lua
  interpreter (lua.c), therefore the copyright notice in lua.h also applies.
*/

#include <assert.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "lauxlib.h"
#include "libtecla.h"
#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtlua/interactive.h"

static lua_State *globalL = NULL;

static void lstop(lua_State *L, UNUSED lua_Debug *ar) {
  lua_sethook(L, NULL, 0, 0);
  luaL_error(L, "interrupted!");
}

static void laction (int i) {
  signal(i, SIG_DFL); /* if another SIGINT happens before lstop,
                         terminate process (default action) */
  lua_sethook(globalL, lstop, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}

static int report(lua_State *L, int status) {
  if (status && !lua_isnil(L, -1)) {
    const char *msg = lua_tostring(L, -1);
    if (msg == NULL) msg = "(error object is not a string)";
    fprintf(stderr, "%s\n", msg);
    fflush(stderr);
    lua_pop(L, 1);
  }
  return status;
}

static int traceback(lua_State *L) {
  lua_getfield(L, LUA_GLOBALSINDEX, "debug");
  if (!lua_istable(L, -1)) {
    lua_pop(L, 1);
    return 1;
  }
  lua_getfield(L, -1, "traceback");
  if (!lua_isfunction(L, -1)) {
    lua_pop(L, 2);
    return 1;
  }
  lua_pushvalue(L, 1);  /* pass error message */
  lua_pushinteger(L, 2);  /* skip this function and traceback */
  lua_call(L, 2, 1);  /* call debug.traceback */
  return 1;
}

static int docall(lua_State *L, int narg, int clear) {
  int status;
  int base = lua_gettop(L) - narg;  /* function index */
  lua_pushcfunction(L, traceback);  /* push traceback function */
  lua_insert(L, base);  /* put it under chunk and args */
  signal(SIGINT, laction);
  status = lua_pcall(L, narg, (clear ? 0 : LUA_MULTRET), base);
  signal(SIGINT, SIG_DFL);
  lua_remove(L, base);  /* remove traceback function */
  /* force a complete garbage collection in case of errors */
  if (status != 0) lua_gc(L, LUA_GCCOLLECT, 0);
  return status;
}

static const char* get_prompt(lua_State *L, bool firstline) {
  const char *p;
  lua_getfield(L, LUA_GLOBALSINDEX, firstline ? "_PROMPT" : "_PROMPT2");
  p = lua_tostring(L, -1);
  if (p == NULL) p = (firstline ? "> " : ">> ");
  lua_pop(L, 1);  /* remove global */
  return p;
}

static int incomplete(lua_State *L, int status) {
  if (status == LUA_ERRSYNTAX) {
    size_t lmsg;
    const char *msg = lua_tolstring(L, -1, &lmsg);
    const char *tp = msg + lmsg - (sizeof (LUA_QL("<eof>")) - 1);
    if (strstr(msg, LUA_QL("<eof>")) == tp) {
      lua_pop(L, 1);
      return 1;
    }
  }
  return 0;  /* else... */
}

static int pushline(lua_State *L, bool firstline, UNUSED GetLine *gl) {
  char buffer[BUFSIZ];
  char *b = buffer;
  size_t l;
  const char *prmt = get_prompt(L, firstline);
#ifdef CURSES
  if (!(b = gl_get_line(gl, prmt, NULL, 0)))
    return 0; /* no input */
#else
  b = buffer; /* use static buffer */
  fputs(prmt, stdout); fflush(stdout);/* show prompt */
  if (!fgets(b, BUFSIZ, stdin))  /* get line */
    return 0; /* no input */
#endif
  l = strlen(b);
  if (l > 0 && b[l-1] == '\n')  /* line ends with newline? */
    b[l-1] = '\0';  /* remove it */
  if (firstline && b[0] == '=')  /* first line starts with `=' ? */
    lua_pushfstring(L, "return %s", b+1);  /* change it to `return' */
  else
    lua_pushstring(L, b);
  return 1;
}

static int loadline(lua_State *L, GetLine *gl) {
  int status;
  lua_settop(L, 0);
  if (!pushline(L, true, gl))
    return -1;  /* no input */
  for (;;) {  /* repeat until gets a complete line */
    status = luaL_loadbuffer(L, lua_tostring(L, 1), lua_strlen(L, 1), "=stdin");
    if (!incomplete(L, status)) break;  /* cannot try to add lines? */
    if (!pushline(L, false, gl))  /* no more input? */
      return -1;
    lua_pushliteral(L, "\n");  /* add a new line... */
    lua_insert(L, -2);  /* ...between the two lines */
    lua_concat(L, 3);  /* join them */
  }
#ifdef CURSES
  if (lua_strlen(L, 1) > 0) { /* non-empty line? */
    char *line;
    int rval;
    /* save complete line in history */
    line = cstr_dup(lua_tostring(L, 1));
    cstr_rep(line, '\n', ' '); /* replace all newlines in <line> with blanks,
                                  because otherwise gl_append_history() would
                                  truncate <line> at the first newline */
    rval = gl_append_history(gl, line);
    assert(!rval);
    ma_free(line);
  }
  lua_remove(L, 1);  /* remove line */
#endif
  return status;
}

static void dotty(lua_State *L, GetLine *gl) {
  int status;
  while ((status = loadline(L, gl)) != -1) {
    if (status == 0) status = docall(L, 0, 0);
    report(L, status);
    if (status == 0 && lua_gettop(L) > 0) {  /* any result to print? */
      lua_getglobal(L, "print");
      lua_insert(L, 1);
      if (lua_pcall(L, lua_gettop(L)-1, 0, 0) != 0)
        fprintf(stderr, "%s\n", lua_pushfstring(L,
                               "error calling " LUA_QL("print") " (%s)",
                               lua_tostring(L, -1)));
        fflush(stderr);
    }
  }
  lua_settop(L, 0);  /* clear stack */
  fputs("\n", stdout);
  fflush(stdout);
}

void run_interactive_lua_interpreter(lua_State *L)
{
#ifdef CURSES
  GetLine *gl;
  gl = new_GetLine(2048, 8096);
  gl_automatic_history(gl, 0); /* disable automatic history saving, we save
                                  complete input lines explicitly */
  if (!gl) {
    fprintf(stderr, "cannot create GetLine object\n");
    exit(EXIT_FAILURE);
  }
#endif
  globalL = L; /* for signal handling */
#ifdef CURSES
  dotty(L, gl);
  del_GetLine(gl);
#else
  dotty(L, NULL);
#endif
}
