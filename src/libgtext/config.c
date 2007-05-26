/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
#include <assert.h>
#include <libgtcore/str.h>
#include <libgtcore/ensure.h>
#include <libgtext/config.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

struct Config {
	lua_State *L;
  Str *fn;
};

#define COLORS_ARE_EQUAL(c1, c2) \
  ((c1.red == c2.red) && (c1.green == c2.green) && (c1.blue == c2.blue))

#define NEWTABLE(key) \
  lua_pushstring(cfg->L, key); \
	lua_newtable(cfg->L); \
	lua_settable(cfg->L, -3); \

#define ADDCOLOR(rgb) \
	lua_pushstring(cfg->L, #rgb); \
	lua_pushnumber(cfg->L, color.rgb); \
	lua_settable(cfg->L, -3);

#define GET_AND_SET_COLOR_VALUE(rgb) \
  lua_getfield(cfg->L, -1, #rgb); \
  if(lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) { \
    env_log_log(env, "%s value for type '%s' is undefined or not numeric, " \
		                 "using default",#rgb, key); \
  } \
  else \
  { \
    color.rgb = lua_tonumber(cfg->L,-1); \
  } \
  lua_pop(cfg->L, 1); \

Config* config_new(Env *env)
{
	Config *cfg;
  env_error_check(env);
  cfg = env_ma_malloc(env, sizeof (Config));
	cfg->fn = NULL;
  cfg->L = luaL_newstate();
  if (!cfg->L) 
  {
    env_error_set(env, "out of memory (cannot create new lua state)");
  }
	else luaL_openlibs(cfg->L);
  return cfg;
}

/*!
Deletes a Config object.
\param cfg Pointer to Config object to delete.
\param env Pointer to Environment object.
*/
void config_delete(Config *cfg, Env *env)
{
  assert(cfg);
	if(cfg->L) lua_close(cfg->L);
	if(cfg->fn) str_delete(cfg->fn,env);
  env_ma_free(cfg, env);
}

void config_load_file(Config *cfg, Str *fn, Env* env)
{
  int has_err = 0;
  env_error_check(env);
  assert(cfg && cfg->L && fn);
  cfg->fn = str_ref(fn);
  if (luaL_loadfile(cfg->L, str_get(fn)) ||
      lua_pcall(cfg->L, 0, 0, 0)) 
	{
    env_error_set(env, "cannot run configuration file: %s",
                  lua_tostring(cfg->L, -1));
    has_err = -1;
  }
  if (!has_err) 
	{
    lua_getglobal(cfg->L, "config");
    if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) 
		{
      env_error_set(env, "'config' is not defined or not a table in \"%s\"",
                str_get(fn));
      has_err = -1;
    }
		lua_pop(cfg->L, 1);
  }
}

void config_set_color(Config *cfg, const char *key, Color color, Env* env)
{
  assert(cfg && key);
	env_error_check(env);
	lua_getglobal(cfg->L, "config");
	if(lua_isnil(cfg->L, -1))
	{
	  lua_newtable(cfg->L);
		lua_setglobal(cfg->L, "config");
		lua_getglobal(cfg->L, "config");
	}
	lua_getfield(cfg->L, -1, "colors");
	if(lua_isnil(cfg->L, -1))
	{
	  lua_pop(cfg->L, 1);
    NEWTABLE("colors")
		lua_getfield(cfg->L, -1, "colors");
	}
	lua_getfield(cfg->L, -1, key);
	if(lua_isnil(cfg->L, -1))
	{
	  lua_pop(cfg->L, 1);
    NEWTABLE(key);
		lua_getfield(cfg->L, -1, key);
	}
  ADDCOLOR(red);
	ADDCOLOR(green);
	ADDCOLOR(blue);
	lua_pop(cfg->L, 3);
}

Color config_get_color(Config *cfg, const char *key, Env* env)
{
  assert(cfg && key);
  Color color;
	int i=0;
	env_error_check(env);
	/* set default colors */
	color.red=0.8; color.green = 0.8; color.blue=0.8;
	/* get config table */
	lua_getglobal(cfg->L, "config");
	if(lua_isnil(cfg->L, -1))
	{
	  env_log_log(env, "'config' is not defined or not a table,"
		                 " will use defaults");
	  return color;
	}
	i++;
	/* get colors section */
	lua_getfield(cfg->L, -1, "colors");
	if(lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) 
	{
    env_log_log(env, "section 'colors' is not defined (or not a table) "
		                 " in configuration file, will use defaults");
		lua_pop(cfg->L, 1);
		return color;
	} else i++;
	/* lookup color entry for given feature */
  lua_getfield(cfg->L, -1, key);
	if(lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1)) 
	{
    env_log_log(env, "no colors are defined for type '%s', will use defaults",
		            key);
		lua_pop(cfg->L, 1);
		return color;
  } else i++;
	/* update color struct */
  GET_AND_SET_COLOR_VALUE(red);
	GET_AND_SET_COLOR_VALUE(green);
	GET_AND_SET_COLOR_VALUE(blue);
	/* reset stack to original state for subsequent calls */
	lua_pop(cfg->L, i);
	return color;
}

void config_reload(Config *cfg, Env *env)
{
  config_load_file(cfg, cfg->fn, env);
}

int config_unit_test(Env* env)
{
  int has_err = 0;
	Config *cfg;
	Str *luafile = str_new_cstr("test.lua",env);
	Color col1, col2, col, defcol, tmpcol; 

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;
  col2.red=.4;col2.green=.5;col2.blue=.6;
  col.red=1;col.green=1;col.blue=1;
	defcol.red=.8;defcol.green=.8;defcol.blue=.8;

  /* instantiate new config object */
	cfg = config_new(env);
	
	/* at the beginning, all colors are defaults, since nothing is defined */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "cds", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,defcol));

  /* execute the lua test file */
	config_load_file(cfg, luafile, env);
 
  /* now we expect the colors to exist and have certain values */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, !COLORS_ARE_EQUAL(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,col1));
  tmpcol = config_get_color(cfg, "cds", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,col2));

  /* change a color... */
  config_set_color(cfg, "exon", col, env);
	
	/* is it saved correctly? */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, !COLORS_ARE_EQUAL(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,col));

  /* create a new color definition */
  config_set_color(cfg, "foo", col, env);
	
	/* is it saved correctly? */
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, !COLORS_ARE_EQUAL(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, COLORS_ARE_EQUAL(tmpcol,col));

  /* mem cleanup */
  str_delete(luafile, env);
	config_delete(cfg, env);
	
	return has_err;
}
