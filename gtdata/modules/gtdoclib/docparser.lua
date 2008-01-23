--[[
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
]]

module(..., package.seeall)

require 'lpeg'
require 'fileutils'

DocParser = {}

-- Common Lexical Elements
local Any           = lpeg.P(1)
local Newline       = lpeg.S("\n")
local Whitespace    = lpeg.S(" \t\n")
local OptionalSpace = Whitespace^0
local Space         = Whitespace^1

-- Lexical Elements of Lua
local LuaLongCommentStart  = lpeg.P("--[[")
local LuaLongCommentEnd    = lpeg.P("]]")
local LuaLongComment       = lpeg.Cc("long comment") *
                             lpeg.C(LuaLongCommentStart *
                                    (Any - LuaLongCommentEnd)^0 *
                                    LuaLongCommentEnd)
local LuaShortCommentStart = lpeg.P("--")
local LuaShortCommentEnd   = Newline
local LuaShortCommentLine  = OptionalSpace * LuaShortCommentStart *
                             lpeg.C((Any - LuaShortCommentEnd)^0) *
                             LuaShortCommentEnd
local LuaShortComment      = lpeg.Cc("comment") * LuaShortCommentLine^1
local LuaOptionalComment   = LuaShortComment +
                             lpeg.Cc("comment") * lpeg.Cc("undefined")
local LuaCommentStart      = LuaLongCommentStart + LuaShortCommentStart
local LuaEnd               = lpeg.P("end")
local LuaLocalFunction    =  lpeg.P("local") * OptionalSpace *
                             lpeg.P("function") * OptionalSpace *
                             lpeg.P(Any - lpeg.P("("))^1 *
                             lpeg.P("(") * (Any - lpeg.P(")"))^0 *
                             lpeg.P(")")
local LuaGlobalFunction    = -lpeg.P("local") * OptionalSpace *
                             lpeg.C(lpeg.P("function")) * OptionalSpace *
                             lpeg.C(lpeg.P(Any - lpeg.P("("))^1) *
                             lpeg.P("(") * lpeg.C((Any - lpeg.P(")"))^0) *
                             lpeg.P(")")
local ExportLuaMethod      = lpeg.Ct(LuaOptionalComment * LuaGlobalFunction)
local CodeStop             = LuaCommentStart + LuaLocalFunction + LuaGlobalFunction
local LuaCode              = lpeg.Cc("code") * lpeg.C((Any - CodeStop)^1)

-- Lexical Elements of C
local Character      = lpeg.R("AZ", "az")
local CCommentStart  = lpeg.P("/*")
local CCommentEnd    = lpeg.P("*/")
local ExportCComment = CCommentStart * lpeg.P(" exports the ") *
                       lpeg.Ct(lpeg.Cc("class") * lpeg.C(Character^1)) *
                       (Any - lpeg.P("to Lua:"))^1 *
                       lpeg.P("to Lua:") * ExportLuaMethod^0 *
                       (Any - CCommentEnd)^0 * CCommentEnd
local CComment       = CCommentStart * (Any - CCommentEnd)^0 * CCommentEnd
local CCode          = (Any - CCommentStart)^1

-- Lua Grammar
local Elem, Start = lpeg.V"Elem", lpeg.V"Start"
local LuaGrammar = lpeg.P{ Start,
  Start = lpeg.Ct(Elem^0);
  Elem  = ExportLuaMethod + LuaLongComment + LuaShortComment + Space +
          LuaLocalFunction + LuaCode;
}
LuaGrammar = LuaGrammar * -1

-- C Grammar
local CGrammar = lpeg.P{ Start,
 Start = lpeg.Ct(Elem^0);
 Elem  = lpeg.Ct(ExportCComment) + CComment + Space + CCode;
}
CGrammar = CGrammar * -1

function DocParser:new()
  o = {}
  o.c_pattern   = CGrammar
  o.lua_pattern = LuaGrammar
  setmetatable(o, self)
  self.__index = self
  return o
end

function DocParser:parse(filename, be_verbose)
  assert(filename)
  assert(is_header(filename) or is_lua_file(filename))
  if be_verbose then
    print("parsing " .. filename)
  end
  local file, err = io.open(filename, "r")
  assert(file, err)
  local filecontent = file:read("*a")
  if is_header(filename) then
    return lpeg.match(self.c_pattern, filecontent)
  else -- is_lua_file(filename)
    return lpeg.match(self.lua_pattern, filecontent)
  end
end
