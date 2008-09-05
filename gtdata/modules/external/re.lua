-- $Id: re.lua,v 1.26 2008/03/07 14:20:55 roberto Exp $

local m = require"lpeg"
local _G = _G
local tonumber, type, print, error = tonumber, type, print, error
local mt = getmetatable(m.P(0))

module "re"

local any = m.P(1)

-- Pre-defined names
Predef = {
  l = m.R"az",
  u = m.R"AZ",
  d = m.R"09",
  x = m.R("09", "AF", "af"),
  s = m.S"\n\r\t\v\f ",
  p = m.S"!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~",
  nl = m.P"\n", 
}

Predef.a = Predef.l + Predef.u
Predef.w = Predef.a + Predef.d
Predef.g = Predef.w + Predef.p
Predef.c = any - (Predef.g + ' ')

Predef.L = any - Predef.l
Predef.U = any - Predef.u
Predef.D = any - Predef.d
Predef.X = any - Predef.x
Predef.S = any - Predef.s
Predef.P = any - Predef.p
Predef.A = any - Predef.a
Predef.W = any - Predef.w
Predef.G = any - Predef.g
Predef.C = any - Predef.c


local I = m.P(function (s,i) print(i, s:sub(1, i-1)); return i end)


local function getdef (id, Defs)
  local c = Defs and Defs[id]
  if not c then error("undefined name: " .. id) end
  return c
end


local function patt_error (s, i)
  local msg = (#s < i + 20) and s:sub(i)
                             or s:sub(i,i+20) .. "..."
  msg = ("pattern error near '%s'"):format(msg)
  error(msg, 2)
end

local function mult (p, n)
  local np = m.P(true)
  while n >= 1 do
    if n%2 >= 1 then np = np * p end
    p = p * p
    n = n/2
  end
  return np
end

local function equalcap (s, i, c)
  if type(c) ~= "string" then return nil end
  local e = #c + i
  if s:sub(i, e - 1) == c then return e else return nil end
end


local S = (m.S(" \t\n") + "--" * (any - m.S"\n")^0)^0

-- identifiers allways need an environment
local Identifier = m.R("AZ", "az") * m.R("AZ", "az", "09")^0

local exp_follow = m.P"/" + ")" + "}" + "~}" + -1 + "%" + "<" + Identifier

Identifier = m.C(Identifier) * m.Carg(1)

local num = m.C(m.R"09"^1) * S / tonumber

local String = "'" * m.C((any - "'")^0) * "'" +
               '"' * m.C((any - '"')^0) * '"'

local Range = m.Cs(any * (m.P"-"/"") * (any - "]")) / m.R

local Cat = "%" * Identifier / function (c,Defs)
  local cat =  Defs and Defs[c] or Predef[c]
  if not cat then error ("name '" .. c .. "' undefined") end
  return cat
end


local item = Cat + Range + m.C(any)

local Class =
    "["
  * (m.C(m.P"^"^-1))    -- optional complement symbol
  * m.Ca(item * ((item - "]") / mt.__add)^0) /
                          function (c, p) return c == "^" and any - p or p end
  * "]"

local function adddef (t, k, Defs, exp)
  if t[k] then
    error("'"..k.."' already defined as a rule")
  else
    t[k] = exp
  end
  return t
end

local function firstdef (n, Defs, r) return adddef({n}, n, Defs, r) end



local exp = m.P{ "Exp",
  Exp = S * ( m.V"Grammar"
            + m.Ca(m.V"Seq" * ("/" * S * m.V"Seq" / mt.__add)^0) );
  Seq = m.Ca(m.Cc(m.P"") * (m.V"Prefix" / mt.__mul)^0)
        * (#exp_follow + patt_error);
  Prefix = "&" * S * m.V"Prefix" / mt.__len
         + "!" * S * m.V"Prefix" / mt.__unm
         + m.V"Suffix";
  Suffix = m.Ca(m.V"Primary" * S *
          ( ( m.P"+" * m.Cc(1) / mt.__pow        -- patt^1
            + m.P"*" * m.Cc(0) / mt.__pow        -- patt^0
            + m.P"?" * m.Cc(-1) / mt.__pow   -- patt^-1
            + "^" * ( num / mult
                    + "+" * num / mt.__pow
                    + "-" * num / function (patt,n) return patt^-n end )
            + "->" * S * ( String / mt.__div
                         + m.P"{}" / m.Ct
                         + (Identifier / getdef) / mt.__div )
            + "=>" * S * (Identifier / getdef) / m.Cmt
            ) * S
          )^0 );
  Primary = "(" * m.V"Exp" * ")"
            + m.P"{}" / m.Cp
            + "{~" * m.V"Exp" * "~}" / m.Cs
            + "{" * m.V"Exp" * "}" / m.C
            + String / m.P
            + Class
            + Cat
            + "%" * num / function (n) return m.Cmt(m.Cb(n), equalcap) end
            + m.P"." * m.Cc(any)
            + "<" * Identifier * ">" / m.V;
  Definition = Identifier * S * '<-' * m.V"Exp";
  Grammar = m.Ca(m.V"Definition" / firstdef * (m.V"Definition" / adddef)^0) /
                m.P
}

local pattern = S * exp / m.P * (-any + patt_error)
                                   


function compile (p, defs)
  if m.type(p) == "pattern" then return p end   -- already compiled
  local cp = pattern:match(p, 1, defs)
  if not cp then error("incorrect pattern", 3) end
  return cp
end


local mem = {}
local fmem = {}
local gmem = {}
local mt = {__mode = "v"}
_G.setmetatable(mem, mt)
_G.setmetatable(fmem, mt)
_G.setmetatable(gmem, mt)

function match (s, p, i)
  local cp = mem[p]
  if not cp then
    cp = compile(p)
    mem[p] = cp
  end
  return cp:match(s, i or 1)
end

function find (s, p, i)
  local cp = fmem[p]
  if not cp then
    cp = compile(p)
    cp = m.P{ m.Cp() * cp + 1 * m.V(1) }
    fmem[p] = cp
  end
  return cp:match(s, i or 1)
end

function gsub (s, p, rep)
  gmem[p] = gmem[p] or {}
  local cp = gmem[p][rep]
  if not cp then
    cp = compile(p)
    cp = m.Cs((cp / rep + 1)^0)
    gmem[p][rep] = cp
  end
  return cp:match(s)
end
