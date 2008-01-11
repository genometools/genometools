-- $Id: re.lua,v 1.18 2007/10/10 18:53:45 roberto Exp $

local m = require"lpeg"
local _G = _G
local print, error = print, error
local mt = getmetatable(m.P(0))

module "re"

local I = m.P(function (s,i) print(i, s:sub(1, i-1)); return i end)

local any = m.P(1)

local function complement (c, p) return c and any - p or p end


local Defs


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

local function minmax (p, min, max)
  min = _G.tonumber(min)
  if max == true then return p^min
  else
    max = (max == false) and min or _G.tonumber(max)
    if max < min then return m.P(false) end
    local r = mult(p, min)
    if max > min then r = r * p^(min - max) end
    return r
  end
end


local S = (m.S(" \t\n") + "--" * (any - m.S"\n")^0)^0

local String = "'" * m.C((any - "'")^0) * "'" +
               '"' * m.C((any - '"')^0) * '"'

local Range = m.Cs(any * (m.P"-"/"") * (any - "]")) / m.R

local item = Range + m.C(any)

local Class =
    "["
  * ("^" * m.Cc(true) + m.Cc(false))   -- optional complement symbol
  * m.Ca(item * ((item - "]") / mt.__add)^0) / complement
  * "]"

local Identifier = m.R("AZ", "az", "__") * m.R("AZ", "az", "__", "09")^0
local num = m.C(m.R"09"^1) * S

-- {num} or {num,} or {num,num}
local rep = m.P"{" * S * num * ("," * S * (num + m.Cc(true)) + m.Cc(false))
          * "}"

local exp_follow = m.P"/" + ")" + "}" + "~}" + -1 + Identifier

local exp = m.P{ "Exp",
  Exp = S * m.Ca(m.V"Seq" * ("/" * S * m.V"Seq" / mt.__add)^0);
  Seq = m.Ca(m.Cc(m.P"") * (m.V"Prefix" / mt.__mul)^0)
        * (#exp_follow + patt_error);
  Prefix = "&" * S * m.V"Prefix" / mt.__len
         + "!" * S * m.V"Prefix" / mt.__unm
         + m.V"Sufix";
  Sufix = m.Ca(m.V"Primary" * S *
          ( ( m.P"+" * m.Cc(1, true) / minmax    -- 1 -> true
            + m.P"*" * m.Cc(0, true) / minmax  -- 0 -> true
            + m.P"?" * m.Cc(0, 1) / minmax     -- 0 -> 1
            + rep / minmax
            + "->" * S * ( String / mt.__div
                         + m.P"{}" / m.Ct
                         + m.C(Identifier) / function (p,id)
                             local c = Defs and Defs[id]
                             if not c then error("undefined name: "..id) end
                             return p / c
                           end )
            ) * S
          )^0 );
  Primary = "(" * m.V"Exp" * ")"
            + m.P"{}" / m.Cp
            + "{~" * m.V"Exp" * "~}" / m.Cs
            + "{" * m.V"Exp" * "}" / m.C
            + String / m.P
            + Class
            + m.P"." * m.Cc(any)
            + Identifier * -(S * '<-') / function (n)
                return Defs and m.P(Defs[n]) or m.V(n)
              end;
}


local function adddef (t, k, exp)
  if Defs and Defs[k] then
    error("'"..k.."' defined both externally and as a rule")
  elseif t[k] then
    error("'"..k.."' already defined as a rule")
  else
    t[k] = exp
  end
  return t
end


local definition = m.C(Identifier) * S * '<-' * exp

local grammar = m.Ca(
     definition / function (n, r) return adddef({n}, n, r) end
   * (definition / adddef)^0
)

local pattern = S * (grammar + exp) / m.P * (-any + patt_error)
                                   


function compile (p, defs)
  Defs = defs
  local cp = pattern:match(p)
  if not cp then error("incorrect pattern", 3) end
  return cp
end


local mem = {}
_G.setmetatable(mem, {__mode = "v"})

function match (s, p)
  local cp = mem[p]
  if not cp then
    cp = compile(p)
    mem[p] = cp
  end
  return cp:match(s)
end

