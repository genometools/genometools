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

require 'gtdoclib'

function usage()
  io.stderr:write(string.format("Usage: [-html | -v] %s gt_home\n", arg[0]))
  io.stderr:write("Generate documentation for the GenomeTools home directory " ..
                  "gt_home.\n")
  os.exit(1)
end

local be_verbose    = false
local out_mode      = "txt"

if #arg >= 1 then
  while #arg >= 1 and string.match(arg[1], "^-") do
    if string.match(arg[1], "^-h") then
      out_mode = "html"
      table.remove(arg, 1)
    elseif string.match(arg[1], "^-v") then
      be_verbose = true
      table.remove(arg, 1)
    end
  end
  if #arg == 1 then
    gt_home = arg[1]
  else
    usage()
  end
else
  usage()
end

local template_path = gt_home .. "/gtdata/modules/gtdoclib/"

local export = { "src/libgtlua",
                 "gtdata/modules/gtlua.lua",
                 "gtdata/modules/gtlua" }

local doc_parser      = DocParser:new()
local doc_base        = DocBase:new()

local function show_rec_array(array)
  assert(array)
  for _, v in ipairs(array) do
    if type(v) == "table" then
      show_rec_array(v)
    else
      print(v)
    end
  end
end

local function process_file(filename, be_verbose)
  assert(filename)
  if is_header(filename) or is_lua_file(filename) then
    local ast = doc_parser:parse(filename, be_verbose)
    if be_verbose then
      print("showing ast:")
      show_rec_array(ast)
    end
    doc_base:process_ast(ast, be_verbose)
  end
end

for _, v in ipairs(export) do
  local filename = gt_home .. "/" .. v
  if is_dir(filename) then
    for f in lfs.dir(filename) do
      local filename = filename .. "/" .. f
      process_file(filename, be_verbose)
    end
  else
    process_file(filename, be_verbose)
  end
end

local doc_visitor
if out_mode == "txt" then
  doc_visitor = DocVisitorTxt:new()
  doc_base:accept(doc_visitor)
else
  assert(out_mode == "html")
  doc_visitor = DocVisitorHTML:new(template_path)
  doc_visitor:show_header()
  doc_base:accept(doc_visitor)
  doc_visitor:show_footer()
end
