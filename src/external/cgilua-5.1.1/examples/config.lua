-- CGILua example configuration file
-- $Id: config.lua,v 1.2 2007/10/31 20:50:53 carregal Exp $
--
-- You may want to use parts of the this file in your CGILua configuration
--
--

-- Emulating previous CGILua behavior loading file "env.lua" from the script's directory
--[[
cgilua.addopenfunction (function ()
	cgilua.doif ("env.lua")
end)
--]]

-- Optional backward compatibility with the global cgi table
--cgi = {}
--for i,v in pairs(cgilua.POST) do
--       cgilua.urlcode.insertfield (cgi, i, v)
--end
--for i,v in pairs(cgilua.QUERY) do
--       cgilua.urlcode.insertfield (cgi, i, v)
--end

-- Basic configuration for using sessions
-- require"cgilua.session"
-- cgilua.session.setsessiondir (CGILUA_TMP)
-- Add cgilua.enablesession() at the beginning of every script which depends
-- on sessions.

-- Compatibility

cgilua.preprocess = cgilua.handlelp
cgilua.includehtml = cgilua.lp.include

-- Directories for applications' libraries.
-- The following table should be indexed by the virtual path of the application
-- and contain the absolute path of the application's Lua-library directory.
--[[
local app_lib_dir = {
	["/t/"] = "/usr/local/src/cgilua/tests",
}
local package = package
cgilua.addopenfunction (function ()
	local app = app_lib_dir[cgilua.script_vdir]
	if app then
		package.path = app.."/?.lua;"..package.path
	end
end)
--]]