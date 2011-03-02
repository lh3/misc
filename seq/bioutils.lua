#!/usr/bin/env luajit

-- klua routines
function string:split(sep, n)
	local a, start = {}, 1;
	sep = sep or "%s+";
	repeat
		local b, e = self:find(sep, start);
		if b == nil then
			table.insert(a, self:sub(start));
			break
		end
		a[#a+1] = self:sub(start, b - 1);
		start = e + 1;
		if n and #a == n then
			table.insert(a, self:sub(start));
			break
		end
	until start > #self;
	return a;
end

function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end
--

--[[
	BED utilities
]]--

function bed_merge(arg)
	local c, b, e, n;
	local fp = io.xopen(arg[1]);
	for l in fp:lines() do
		local t = l:split('\t', 4);
		t[2], t[3] = tonumber(t[2]), tonumber(t[3]);
		if t[2] < 0 then t[2] = 0 end
		if t[1] ~= c or t[2] > e then -- no overlap
			if c ~= nil then print(c, b, e) end
			c, b, e = t[1], t[2], t[3]
		else
			e = (e < t[3] and  t[3]) or e;
		end
	end
	print(c, b, e);
	fp:close();
end

function bed_cov2(arg)
	local c, e, ob, oe; -- e is the rightmost pos;
	local fp = io.xopen(arg[1]);
	ob, oe = -1, -1; -- ob and oe are the coordinates of the region to output
	for l in fp:lines() do
		local t = l:split('\t', 4);
		t[2], t[3] = tonumber(t[2]), tonumber(t[3]);
		if t[2] < 0 then t[2] = 0 end
		if t[1] ~= c or t[2] >= e then -- no overlap; output
			if ob ~= -1 then print(c, ob, oe) end
			c, e, ob, oe = t[1], t[3], -1, -1;
		else -- overlap
			if ob == -1 or t[2] >= oe then -- not that if ob==-1, oe<=t[2]<e here; to output
				if ob ~= -1 then print(c, ob, oe) end
				ob, oe = t[2], (t[3] < e and t[3]) or e;
			else -- update oe; this only happens if a region is covered by 3 or more times
				local new_oe = (t[3] < e and t[3]) or e;
				if oe < new_oe then oe = new_oe end
			end
			if e < t[3] then e = t[3] end
		end
	end
	if ob ~= -1 then print(c, ob, oe) end
	fp:close();
end

--[[
	PLIST utilities
]]--

function pos_read_list(fn)
	local fp = io.xopen(fn);
	local hash = {};
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		hash[chr .. pos] = 1;
	end
	fp:close();
	return hash;
end

function pos_uniq(arg)
	local fp, hash = io.xopen(arg[1]), pos_read_list(arg[2]);
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		if hash[chr .. pos] == nil then print(l) end
	end
	fp:close();
end

function pos_joint(arg)
	local fp, hash = io.xopen(arg[1]), pos_read_list(arg[2]);
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		if hash[chr .. pos] then print(l) end
	end
	fp:close();
end

--[[
]]--

function misc_N50(arg)
	local col, list = 1, {};
	local fp = io.xopen(arg[1]);
	local x = 0;
	if arg[2] ~= nil then col = tonumber(arg[2]) end
	for l in fp:lines() do
		local t = l:split("%s+", col+1);
		local y = tonumber(t[col]);
		x = x + y;
		table.insert(list, y);
	end
	local sum = x;
	table.sort(list);
	x = 0;
	for i = #list,1,-1 do
		x = x + list[i];
		if x >= sum/2 then
			print(#list, sum, list[i]);
			return
		end
	end
end

-- main function

if #arg == 0 then
	print('');
	print('Usage:   bioutils.lua <command> <arguments>\n');
	print('Command: bedmerge  merge overlapping regions in a sorted BED');
	print('         bedcov2   extract regions covered by >=2 records in a sorted BED');
	print('');
	print('         posuniq   output records unique in the 1st PLIST file');
	print('         posjoint  output records in the 1st PLIST present in the 2nd');
	print('');
	print('         N50       calculate N50');
	print('');
	os.exit(1);
end

local cmd = arg[1];
table.remove(arg, 1);
if cmd == 'bedmerge' then bed_merge(arg)
elseif cmd == 'bedcov2' then bed_cov2(arg)
elseif cmd == 'posuniq' then pos_uniq(arg)
elseif cmd == 'posjoint' then pos_joint(arg)
elseif cmd == 'N50' then misc_N50(arg);
else print('Unrecognized command.') end
