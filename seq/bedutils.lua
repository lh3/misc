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

function bed_merge(fn)
	local c, b, e, n;
	local fp = io.xopen(fn);
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

function bed_cov2(fn)
	local c, e, ob, oe; -- e is the rightmost pos;
	local fp = io.xopen(fn);
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

function plist_read_list(fn)
	local fp = io.xopen(fn2);
	local hash = {};
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		hash[chr .. pos] = 1;
	end
	fp:close();
	return hash;
end

function plist_uniq(fn1, fn2)
	local fp, hash = io.xopen(fn1), plist_read_list(fn2);
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		if hash[chr .. pos] == nil then print(l) end
	end
	fp:close();
end

function plist_joint(fn1, fn2)
	local fp, hash = io.xopen(fn1), plist_read_list(fn2);
	for l in fp:lines() do
		local chr, pos = l:match("^(%S+)%s+(%d+)");
		if hash[chr .. pos] then print(l) end
	end
	fp:close();
end

-- main function

if #arg == 0 then
	print('');
	print('Usage:   bedutils.lua <command> <arguments>\n');
	print('Command: merge     merge overlapping regions in a sorted BED');
	print('         cov2      extract regions covered by >=2 records in a sorted BED');
	print('');
	print('         luniq     output records unique in the 1st PLIST file');
	print('         ljoint    output records in the 1st PLIST present in the 2nd');
	print('');
	os.exit(1);
end

local cmd = arg[1];
table.remove(arg, 1);
if cmd == 'merge' then bed_merge(arg[1])
elseif cmd == 'cov2' then bed_cov2(arg[1])
elseif cmd == 'luniq' then plist_uniq(arg[1], arg[2])
elseif cmd == 'ljoint' then plist_joint(arg[1], arg[2])
else print('Unrecognized command.') end
