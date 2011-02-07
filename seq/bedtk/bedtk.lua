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
	elseif fn:sub(-3) == '.gz' then return io.popen('gzip -dc ' .. fn, mode);
	elseif fn:sub(-4) == '.bz2' then return io.popen('bzip2 -dc ' .. fn, mode);
	else return io.open(fn, mode) end
end
--

function bed_read2(fp1, fp2)
	if fp2 == nil then
		return function()
			local t = fp1:read():split("\t", 4);
			t[2], t[3] = tonumber(t[2]), tonumber(t[3]);
			return t;
		end
	else
		local fp, t, i = {fp1, fp2}, {}, 2;
		t[1] = fp[1]:read():split("\t", 4);
		t[1][2], t[1][3] = tonumber(t[1][2]), tonumber(t[1][3]);
		return function()
			local ti = t[i];
			ti = fp[i]:read():split("\t", 4);
			ti[2], ti[3] = tonumber(ti[2]), tonumber(ti[3]);
			if t[1][1] == t[2][1] then i = (t[1][2] < t[2][2] and 1) or 2;
			else i = (t[1][1] < t[2][1] and 1) or 2 end
			return t[i];
		end
	end
end

function bed_union(fp)
	local c, b, e, n;
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
end

function bed_cov2(fp)
	local c, e, ob, oe; -- e is the rightmost pos;
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
end

-- main function

if #arg == 0 then
	print('Usage: bedtk.lua union|cov2 [in.bed]');
	os.exit(1);
end

local fp = io.xopen(arg[2]);
--[[
for t in bed_read2(fp) do
	print(t[1], t[2], t[3]);
end
os.exit(1);
]]--
if arg[1] == 'union' then bed_union(fp)
elseif arg[1] == 'cov2' then bed_cov2(fp)
else print('Unrecognized command.') end
