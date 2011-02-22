#!/usr/bin/env luajit

-- Description: string split
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

-- Description: intelligent file open
function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

if #arg < 3 then
	print("Usage: bc2rg.lua <barcode.list> <in.sam> <libname>");
	os.exit(1);
end

-- read barcode.list
local hash = {};
local hdr = {};
local nuc = {'A', 'C', 'G', 'T'};
local fp = io.xopen(arg[1]);
for line in fp:lines() do
	local bc, sam = line:match("^(%S+)%s+(%S+)");
	local rg = sam..'-'..arg[3]..'-'..bc;
	hash[bc] = rg;
	table.insert(hdr, "@RG\tID:"..rg.."\tSM:"..sam.."\tLB:"..rg);
	for i = 1, #bc do
		for j = 1, 4 do
			if nuc[j] ~= bc:sub(i, i) then
				local b = bc:sub(1, i-1) .. nuc[j] .. bc:sub(i+1);
				hash[b] = rg;
			end
		end
	end
end
table.insert(hdr, "@RG\tID:N/A\tSM:N/A\tLB:N/A");
fp:close();

-- process the SAM file
fp = io.xopen(arg[2]);
local first = true;
for l in fp:lines() do
	if l:sub(1,1) == '@' then
		print(l);
	else
		if first then
			print(table.concat(hdr, "\n"));
			first = false;
		end
		local bc = l:match("BC:Z:(%S+)"):upper();
		if bc == nil or hash[bc] == nil then
			print(l, "RG:Z:N/A");
		else
			print(l, "RG:Z:" .. hash[bc]);
		end
	end
end
fp:close();
