#!/usr/bin/env luajit

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

if #arg < 2 then
	print('Usage: bgl2vcf.lua <in.phased> <in.gprobs>')
	os.exit(1)
end

local fpp = io.xopen(arg[1]);
local fpg = io.xopen(arg[2]);
for lg in fpg:lines() do
	local lp = fpp:read()
	local tp, tg, a = lp:split('%s'), lg:split('%s', 4), {}
	if tp[1] == 'I' then
		for i = 3, #tp, 2 do a[#a+1] = tp[i] end
		print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', table.concat(a, '\t'))
	else
		local chr, pos = tg[1]:match('(%S+):(%d+)$')
		a = {chr, pos, '.', tg[2], tg[3], 30, '.', '.', 'GT'}
		for i = 3, #tp, 2 do
			a[#a+1] = ((tp[i] == tg[2] and 0) or 1) .. '|' .. ((tp[i+1] == tg[2] and 0) or 1)
		end
		print(table.concat(a, '\t'))
	end
end
fpg:close(); fpp:close();

