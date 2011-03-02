#!/usr/bin/env luajit

local bit = require("bit")

function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

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

local comp_table = {}
comp_table[string.byte('A')] = 'T';
comp_table[string.byte('T')] = 'A';
comp_table[string.byte('C')] = 'G';
comp_table[string.byte('G')] = 'C';

local fp = io.xopen(arg[1])
local last, out1, out2, k
for l in fp:lines() do
	if l:sub(1, 1) ~= '@' then
		local t = l:split('\t', 12)
		if t[10]:find('N') == nil then
			if last ~= t[1] then
				if out2 and k == 2 and out1[3] >= 25 then
					print(table.concat(out1, "\t"), table.concat(out2, "\t"))
					out2 = nil
				end
				k, last = 1, t[1]
			else k = k + 1 end
			if k == 1 then
				out1 = {t[1], t[10], t[11], tonumber(t[2]), t[3], t[4], t[5], t[6]}
			elseif k == 2 then
				out2 = {t[2], t[3], t[4], t[5], t[6]}
				if bit.band(out1[4], 16) ~= 0 then -- reverse complement
					local x, s = {}, out1[2]:reverse()
					for i = 1, #s do x[i] = comp_table[s:byte(i)] end
					out1[2] = table.concat(x)
				end
				-- compute the average base quality
				if out1[2]:find("[CGT]") == nil then out1[3] = 0 -- poly A
				else 
					local x, s = 0, out1[3];
					for i = 1, #s do x = x + s:byte(i) - 33 end
					x = math.floor(x / #s + .499)
					out1[3] = x
				end
			end -- do nothing if there are 3 or more hits
		end
	end
end
