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


if #arg == 0 then
	print("Usage: vcf2bgl.lua <in.vcf>")
	os.exit(1)
end

local lookup = {}
for i = 0, 10000 do lookup[i] = string.format("%.4f", math.pow(10, -i/10)) end

local fp = io.xopen(arg[1])
for l in fp:lines() do
	if l:sub(1, 2) == '##' then -- meta lines; do nothing
	elseif l:sub(1, 1) == '#' then -- sample lines
		local t, s = l:split('\t'), {}
		for i = 10, #t do s[#s+1] = t[i]; s[#s+1] = t[i]; s[#s+1] = t[i] end
		print('marker', 'alleleA', 'alleleB', table.concat(s, '\t'))
	else -- data line
		local t = l:split('\t');
		if t[5] ~= '.' and t[5]:find(",") == nil and #t[5] == 1 and #t[4] == 1 then -- biallic SNP
			local x, z = -1, {};
			if t[9]:find('PL') then
				for i = 10, #t do
					local AA, Aa, aa = t[i]:match('(%d+),(%d+),(%d+)')
					AA = tonumber(AA); Aa = tonumber(Aa); aa = tonumber(aa);
					if AA ~= nil then
						z[#z+1] = lookup[AA]; z[#z+1] = lookup[Aa]; z[#z+1] = lookup[aa];
					else z[#z+1] = 1; z[#z+1] = 1; z[#z+1] = 1; end
				end
				print(t[1]..':'..t[2], t[4], t[5], table.concat(z, '\t'))
			elseif t[9]:find('GL') then
				print('Error: not implemented')
				os.exit(1)
			end
		end
	end
end
fp:close()
