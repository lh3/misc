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

function math.lgamma(z)
	local x;
	x = 0.1659470187408462e-06     / (z+7);
	x = x + 0.9934937113930748e-05 / (z+6);
	x = x - 0.1385710331296526     / (z+5);
	x = x + 12.50734324009056      / (z+4);
	x = x - 176.6150291498386      / (z+3);
	x = x + 771.3234287757674      / (z+2);
	x = x - 1259.139216722289      / (z+1);
	x = x + 676.5203681218835      / z;
	x = x + 0.9999999999995183;
	return math.log(x) - 5.58106146679532777 - z + (z-0.5) * math.log(z+6.5);
end

function math.igamma(s, z, complement)

	local function _kf_gammap(s, z)
		local sum, x = 1, 1;
		for k = 1, 100 do
			x = x * z / (s + k);
			sum = sum + x;
			if x / sum < 1e-14 then break end
		end
		return math.exp(s * math.log(z) - z - math.lgamma(s + 1.) + math.log(sum));
	end

	local function _kf_gammaq(s, z)
		local C, D, f, TINY;
		f = 1. + z - s; C = f; D = 0.; TINY = 1e-290;
		-- Modified Lentz's algorithm for computing continued fraction. See Numerical Recipes in C, 2nd edition, section 5.2
		for j = 1, 100 do
			local d;
			local a, b = j * (s - j), j*2 + 1 + z - s;
			D = b + a * D;
			if D < TINY then D = TINY end
			C = b + a / C;
			if C < TINY then C = TINY end
			D = 1. / D;
			d = C * D;
			f = f * d;
			if math.abs(d - 1) < 1e-14 then break end
		end
		return math.exp(s * math.log(z) - z - math.lgamma(s) - math.log(f));
	end

	if complement then
		return ((z <= 1 or z < s) and 1 - _kf_gammap(s, z)) or _kf_gammaq(s, z);
	else 
		return ((z <= 1 or z < s) and _kf_gammap(s, z)) or (1 - _kf_gammaq(s, z));
	end
end

matrix = {}
function matrix.chi2(a)
	if #a == 2 and #a[1] == 2 then -- 2x2 table
		local x, z
		x = (a[1][1] + a[1][2]) * (a[2][1] + a[2][2]) * (a[1][1] + a[2][1]) * (a[1][2] + a[2][2])
		if x == 0 then return 0, 1 end
		z = a[1][1] * a[2][2] - a[1][2] * a[2][1]
		z = (a[1][1] + a[1][2] + a[2][1] + a[2][2]) * z * z / x
		return z, math.igamma(.5, .5 * z, true)
	else -- generic table
		local rs, cs, n, m, N, z = {}, {}, #a, #a[1], 0, 0
		for i = 1, n do rs[i] = 0 end
		for j = 1, m do cs[j] = 0 end
		for i = 1, n do -- compute column sum and row sum
			for j = 1, m do cs[j], rs[i] = cs[j] + a[i][j], rs[i] + a[i][j] end
		end
		for i = 1, n do N = N + rs[i] end
		for i = 1, n do -- compute the chi^2 statistics
			for j = 1, m do
				local E = rs[i] * cs[j] / N;
				z = z + (a[i][j] - E) * (a[i][j] - E) / E
			end
		end
		return z, math.igamma(.5 * (n-1) * (m-1), .5 * z, true);
	end
end

-- main rountine

if #arg < 3 then
	print("Usage: vcf2chi2.lua <in.vcf> <group1.list> <group2.list>");
	os.exit(1)
end

local g = {};

-- read the list of groups
local fp = io.xopen(arg[2]);
for l in fp:lines() do local x = l:match("^(%S+)"); g[x] = 1 end -- FIXME: check duplicate
fp:close()
fp = io.xopen(arg[3]);
for l in fp:lines() do local x = l:match("^(%S+)"); g[x] = 2 end
fp:close()

-- process VCF
fp = io.xopen(arg[1])
local h = {{}, {}}
for l in fp:lines() do
	if l:sub(1, 2) == '##' then print(l) -- meta lines; do nothing
	elseif l:sub(1, 1) == '#' then -- sample lines
		local t = l:split('\t');
		for i = 10, #t do
			if g[t[i]] == 1 then table.insert(h[1], i)
			elseif g[t[i]] == 2 then table.insert(h[2], i) end
		end
		while #t > 8 do table.remove(t) end
		print(table.concat(t, "\t"))
	else -- data line
		local t = l:split('\t');
		if t[5] ~= '.' and t[5]:find(",") == nil then -- biallic
			local a = {{0, 0}, {0, 0}}
			for i = 1, 2 do
				for _, k in pairs(h[i]) do
					if t[k]:find("^0.0") then a[i][1] = a[i][1] + 2
					elseif t[k]:find("^1.1") then a[i][2] = a[i][2] + 2
					elseif t[k]:find("^0.1") or t[k]:find("^1.0") then
						a[i][1], a[i][2] = a[i][1] + 1, a[i][2] + 1
					end
				end
			end
			local chi2, p = matrix.chi2(a)
			while #t > 8 do table.remove(t) end
			--print(a[1][1], a[1][2], a[2][1], a[2][2], chi2, p);
			print(table.concat(t, "\t") .. ";PCHI2=" .. string.format("%.3g", p))
		end
	end
end
fp:close()
