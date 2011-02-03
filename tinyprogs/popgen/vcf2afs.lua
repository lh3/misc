#!/usr/bin/env luajit

--[[
  Author: lh3
  Description: From a VCF file, get the site allele frequency or the allele
               frequency spectrum in each population.
]]--

--------- routines from klua.lua ------------
function os.getopt(args, ostr)
	local arg, place = nil, 0;
	return function ()
		if place == 0 then -- update scanning pointer
			if #args == 0 or args[1]:sub(1, 1) ~= '-' then return nil end
			if #args[1] >= 2 then
				if args[1]:sub(2, 2) == '-' then -- found "--"
					table.remove(args, 1);
					return nil;
				end
			end
			place = 2
		end
		local optopt = args[1]:sub(place, place);
		place = place + 1;
		local oli = ostr:find(optopt);
		if optopt == ':' or oli == nil then -- unknown option
			if optopt == '-' then return nil end
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
			return '?';
		end
		oli = oli + 1;
		if ostr:sub(oli, oli) ~= ':' then -- do not need argument
			arg = nil;
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
		else -- need an argument
			if place <= #args[1] then  -- no white space
				arg = args[1]:sub(place);
			else
				table.remove(args, 1);
				if #args == 0 then -- an option requiring argument is the last one
					place = 0;
					if ostr:sub(1, 1) == ':' then return ':' end
					return '?';
				else arg = args[1] end
			end
			table.remove(args, 1);
			place = 0;
		end
		return optopt, arg;
	end
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

----------- klua routines end here ------------

-- parse the command line
local site_only = false; -- print site allele frequency or not
for c in os.getopt(arg, 's') do
	if c == 's' then site_only = true end
end
if #arg == 0 then
	print("Usage: luajit vcf2afs.lua [-s] <samples.txt> <in.vcf>")
	os.exit(1)
end

-- read the sample-population pairs
local pop, sample = {}, {}
local fp = io.open(arg[1]);
for l in fp:lines() do
	local s, p = l:match("^(%S+)%s+(%S+)"); -- sample, population pair
	sample[s] = p; -- FIXME: check duplications
	if pop[p] then table.insert(pop[p], s);
	else pop[p] = {s} end
end
fp:close();

-- parse VCF
fp = (#arg >= 2 and io.open(arg[2])) or io.stdin;
local col, cnt = {}, {};
for k in pairs(pop) do
	col[k], cnt[k] = {}, {[0]=0};
end
for l in fp:lines() do
	if l:sub(1, 2) == '##' then -- meta lines; do nothing
	elseif l:sub(1, 1) == '#' then -- the sample line
		local t = l:split('\t');
		for i = 10, #t do
			local k = sample[t[i]];
			if k ~= nil then
				table.insert(col[k], i);
				table.insert(cnt[k], 0);
				table.insert(cnt[k], 0);
			end
		end
	else -- data lines
		local t = l:split('\t');
		if t[5] ~= '.' and t[5]:find(",") == nil then -- biallic
			if site_only == true then io.write(t[1], '\t', t[2], '\t', t[4], '\t', t[5]) end
			for k, v in pairs(col) do
				local ac, an = 0, 0;
				for i = 1, #v do
					local a1, a2 = t[v[i]]:match("^(%d).(%d)");
					if a1 ~= nil then ac, an = ac + a1 + a2, an + 2 end
				end
				if site_only == true then io.write('\t', k, ':', an, ':', ac) end
				if an == #cnt[k] then cnt[k][ac] = cnt[k][ac] + 1 end
			end
			if site_only == true then io.write('\n') end
		end
	end
end
fp:close();

-- print
if site_only == false then
	for k, v in pairs(cnt) do
		io.write(k .. "\t" .. #v);
		for i = 0, #v do io.write("\t" .. v[i]) end
		io.write('\n');
	end
end
