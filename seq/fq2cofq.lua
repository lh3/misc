#!/usr/bin/env luajit

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

function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

if #arg < 1 then
	print('Usage: fq2cgfq.lua [-IF] [-b barcode] <in1.fq> [in2.fq]')
	os.exit(1)
end

local is_il13, is_flt, acgtn, bc_len, bc_hash = false, true, {'A', 'C', 'G', 'T', 'N'}, 0
for o, a in os.getopt(arg, 'IFb:') do
	if o == 'I' then is_il13 = true
	elseif o == 'b' then
		local t = a:split(',')
		bc_len, bc_hash = #t[1], {}
		for _, v in ipairs(t) do
			assert(#v == bc_len)
			for i = 1, bc_len do
				for _, j in ipairs(acgtn) do
					w = v:sub(1, i - 1) .. j .. v:sub(i + 1)
					bc_hash[w] = true
				end
			end
		end
	elseif o == 'F' then is_flt = true
	end
end

local fp = {io.xopen(arg[1]), nil}
if arg[2] then fp[2] = io.xopen(arg[2]) end
while true do
	local l = {{}, {}}
	for i = 1, 4 do l[1][i] = fp[1]:read() end
	if fp[2] then
		for i = 1, 4 do l[2][i] = fp[2]:read() end
	end
	if l[1][4] == nil or (fp[2] and l[2][4] == nil) then break end
	local bc = bc_len ~= 0 and l[1][2]:sub(1, bc_len) or nil
	if (is_flt == false or l[1][1]:find('^%S+:%d+:%d+:%d+:%d+:0#0') == nil) and (bc == nil or bc_hash[bc]) then
		if is_il13 then
			local q = {{}, {}}
			for i = 1, #l[1][4] do q[1][i] = string.char(l[1][4]:byte(i) - 31) end
			l[1][4] = table.concat(q[1])
			if fp[2] then
				for i = 1, #l[2][4] do q[2][i] = string.char(l[2][4]:byte(i) - 31) end
				l[2][4] = table.concat(q[2])
			end
		end
		print(l[1][1], bc) print(l[1][2]:sub(bc_len+1)) print('+') print(l[1][4]:sub(bc_len+1))
		if fp[2] then
			print(l[2][1]) print(l[2][2]) print('+') print(l[2][4])
		end
	end
end
fp[1]:close()
if fp[2] then fp[2]:close() end
