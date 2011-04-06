#!/usr/bin/env luajit

function table.shuffle(a)
	for i = #a, 1, -1 do
		local j = math.random(i)
		a[j], a[i] = a[i], a[j]
	end
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
	print("Usage: permlines.lua <in.txt>")
	os.exit(1)
end

math.randomseed(os.time())
local fp, a = io.xopen(arg[1]), {}
for l in fp:lines() do a[#a+1] = l end
fp:close()
table.shuffle(a)
for _, l in pairs(a) do
	print(l)
end
