#!/usr/bin/env luajit

local a = {}
for l in io.lines() do
	a[#a+1] = tonumber(l);
end
table.sort(a);
for i = 1, #a do
	local x = (a[i] == 0 and 1e-20) or a[i];
	print(-math.log10(i/#a), -math.log10(a[i]))
end
