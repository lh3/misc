#!/usr/bin/env luajit

require('klib')

function logit(a, b, x) return 1 / (1 + a * math.exp(-b * x)) end -- logit (well close to)

local min_drop, fit = 100, { 500, 0.95, 20000, 0.99 }
local a, b
b = math.log((1 / fit[2] - 1) / (1 / fit[4] - 1)) / (fit[3] - fit[1])
a = (1 / fit[2] - 1) / math.exp(-b * fit[1])

if #arg == 0 then
	print('Usage: psldrop.psl <all-vs-all.psl>')
	os.exit(1)
end

local fp = io.xopen(arg[1])
for l in fp:lines() do
	if l:find('^%d') ~= nil then
		local t = l:split('\t')
		if #t > 18 and t[10]:find(':') ~= nil and t[14]:find(':') ~= nil then
			for i = 1, 8 do t[i] = tonumber(t[i]) end
			for i = 11, 13 do t[i] = tonumber(t[i]) end
			for i = 15, 18 do t[i] = tonumber(t[i]) end
			local r = { t[13] - t[12], t[17] - t[15] }
			local mm = t[2] + t[6] + t[8]
			local err = mm / (t[1] + mm)
			if t[10] ~= t[14] and (err < 1 - logit(a, b, t[1] + mm) or err < 1 - fit[4]) and (t[11] - r[1] < min_drop or t[15] - r[2] < min_drop) then
			--if t[10] ~= t[14] and err < 1 - fit[4] and (t[11] - r[1] < min_drop or t[15] - r[2] < min_drop) then
				if t[11] < t[15] then print(t[10], t[14], t[1], mm)
				elseif t[11] > t[15] then print(t[14], t[10], t[1], mm)
				else
					if t[10] > t[14] then print(t[10], t[14], t[1], mm)
					else print(t[14], t[10], t[1], mm) end
				end
			end
		end
	end
end
fp:close()
