#!/usr/bin/env luajit

require('klib')

local score, fit, doflt, printsam = { 1, 3, 5, 2 }, { 500, 0.95, 20000, 0.99 }, false, false

for o, a in os.getopt(arg, 'a:b:q:r:l:L:fp') do
	if o == 'a' then score[1] = a
	elseif o == 'b' then score[2] = a
	elseif o == 'q' then score[3] = a
	elseif o == 'r' then score[4] = a
	elseif o == 'l' then fit[1] = a
	elseif o == 'L' then fit[3] = a
	elseif o == 'f' then doflt = true
	elseif o == 'p' then printsam = true
	end
end

if #arg == 0 then
	print('Usage: fltreg.lua [-fp] [-a match] [-b mm] [-q gopen] [-r gext] [-l minl] [-L maxl] <in.sam>')
	os.exit(1)
end

function logit(a, b, x) return 1 / (1 + a * math.exp(-b * x)) end

local fp = io.xopen(arg[1])
local a, b

b = math.log((1 / fit[2] - 1) / (1 / fit[4] - 1)) / (fit[3] - fit[1])
a = (1 / fit[2] - 1) / math.exp(-b * fit[1])

local all, drop = {}, {}
for l in fp:lines() do
	local t, x = l:split('\t'), { 0, 0, 0, 0, 0, 0 } -- match, mismatch, insOpen, insGap, delOpen, delGap (Gap=Open+Ext)
	if t[1]:sub(1, 1) ~= '@' then all[t[1]] = 1 end
	if t[1]:sub(1, 1) ~= '@' and t[3] ~= '*' then
		local as, h = l:match('\tAS:i:(%d+)'), 0
		as = tonumber(as)
		for j, o in t[6]:gmatch('(%d+)(%S)') do
			j = tonumber(j)
			if o == 'M' then x[1] = x[1] + j
			elseif o == 'I' then x[3], x[4] = x[3] + 1, x[4] + j
			elseif o == 'D' then x[5], x[6] = x[5] + 1, x[6] + j
			elseif o == 'H' then h = h + j
			end
		end
		local qlen, tlen = x[1] + x[4], x[1] + x[6]
		local tmp = qlen - score[3] * (x[3] + x[5]) - score[4] * x[6] - (score[4] + score[1]) * x[4]
		x[2] = (tmp - as) / (score[1] + score[2])
		x[1] = x[1] - x[2]
		if (x[1] / qlen > fit[4] or x[1] / qlen > logit(a, b, qlen)) and h < fit[1] then
			drop[t[1]] = 1
		end
		if doflt == false then
			if printsam then
				if x[1]/qlen >= logit(a, b, qlen) then print(l) end
			else
				print(t[1], qlen + h, qlen, string.format("%.4f", x[1] / qlen), string.format("%.4f", logit(a, b, qlen)), t[3], t[4]-1, t[4]-1+tlen, t[5])
			end
		end
	end
end

if doflt == true then
	for k, _ in pairs(all) do
		if drop[k] == nil then print(k) end
	end
end
