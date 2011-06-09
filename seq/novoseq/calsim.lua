#!/usr/bin/env luajit

require('klib')

local score, fit, doflt, bed = { 1, 3, 5, 2 }, { 500, 0.95, 20000, 0.99 }, false, nil

for o, a in os.getopt(arg, 'B:a:b:q:r:l:L:fp') do
	if o == 'a' then score[1] = a
	elseif o == 'b' then score[2] = a
	elseif o == 'q' then score[3] = a
	elseif o == 'r' then score[4] = a
	elseif o == 'l' then fit[1] = a
	elseif o == 'L' then fit[3] = a
	elseif o == 'f' then doflt = true
	elseif o == 'B' then -- read BED
		local fp = io.xopen(a)
		bed = {}
		for l in fp:lines() do
			local t = l:split('\t')
			local s, b, e = t[1], tonumber(t[2]), tonumber(t[3])
			local s0, b0, e0 = s:match('^(%S+):(%d+)-(%d+)')
			if s0 ~= nil then
				s, b, e = s0, b0-1+b, b0-1+e
			end
			if bed[s] == nil then bed[s] = {} end
			table.insert(bed[s], {b, e, true})
		end
		fp:close()
	end
end

if #arg == 0 then
	print('Usage: calsim.lua [-f] [-B in.bed] [-a match] [-b mm] [-q gopen] [-r gext] [-l minl] [-L maxl] <in.sam>')
	os.exit(1)
end

function logit(a, b, x) return 1 / (1 + a * math.exp(-b * x)) end -- logit (well close to)

local a, b
b = math.log((1 / fit[2] - 1) / (1 / fit[4] - 1)) / (fit[3] - fit[1])
a = (1 / fit[2] - 1) / math.exp(-b * fit[1])

local last, fp = {'', 0, false}, io.xopen(arg[1])
for l in fp:lines() do
	local t, x = l:split('\t'), { 0, 0, 0, 0, 0, 0 } -- match, mismatch, insOpen, insGap, delOpen, delGap (Gap=Open+Ext)
	if t[1]:sub(1, 1) ~= '@' then -- skip header lines
		local qlen, tlen, h = #t[10], 0, 0 -- query length in the alignment, target length, length of H
		if t[3] ~= '*' then -- have CIGAR; compute identity, qlen and tlen
			local as = l:match('\tAS:i:(%d+)')
			as = tonumber(as)
			for j, o in t[6]:gmatch('(%d+)(%S)') do
				j = tonumber(j)
				if o == 'M' then x[1] = x[1] + j
				elseif o == 'I' then x[3], x[4] = x[3] + 1, x[4] + j
				elseif o == 'D' then x[5], x[6] = x[5] + 1, x[6] + j
				elseif o == 'H' then h = h + j
				end
			end
			qlen, tlen = x[1] + x[4], x[1] + x[6]
			x[2] = (qlen - score[3] * (x[3] + x[5]) - score[4] * x[6] - (score[4] + score[1]) * x[4] - as) / (score[1] + score[2])
			x[1] = x[1] - x[2]
			if bed ~= nil and bed[t[1]] ~= nil then
				b = t[6]:match('^(%d+)H') or 0
				e = t[6]:match('(%d+)H$') or 0
				b, e = tonumber(b), qlen + h - tonumber(e)
				for _, v in ipairs(bed[t[1]]) do
					if v[1] - b > fit[1] and e - v[2] > fit[1] and x[1]/qlen >= logit(a, b, qlen) then v[3] = false end
				end
			end
		end
		if t[1] ~= last[1] then
			if doflt and last[3] then print(last[1], "4\t*\t0\t0\t*\t*\t0\t0", string.rep('N', last[2]), '*') end
			last = {t[1], qlen + h, true}
		end
		if doflt then
			if x[1]/qlen >= logit(a, b, qlen) then print(l); last[3] = false; end
		elseif bed == nil then
			print(t[1], qlen + h, qlen, string.format("%.4f", x[1] / qlen), string.format("%.4f", logit(a, b, qlen)), t[3], t[4]-1, t[4]-1+tlen, t[5])
		end
	elseif doflt then print(l) end
end
if doflt and last[3] then print(last[1], "4\t*\t0\t0\t*\t*\t0\t0", string.rep('N', last[2]), '*') end
if bed ~= nil then
	for k, v in pairs(bed) do
		for _, r in ipairs(v) do
			if r[3] then print(k, r[1], r[2]) end
		end
	end
end
