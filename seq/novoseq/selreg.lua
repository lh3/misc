#!/usr/bin/env luajit

bit = require('bit')
require('klib')

local min_s, min_l = 500, 1000

for o, a in os.getopt(arg, 's:l:') do
	if o == 's' then min_s = tonumber(a)
	elseif o == 'l' then min_l = tonumber(a)
	end
end

if #arg == 0 then
	print('Usage: luajit selreg.lua [-s mins] [-l minl] <in.sam>')
	os.exit(1)
end

function process_regs(regs, len, qname)
	table.sort(regs, function(a,b) return a[1] < b[1] end)
	local le, ret = 0, {}
	for i = 1, #regs do
		--print('*'..qname, len, regs[i][1], regs[i][2])
		local b, e = le, regs[i][1]
		if e - b >= min_l then ret[#ret+1] = {b, e} end
		le = le > regs[i][2] and le or regs[i][2]
	end
	if len - le >= min_l then ret[#ret+1] = {le, len} end
	-- print
	for i = 1, #ret do
		print(qname, ret[i][1], ret[i][2])
	end
end

local last, lastlen, regs = '', 0, {}
local fp = io.xopen(arg[1])
for l in fp:lines() do
	if l:sub(1, 1) == '@' then
	else
		local t = l:split('\t')
		if t[1] ~= last then
			if last ~= '' then
				process_regs(regs, lastlen, last)
			end
			last, lastlen, regs = '', 0, {}
		end
		local x, y, qbeg, qend = tonumber(t[4]) - 1, 0, 0, 0
		-- parse CIGAR
		if t[6] ~= '*' then
			for j, o in t[6]:gmatch('(%d+)(%S)') do
				if o == 'S' or o == 'H' then
					if y == 0 then
						qbeg = j -- this is the number for the first 'H'
					else
						qend = j -- this is for the last 'H'
					end
					y = y + j
				elseif o == 'M' then
					x, y = x + j, y + j
				elseif o == 'I' then
					y = y + j
				elseif o == 'D' then
					x = x + j
				end
			end
			qbeg, qend = tonumber(qbeg), tonumber(qend)
			-- parse score
			local score = l:match('\tAS:i:(%d+)')
			score = tonumber(score)
			if score >= min_s or (qbeg < min_s and qbeg < #t[10]) or (qend < min_s and qend < #t[10]) then -- do not filter hits at the ends
				-- parse strand
				if bit.band(tonumber(t[2]), 0x10) == 0x10 then
					qbeg, qend = qend, y - qbeg
				else qend = y - qend end
				regs[#regs+1] = {qbeg, qend}
			end
		elseif #t[10] >= min_l then -- unmapped query
			print(t[1], 0, #t[10])
		end
		-- set last
		last, lastlen = t[1], y
	end
end
if #regs > 0 then
	process_regs(regs, lastlen, last)
end
fp:close()
